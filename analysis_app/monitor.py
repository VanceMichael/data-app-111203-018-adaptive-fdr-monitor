"""持续筛查：在线 FDR 监控器的主编排。

处理顺序（确定性）：每次摄入按 (visible_at, machine, family, tick,
observation_id) 稳定排序后顺序处理。已花费预算的位置只增不改；
visible_at 低于已见水位线的迟到记录被隔离（quarantine），不会插入
历史位置；重复 observation_id 直接跳过，不重复消费预算。

层级错误预算：

1. 全局预算 ``global_alpha`` 按固定 family 权重切分，
   family 包络 alpha_f = global_alpha * w_f，权重在监控器初始化时
   冻结（状态恢复后不可变）；
2. 每个 family 内部再按固定比例切成两个流：LORD3 流
   （满足独立/PRDS 与超均匀假设时使用）与在线 Bonferroni 流
   （任意依赖下保守降级，按并集界控制 FWER）；
3. family 内任一发现必须通过对应流的阈值；层级并集界给出
   全局 FDR <= global_alpha（详见 docs 中的规则说明）。

参考分布由 ``references.ReferenceBook`` 按检修 epoch 独立维护，
检修只触发新 epoch 的参考重建，不触碰 family/全局错误账本。
"""
import json
import math
from dataclasses import dataclass, field, asdict

from .fdr import OnlineFDRController, Stream
from .references import ReferenceBook
from .stats import lag1_autocorrelation

STATE_VERSION = 2

DISPOSITION_TESTED = 'tested'
DISPOSITION_COLD_START = 'cold_start'
DISPOSITION_REFERENCE_MEMBER = 'reference_member'
DISPOSITION_MAINTENANCE = 'maintenance'
DISPOSITION_DUPLICATE = 'duplicate'
DISPOSITION_LATE = 'late_arrival'
DISPOSITION_QUARANTINED = 'quarantined'

RULE_ASSUMPTIONS = {
    'primary_rule': 'lord3',
    'gamma_sequence': 'gamma_t = (6/pi^2) / t^2，t>=1，可求和，sum gamma_t = 1',
    'lord3_initial_wealth': 'W0 = envelope_alpha / 2',
    'lord3_conditions': [
        '同一流内各次检验的 p 值相互独立，或满足 PRDS 正向依赖',
        '每个真零假设下 p 值为超均匀（super-uniform）分布',
        '参考分布从指定历史稳定窗口冻结估计，窗口外观测不回写参考',
    ],
    'fallback_rule': 'online_bonferroni',
    'fallback_conditions': [
        '仅要求真零假设下 p 值超均匀，不限制序列依赖结构',
        '阈值 beta_t = envelope_alpha * gamma_t，并集界给出 FWER <= envelope_alpha',
    ],
    'fallback_triggers': [
        '检验残差（按 epoch 去均值、滚动窗口）的滞后 1 阶自相关绝对值连续 3 次超过 ar_cutoff 触发降级；连续 3 次回到阈值内可恢复',
        '配置 force_conservative=true 时所有检验直接进入 Bonferroni 流',
    ],
    'global_bound': ('全局 alpha 按固定 family 权重切分，family 内按固定比例切流；'
                     '由层级并集界，全局 FDR <= global_alpha'),
}


@dataclass(frozen=True)
class MonitorConfig:
    global_alpha: float = 0.05
    reference_size: int = 4
    lord_share: float = 0.8
    ar_cutoff: float = 0.5
    ar_min_samples: int = 30
    family_weights: dict = field(default_factory=dict)
    force_conservative: bool = False

    def validate(self):
        if not 0.0 < self.global_alpha < 1.0:
            raise ValueError('global_alpha 必须位于 (0, 1)')
        if self.reference_size < 2:
            raise ValueError('reference_size 至少为 2')
        if not 0.0 <= self.lord_share <= 1.0:
            raise ValueError('lord_share 必须位于 [0, 1]')
        if not 0.0 < self.ar_cutoff < 1.0:
            raise ValueError('ar_cutoff 必须位于 (0, 1)')
        if self.ar_min_samples < 3:
            raise ValueError('ar_min_samples 至少为 3')
        if self.family_weights:
            weights = list(self.family_weights.values())
            if any(weight < 0 for weight in weights):
                raise ValueError('family 权重不能为负')
            total = math.fsum(weights)
            if not 0.999 <= total <= 1.001:
                raise ValueError('family 权重之和必须为 1')

    def to_dict(self):
        return asdict(self)


@dataclass
class FamilyRuntime:
    weight: float
    controllers: dict

    def to_dict(self):
        return {'weight': self.weight,
                'controllers': {key: controller.status().to_dict()
                                for key, controller in self.controllers.items()}}


class Monitor:
    def __init__(self, config=None):
        self.config = config or MonitorConfig()
        self.config.validate()
        self.reference_book = ReferenceBook(self.config.reference_size)
        self.families = {}
        self.processed_ids = set()
        self.watermark = None
        self.global_position = 0
        self.budget_ledger = []
        self.discoveries = []
        self.dispositions = []
        self.pending = []
        self._diagnostics = {}
        self._families_frozen = False
        self._equal_weights = {}

    # ------------------------------------------------------------------
    # 层级预算初始化
    # ------------------------------------------------------------------
    def _ensure_family(self, family):
        runtime = self.families.get(family)
        if runtime is not None:
            return runtime
        if self.config.family_weights:
            if family not in self.config.family_weights:
                return None
            weight = self.config.family_weights[family]
        else:
            # 未显式配置权重时，family 集合必须在首次摄入前通过
            # prime_families 冻结为等分权重；运行中出现的新 family
            # 不能事后切分既有预算
            if not self._families_frozen or family not in self._equal_weights:
                return None
            weight = self._equal_weights[family]
        envelope = self.config.global_alpha * weight
        runtime = FamilyRuntime(
            weight=weight,
            controllers={
                Stream.LORD3.value: OnlineFDRController(
                    envelope_alpha=envelope * self.config.lord_share,
                    stream=Stream.LORD3),
                Stream.ONLINE_BONFERRONI.value: OnlineFDRController(
                    envelope_alpha=envelope * (1.0 - self.config.lord_share),
                    stream=Stream.ONLINE_BONFERRONI),
            })
        self.families[family] = runtime
        return runtime

    def prime_families(self, family_names):
        """首次摄入前冻结 family 集合；显式权重按配置，否则等分。"""
        if self._families_frozen:
            return
        names = sorted(set(family_names))
        if self.config.family_weights:
            # 未登记 family 不进入权重表，由逐记录路径隔离为 quarantined
            self._equal_weights = dict(self.config.family_weights)
        else:
            if not names:
                self._equal_weights = {}
            else:
                share = 1.0 / len(names)
                self._equal_weights = {name: share for name in names}
        self._families_frozen = True

    # ------------------------------------------------------------------
    # 依赖诊断（仅作为降级触发，不用于宣称 FDR 成立）
    # ------------------------------------------------------------------
    def _diagnostic_state(self, machine, epoch):
        """依赖诊断按 (设备, epoch) 独立累计；检修开启新 epoch 后重新开始，
        水平漂移不会跨 epoch 污染自相关估计。"""
        key = (machine, epoch)
        state = self._diagnostics.get(key)
        if state is None:
            state = {'values': [], 'degraded': False, 'last_ar': None,
                     'consecutive_breaches': 0, 'consecutive_clear': 0}
            self._diagnostics[key] = state
        return state

    def _route_stream(self, machine, epoch, residual, reference):
        state = self._diagnostic_state(machine, epoch)
        # 用相对冻结参考均值中心化后的残差进入滚动诊断窗口；窗口长度
        # 固定为 ar_min_samples，旧观测淡出，历史突发/漂移不会永久
        # 污染自相关估计。估计前再对窗口去均值，使水平均值跳变不被
        # 误判为序列相关，只有持续的滞后结构才触发降级。
        centered = residual - reference.mean
        state['values'].append(centered)
        if len(state['values']) > self.config.ar_min_samples:
            del state['values'][:-self.config.ar_min_samples]
        if self.config.force_conservative:
            return Stream.ONLINE_BONFERRONI.value, state
        if len(state['values']) >= self.config.ar_min_samples:
            window_mean = math.fsum(state['values']) / len(state['values'])
            ar1 = lag1_autocorrelation([value - window_mean
                                        for value in state['values']])
            state['last_ar'] = ar1
            if ar1 is not None and abs(ar1) > self.config.ar_cutoff:
                state['consecutive_breaches'] += 1
                state['consecutive_clear'] = 0
            else:
                state['consecutive_clear'] += 1
                state['consecutive_breaches'] = 0
            # 连续 3 次越界黏性降级；降级后连续 3 次回到阈值内则恢复，
            # 避免一次性结构变化导致长期浪费 LORD 预算
            if not state['degraded'] and state['consecutive_breaches'] >= 3:
                state['degraded'] = True
            elif state['degraded'] and state['consecutive_clear'] >= 3:
                state['degraded'] = False
        return (Stream.ONLINE_BONFERRONI.value if state['degraded']
                else Stream.LORD3.value), state

    # ------------------------------------------------------------------
    # 摄入
    # ------------------------------------------------------------------
    @staticmethod
    def _order_key(row):
        return (row['visible_at'], row['machine'], row['family'],
                row['tick'], row['observation_id'])

    def ingest(self, records, input_fingerprint=None):
        """摄入一个批次的记录（可重复调用），返回本批次处置清单。

        记录先按当前水位线预筛：低于水位线的迟到记录立即隔离，不进入
        待处理缓冲；其余记录进入按 (visible_at, machine, family, tick,
        observation_id) 排序的缓冲并与历史缓冲合并排空，使跨批次仍保持
        全局检验顺序。迟到记录永远不会插入已经花费预算的历史位置。
        """
        batch_dispositions = []
        fresh = []
        for row in records:
            if row['observation_id'] in self.processed_ids:
                # 重复输入仍记录审计处置，但绝不进入待处理缓冲、不消费预算
                duplicate = {
                    'observation_id': row['observation_id'],
                    'machine': row['machine'], 'family': row['family'],
                    'tick': row['tick'], 'visible_at': row['visible_at'],
                    'status': DISPOSITION_DUPLICATE,
                    'reason': 'observation_id 已处理，跳过且不重复消费预算',
                    'input_fingerprint': input_fingerprint}
                self.dispositions.append(duplicate)
                batch_dispositions.append(duplicate)
            else:
                fresh.append(row)
        if fresh:
            self.prime_families([row['family'] for row in fresh])
        for row in fresh:
            # 入缓冲前按当前水位线预筛：已经晚于水位线的记录立即隔离，
            # 避免缓冲重排把它推到高水位线之后误判，也保证缓冲内记录
            # 全部仍有资格参与顺序检验
            if self.watermark is not None and row['visible_at'] < self.watermark:
                self.processed_ids.add(row['observation_id'])
                late = {
                    'observation_id': row['observation_id'],
                    'machine': row['machine'], 'family': row['family'],
                    'tick': row['tick'], 'visible_at': row['visible_at'],
                    'status': DISPOSITION_LATE,
                    'reason': (f"visible_at={row['visible_at']} 低于水位线 "
                               f"{self.watermark}，隔离且不插入历史预算位置"),
                    'input_fingerprint': input_fingerprint}
                self.dispositions.append(late)
                batch_dispositions.append(late)
            else:
                self.pending.append(row)
        self.pending.sort(key=self._order_key)
        while self.pending:
            row = self.pending.pop(0)
            disposition = self._ingest_one(row)
            disposition['input_fingerprint'] = input_fingerprint
            self.dispositions.append(disposition)
            batch_dispositions.append(disposition)
        return batch_dispositions

    def _ingest_one(self, row):
        observation_id = row['observation_id']
        base = {'observation_id': observation_id, 'machine': row['machine'],
                'family': row['family'], 'tick': row['tick'],
                'visible_at': row['visible_at']}
        if observation_id in self.processed_ids:
            return {**base, 'status': DISPOSITION_DUPLICATE,
                    'reason': 'observation_id 已处理，跳过且不重复消费预算'}
        if self.watermark is not None and row['visible_at'] < self.watermark:
            self.processed_ids.add(observation_id)
            return {**base, 'status': DISPOSITION_LATE,
                    'reason': f"visible_at={row['visible_at']} 低于水位线 {self.watermark}，隔离且不插入历史预算位置"}
        runtime = self._ensure_family(row['family'])
        if runtime is None:
            self.processed_ids.add(observation_id)
            return {**base, 'status': DISPOSITION_QUARANTINED,
                    'reason': 'family 未在固定权重配置中登记，不能事后切分预算'}
        self.watermark = row['visible_at'] if self.watermark is None else max(
            self.watermark, row['visible_at'])
        epoch, reference_status = self.reference_book.observe(row)
        if reference_status == 'maintenance':
            self.processed_ids.add(observation_id)
            closed_epoch = self.reference_book.epoch_state(row['machine'], epoch)
            reclassified = []
            if closed_epoch is not None and closed_epoch.closed:
                for member_id in closed_epoch.unfrozen_member_ids:
                    original = next((item for item in self.dispositions
                                     if item['observation_id'] == member_id), None)
                    if original is not None:
                        original['status'] = DISPOSITION_COLD_START
                        original['reason'] = (
                            f'epoch {epoch} 在参考窗口收齐前因检修关闭，'
                            '参考不足，冷启动不检验，不消费预算')
                        reclassified.append(member_id)
            return {**base, 'status': DISPOSITION_MAINTENANCE, 'epoch': epoch,
                    'cold_start_reclassified': reclassified,
                    'reason': '检修记录：关闭当前 epoch，下一 tick 起新建参考，错误账本保持不变'}
        if reference_status == 'reference_member':
            self.processed_ids.add(observation_id)
            return {**base, 'status': DISPOSITION_REFERENCE_MEMBER, 'epoch': epoch,
                    'reason': '进入本 epoch 的历史稳定参考窗口，不做检验'}
        if reference_status == 'late_closed':
            self.processed_ids.add(observation_id)
            return {**base, 'status': DISPOSITION_COLD_START, 'epoch': epoch,
                    'reason': '迟到记录落入未冻结即关闭的 epoch，参考不足，冷启动不检验'}
        reference = self.reference_book.reference_for(row['machine'], epoch)
        if reference is None:
            self.processed_ids.add(observation_id)
            return {**base, 'status': DISPOSITION_COLD_START, 'epoch': epoch,
                    'reason': f'epoch {epoch} 稳定窗口不足 {self.config.reference_size} 条，冷启动不检验'}
        return self._test_record(row, epoch, reference, runtime, base)

    def _test_record(self, row, epoch, reference, runtime, base):
        p_value = reference.p_value(row['residual'])
        stream_name, diagnostic = self._route_stream(
            row['machine'], epoch, row['residual'], reference)
        controller = runtime.controllers[stream_name]
        local_position = controller.position + 1
        self.global_position += 1
        entry = controller.submit(p_value)
        ledger_entry = {
            'global_position': self.global_position,
            'family': row['family'], 'family_weight': runtime.weight,
            'stream': stream_name,
            'family_stream_position': local_position,
            'observation_id': row['observation_id'],
            'machine': row['machine'], 'tick': row['tick'],
            'visible_at': row['visible_at'], 'epoch': epoch,
            'envelope_alpha': entry['envelope_alpha'],
            'alpha': entry['alpha'], 'p_value': p_value,
            'rejected': entry['rejected'],
            'null_hypothesis': ('残差服从该设备 epoch {epoch} 冻结参考的标准化分布，'
                                'Student-t 双侧检验，自由度 {df}').format(
                epoch=epoch, df=reference.size - 1),
            'validity_conditions': {
                'reference_epoch': epoch,
                'reference_start_tick': reference.start_tick,
                'reference_end_tick': reference.end_tick,
                'reference_size': reference.size,
                'reference_observation_ids': list(reference.observation_ids),
                'frozen_mean': reference.mean, 'frozen_stddev': reference.stddev,
                'stream_assumptions': (RULE_ASSUMPTIONS['lord3_conditions']
                                       if stream_name == Stream.LORD3.value
                                       else RULE_ASSUMPTIONS['fallback_conditions']),
                'ar1_lag1': diagnostic['last_ar'],
                'degraded': diagnostic['degraded'],
                'forced_conservative': self.config.force_conservative,
            }}
        self.budget_ledger.append(ledger_entry)
        self.processed_ids.add(row['observation_id'])
        disposition = {**base, 'status': DISPOSITION_TESTED, 'epoch': epoch,
                       'stream': stream_name, 'p_value': p_value,
                       'alpha': entry['alpha'], 'rejected': entry['rejected'],
                       'global_position': self.global_position}
        if entry['rejected']:
            self.discoveries.append({
                'global_position': self.global_position,
                'observation_id': row['observation_id'],
                'machine': row['machine'], 'family': row['family'],
                'tick': row['tick'], 'visible_at': row['visible_at'],
                'epoch': epoch, 'stream': stream_name,
                'p_value': p_value, 'alpha': entry['alpha'],
                'reference': reference.to_dict()})
        return disposition

    def finalize(self):
        """摄入结束后的收尾：仍打开但未冻结 epoch 的窗口成员按冷启动处理。

        这些观测的参考窗口始终未收齐，不能检验、不消费预算。收尾动作
        幂等；状态恢复后再次调用不会重复标记。
        """
        for epoch in self.reference_book.iter_epochs():
            if epoch.frozen or epoch.closed:
                continue
            for member_id in epoch.window_ids:
                original = next((item for item in self.dispositions
                                 if item['observation_id'] == member_id), None)
                if original is not None and original['status'] == DISPOSITION_REFERENCE_MEMBER:
                    original['status'] = DISPOSITION_COLD_START
                    original['reason'] = (
                        f'epoch {epoch.epoch} 在数据结束时参考窗口仍不足 '
                        f'{self.config.reference_size} 条，冷启动不检验，不消费预算')
            epoch.finalized = True
            epoch.unfrozen_member_ids = list(epoch.window_ids)

    # ------------------------------------------------------------------
    # 状态快照与恢复
    # ------------------------------------------------------------------
    def state_snapshot(self):
        return {
            'state_version': STATE_VERSION,
            'config': self.config.to_dict(),
            'watermark': self.watermark,
            'global_position': self.global_position,
            'families_frozen': self._families_frozen,
            'equal_weights': dict(self._equal_weights),
            'processed_observation_ids': sorted(self.processed_ids),
            'pending': self.pending,
            'references': self.reference_book.snapshot(),
            'families': {name: runtime.to_dict()
                         for name, runtime in sorted(self.families.items())},
            'diagnostics': {f'{machine}|{epoch}': {
                'values': state['values'], 'degraded': state['degraded'],
                'last_ar': state['last_ar'],
                'consecutive_breaches': state['consecutive_breaches'],
                'consecutive_clear': state['consecutive_clear']}
                for (machine, epoch), state in sorted(self._diagnostics.items())},
            'budget_ledger': self.budget_ledger,
            'discoveries': self.discoveries,
            'dispositions': self.dispositions,
        }

    def save_state(self, path):
        from .cli import save_json
        save_json(path, self.state_snapshot())

    @classmethod
    def restore(cls, payload):
        if payload.get('state_version') != STATE_VERSION:
            raise ValueError('状态版本不兼容，无法恢复')
        config = MonitorConfig(**payload['config'])
        config.validate()
        monitor = cls(config=config)
        monitor.watermark = payload['watermark']
        monitor.global_position = payload['global_position']
        monitor._families_frozen = payload.get('families_frozen', False)
        monitor._equal_weights = dict(payload.get('equal_weights', {}))
        monitor.processed_ids = set(payload['processed_observation_ids'])
        monitor.pending = list(payload.get('pending', []))
        monitor.reference_book = ReferenceBook.restore(payload['references'])
        monitor.families = {}
        for name, item in payload['families'].items():
            controllers = {key: OnlineFDRController.restore(snapshot)
                           for key, snapshot in item['controllers'].items()}
            monitor.families[name] = FamilyRuntime(weight=item['weight'],
                                                   controllers=controllers)
        for composite, state in payload['diagnostics'].items():
            machine, epoch_text = composite.split('|')
            monitor._diagnostics[(machine, int(epoch_text))] = {
                'values': list(state['values']),
                'degraded': state['degraded'],
                'last_ar': state['last_ar'],
                'consecutive_breaches': state.get('consecutive_breaches', 0),
                'consecutive_clear': state.get('consecutive_clear', 0)}
        monitor.budget_ledger = list(payload['budget_ledger'])
        monitor.discoveries = list(payload['discoveries'])
        monitor.dispositions = list(payload['dispositions'])
        return monitor

    # ------------------------------------------------------------------
    # 报告
    # ------------------------------------------------------------------
    def current_state(self):
        stream_counts = {}
        for name, runtime in sorted(self.families.items()):
            for stream_name, controller in runtime.controllers.items():
                status = controller.status()
                stream_counts[f'{name}:{stream_name}'] = {
                    'position': status.position,
                    'reject_count': status.reject_count,
                    'wealth': status.wealth,
                    'envelope_alpha': status.envelope_alpha}
        return {
            'watermark_visible_at': self.watermark,
            'global_position': self.global_position,
            'processed_observations': len(self.processed_ids),
            'pending_observations': len(self.pending),
            'families': {name: {
                'weight': runtime.weight,
                'streams': {key: controller.status().to_dict()
                            for key, controller in runtime.controllers.items()}}
                for name, runtime in sorted(self.families.items())},
            'stream_counters': stream_counts,
            'references': self.reference_book.snapshot(),
            'degraded_machine_epochs': sorted(
                f'{machine}|epoch-{epoch}' for (machine, epoch), state
                in self._diagnostics.items() if state['degraded']),
            'disposition_counts': self._count_dispositions(),
        }

    def _count_dispositions(self):
        counts = {}
        for item in self.dispositions:
            counts[item['status']] = counts.get(item['status'], 0) + 1
        return dict(sorted(counts.items()))

    def report(self, input_fingerprint=None, rows_total=0):
        return {
            'product': 'adaptive-fdr-monitor',
            'schema_version': STATE_VERSION,
            'config': self.config.to_dict(),
            'input_fingerprint': input_fingerprint,
            'input_rows': rows_total,
            'rule_assumptions': RULE_ASSUMPTIONS,
            'budget_ledger': self.budget_ledger,
            'discoveries': self.discoveries,
            'dispositions': self.dispositions,
            'current_state': self.current_state(),
        }


def load_state(path):
    with open(path, encoding='utf-8') as handle:
        return json.load(handle)
