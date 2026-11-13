"""按检修边界切分的设备参考分布（参考窗口）管理。

每台设备的时间轴按 maintenance=1 的观测切成若干 epoch：

* epoch 0 从首个观测开始；
* 每遇到一条 maintenance=1 的观测，当前 epoch 在该 tick 关闭，
  下一个 tick 开启新的 epoch；
* 每个 epoch 的参考分布只能由该 epoch 自身稳定窗口（开头若干条
  maintenance=0 观测）建立，检修前与检修后的窗口永远不会混用；
* 参考窗口不足 ``reference_size`` 条，或窗口内残差无波动（标准差为 0）
  时，该 epoch 处于冷启动（cold_start），不产生 p 值、不消费错误预算；
* 参考在窗口收齐那一刻冻结，之后不可变；检修只会新建 epoch，
  不会修改或删除任何已冻结参考，也不触碰全局错误账本。

检验使用冻结的样本均值与样本标准差，对窗口外观测做 Student-t 检验
（自由度 n-1，含参考估计不确定性），输出双侧 p 值。
"""
import math
from dataclasses import dataclass, field

from .stats import t_two_sided_p


@dataclass
class FrozenReference:
    epoch: int
    machine: str
    start_tick: int
    end_tick: int
    size: int
    mean: float
    stddev: float
    observation_ids: list

    def p_value(self, residual):
        statistic = (residual - self.mean) / self.stddev
        return t_two_sided_p(statistic, self.size - 1)

    def to_dict(self):
        return {'epoch': self.epoch, 'machine': self.machine,
                'start_tick': self.start_tick, 'end_tick': self.end_tick,
                'size': self.size, 'mean': self.mean,
                'stddev': self.stddev, 'observation_ids': list(self.observation_ids)}

    @classmethod
    def from_dict(cls, payload):
        return cls(epoch=payload['epoch'], machine=payload['machine'],
                   start_tick=payload['start_tick'], end_tick=payload['end_tick'],
                   size=payload['size'], mean=payload['mean'],
                   stddev=payload['stddev'],
                   observation_ids=list(payload['observation_ids']))


@dataclass
class EpochState:
    epoch: int
    machine: str
    start_tick: int
    window_values: list = field(default_factory=list)
    window_ids: list = field(default_factory=list)
    window_ticks: list = field(default_factory=list)
    reference: object = None
    closed: bool = False
    finalized: bool = False
    unfrozen_member_ids: list = field(default_factory=list)

    @property
    def frozen(self):
        return self.reference is not None


class ReferenceBook:
    """逐设备维护 epoch 与冻结参考，可快照恢复。"""

    def __init__(self, reference_size):
        if reference_size < 2:
            raise ValueError('参考窗口至少需要 2 条观测')
        self.reference_size = reference_size
        self._all = {}
        self._current = {}

    @staticmethod
    def key(machine, epoch):
        return machine, epoch

    def epoch_state(self, machine, epoch):
        return self._all.get(self.key(machine, epoch))

    def iter_epochs(self):
        """按 (machine, epoch) 顺序遍历全部 epoch 状态。"""
        for key in sorted(self._all):
            yield self._all[key]

    def active_epoch(self, machine):
        current = self._current.get(machine)
        return None if current is None else current.epoch

    def reference_for(self, machine, epoch):
        epoch_state = self.epoch_state(machine, epoch)
        if epoch_state is None:
            return None
        return epoch_state.reference

    def observe(self, record):
        """录入一条观测并维护 epoch 结构。

        返回 ``(epoch, status)``，status 取值：

        * ``reference_member``：属于当前 epoch 的参考窗口；
        * ``maintenance``：检修记录本身，同时开启下一 epoch；
        * ``testable``：可对已冻结参考做检验；
        * ``cold_start``：当前 epoch 参考尚未冻结；
        * ``late_closed``：迟到记录落入已关闭且始终未冻结的 epoch，
          无法补建参考。
        """
        machine = record['machine']
        tick = record['tick']
        current = self._current.get(machine)
        if current is None:
            current = EpochState(epoch=0, machine=machine, start_tick=tick)
            self._current[machine] = current
            self._all[self.key(machine, 0)] = current
        if record['maintenance'] == 1:
            if not current.frozen:
                current.closed = True
                current.unfrozen_member_ids = list(current.window_ids)
            successor = EpochState(epoch=current.epoch + 1, machine=machine,
                                   start_tick=tick + 1)
            self._current[machine] = successor
            self._all[self.key(machine, successor.epoch)] = successor
            return current.epoch, 'maintenance'
        if tick < current.start_tick:
            return self._route_late(machine, tick)
        if not current.frozen:
            current.finalized = False
            current.window_values.append(record['residual'])
            current.window_ids.append(record['observation_id'])
            current.window_ticks.append(tick)
            if len(current.window_values) >= self.reference_size:
                self._freeze(current)
            return current.epoch, 'reference_member'
        return current.epoch, 'testable'

    def _route_late(self, machine, tick):
        for epoch_index in range(self._current[machine].epoch - 1, -1, -1):
            epoch = self.epoch_state(machine, epoch_index)
            if epoch is not None and tick >= epoch.start_tick:
                if epoch.frozen:
                    return epoch.epoch, 'testable'
                return epoch.epoch, 'late_closed'
        return 0, 'late_closed'

    def _freeze(self, epoch):
        values = epoch.window_values
        mean = math.fsum(values) / len(values)
        variance = math.fsum((value - mean) ** 2 for value in values) / (len(values) - 1)
        stddev = variance ** 0.5
        if stddev <= 0.0:
            return
        epoch.reference = FrozenReference(
            epoch=epoch.epoch, machine=epoch.machine,
            start_tick=epoch.start_tick, end_tick=epoch.window_ticks[-1],
            size=len(values), mean=mean, stddev=stddev,
            observation_ids=list(epoch.window_ids))

    def snapshot(self):
        payload = {'reference_size': self.reference_size, 'machines': {}}
        for (machine, index), epoch in sorted(self._all.items()):
            bucket = payload['machines'].setdefault(machine, [])
            bucket.append({'epoch': epoch.epoch, 'start_tick': epoch.start_tick,
                           'window_values': list(epoch.window_values),
                           'window_ids': list(epoch.window_ids),
                           'window_ticks': list(epoch.window_ticks),
                           'closed': epoch.closed,
                           'finalized': epoch.finalized,
                           'unfrozen_member_ids': list(epoch.unfrozen_member_ids),
                           'reference': None if epoch.reference is None
                           else epoch.reference.to_dict()})
        return payload

    @classmethod
    def restore(cls, payload):
        book = cls(reference_size=payload['reference_size'])
        for machine, epochs in payload['machines'].items():
            for item in sorted(epochs, key=lambda part: part['epoch']):
                epoch = EpochState(epoch=item['epoch'], machine=machine,
                                   start_tick=item['start_tick'],
                                   window_values=list(item['window_values']),
                                   window_ids=list(item['window_ids']),
                                   window_ticks=list(item['window_ticks']),
                                   closed=item.get('closed', False),
                                   finalized=item.get('finalized', False),
                                   unfrozen_member_ids=list(
                                       item.get('unfrozen_member_ids', [])))
                if item['reference'] is not None:
                    epoch.reference = FrozenReference.from_dict(item['reference'])
                book._all[book.key(machine, item['epoch'])] = epoch
                current = book._current.get(machine)
                if current is None or item['epoch'] > current.epoch:
                    book._current[machine] = epoch
        return book
