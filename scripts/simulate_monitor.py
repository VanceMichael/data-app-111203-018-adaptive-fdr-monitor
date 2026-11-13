"""adaptive-fdr-monitor 的长期模拟验证（仅使用标准库，确定性随机源）。

三个场景：

1. ``iid_null``：所有设备残差独立同分布 N(0, 1)，参考窗口从同一稳定
   历史窗口冻结。零假设全真，长期 FDR（这里等价于误报率）应被
   global_alpha 控制；同时给出固定 3 倍标准差规则的对照，展示固定
   阈值在高频检验下的累计误报问题。
2. ``ar_null``：零假设仍全真，但设备残差为 AR(1) 相关过程。滞后 1 阶
   自相关超过 ar_cutoff 后设备黏性降级到在线 Bonferroni 流，验证降级
   路径被触发且整体误报仍有界。
3. ``alternating``：部分设备在检修后进入均值漂移的备择状态，验证
   发现能力（power），以及检修只重建参考、不清空全局错误账本。

用法：
    python -m scripts.simulate_monitor --trials 30 --ticks 400 --output output/simulation.json
"""
import argparse
import json
import math
import random
from pathlib import Path

from analysis_app.monitor import Monitor, MonitorConfig

MACHINES = 12
FAMILIES = ['F0', 'F1']
REFERENCE_SIZE = 20
AR_MIN_SAMPLES = 60
AR_CUTOFF = 0.5


def make_rows(rng, ticks, mode, maintenance_tick=None):
    """生成合成观测。

    * ``iid_null`` / ``ar_null``：零假设全真，后者为 AR(1) 相关过程；
    * ``early_burst``：参考窗口冻结后（tick 30-99）部分设备存在强突发
      信号（真实备择），之后全部回到零假设，用于验证 LORD 早期拒绝
      补充财富后，长期误报仍然受控；
    * ``alternating``：三分之一设备在检修后产生持续 4σ 均值漂移，
      另三分之一为 2.5σ 中等漂移；漂移从新 epoch 参考窗口冻结后开始，
      参考不会吸收备择。
    """
    rows = []
    ar_state = {f'M{m:02d}': 0.0 for m in range(MACHINES)}
    ar_phi = 0.6 if mode == 'ar_null' else 0.0
    drift_machines = {f'M{m:02d}' for m in range(MACHINES) if m % 3 in (0, 1)}
    for tick in range(ticks):
        for m in range(MACHINES):
            machine = f'M{m:02d}'
            family = FAMILIES[m % 2]
            maintenance = 1 if (maintenance_tick is not None
                                and tick == maintenance_tick) else 0
            if mode == 'ar_null':
                noise = rng.gauss(0.0, math.sqrt(1.0 - ar_phi ** 2))
                ar_state[machine] = ar_phi * ar_state[machine] + noise
                residual = ar_state[machine]
            else:
                # iid_null、early_burst 与 alternating 均使用独立高斯噪声；
                # 漂移是确定性均值项，不引入序列相关
                residual = rng.gauss(0.0, 1.0)
            if mode == 'early_burst' and 30 <= tick < 100 and m % 4 == 0:
                residual += 5.0
            if (mode == 'alternating' and maintenance_tick is not None
                    and tick > maintenance_tick + REFERENCE_SIZE):
                if m % 3 == 0:
                    residual += 4.0      # 强持续漂移：应在预算允许时被发现
                elif m % 3 == 1:
                    residual += 2.5      # 中等漂移：边缘可检出，展示功效取舍
            rows.append({'observation_id': f'{machine}-{tick}',
                         'machine': machine, 'family': family,
                         'line': f'L{m % 3}', 'tick': tick,
                         'residual': residual, 'maintenance': maintenance,
                         'visible_at': tick})
    return rows


def run_trial(seed, mode, ticks, global_alpha):
    rng = random.Random(seed)
    maintenance_tick = ticks // 2 if mode == 'alternating' else None
    rows = make_rows(rng, ticks, mode, maintenance_tick)
    config = MonitorConfig(global_alpha=global_alpha,
                           reference_size=REFERENCE_SIZE,
                           ar_min_samples=AR_MIN_SAMPLES,
                           ar_cutoff=AR_CUTOFF,
                           family_weights={'F0': 0.5, 'F1': 0.5})
    monitor = Monitor(config=config)
    # 模拟中断恢复：分三批摄入，每批之间做快照/恢复
    boundaries = [len(rows) // 3, 2 * len(rows) // 3]
    chunks = [rows[:boundaries[0]], rows[boundaries[0]:boundaries[1]],
              rows[boundaries[1]:]]
    for chunk in chunks:
        monitor.ingest(chunk, f'batch-{seed}')
        monitor = Monitor.restore(monitor.state_snapshot())
    monitor.finalize()
    tested = monitor.global_position
    burst_machines = {f'M{m:02d}' for m in range(MACHINES) if m % 4 == 0}
    if mode == 'alternating':
        drift_machines = {f'M{m:02d}' for m in range(MACHINES) if m % 3 in (0, 1)}
        false_discoveries = sum(
            1 for item in monitor.discoveries
            if item['tick'] <= maintenance_tick
            or item['machine'] not in drift_machines)
        true_discoveries = len(monitor.discoveries) - false_discoveries
        alternative_tests = sum(
            1 for entry in monitor.budget_ledger
            if entry['machine'] in drift_machines
            and entry['tick'] > maintenance_tick + REFERENCE_SIZE)
    elif mode == 'early_burst':
        # tick 30-99 且突发机上的发现为真发现，其余全部为误报
        false_discoveries = sum(
            1 for item in monitor.discoveries
            if not (30 <= item['tick'] < 100 and item['machine'] in burst_machines))
        true_discoveries = len(monitor.discoveries) - false_discoveries
        alternative_tests = sum(
            1 for entry in monitor.budget_ledger
            if 30 <= entry['tick'] < 100 and entry['machine'] in burst_machines)
    else:
        false_discoveries = len(monitor.discoveries)
        true_discoveries = 0
        alternative_tests = 0
    degraded = sorted(monitor.current_state()['degraded_machine_epochs'])
    bonferroni_tests = sum(1 for entry in monitor.budget_ledger
                           if entry['stream'] == 'online_bonferroni')
    lord_tests = sum(1 for entry in monitor.budget_ledger
                     if entry['stream'] == 'lord3')
    # 固定 3 倍标准差对照（用每机 epoch0 前 REFERENCE_SIZE 条做参考）
    ref_stats = {}
    for machine in {row['machine'] for row in rows}:
        window = [row['residual'] for row in rows
                  if row['machine'] == machine and row['tick'] < REFERENCE_SIZE]
        mean = math.fsum(window) / len(window)
        var = math.fsum((v - mean) ** 2 for v in window) / (len(window) - 1)
        ref_stats[machine] = (mean, math.sqrt(var))
    fixed_alarms, fixed_tests = 0, 0
    for row in rows:
        if row['maintenance'] == 1 or row['tick'] < REFERENCE_SIZE:
            continue
        mean, sd = ref_stats[row['machine']]
        fixed_tests += 1
        if abs((row['residual'] - mean) / sd) > 3.0:
            fixed_alarms += 1
    return {
        'seed': seed, 'mode': mode, 'tested': tested,
        'discoveries': len(monitor.discoveries),
        'false_discoveries': false_discoveries,
        'true_discoveries': true_discoveries,
        'alternative_tests': alternative_tests,
        'fdp': (false_discoveries / len(monitor.discoveries)
                if monitor.discoveries else 0.0),
        'any_false_alarm': false_discoveries > 0,
        'lord_tests': lord_tests, 'bonferroni_tests': bonferroni_tests,
        'degraded_machines': degraded,
        'fixed_3sigma_alarms': fixed_alarms,
        'fixed_3sigma_tests': fixed_tests,
        'watermark': monitor.watermark,
    }


def aggregate(trials, mode):
    count = len(trials)
    tested = sum(t['tested'] for t in trials)
    false_total = sum(t['false_discoveries'] for t in trials)
    disc_total = sum(t['discoveries'] for t in trials)
    result = {
        'mode': mode, 'trials': count,
        'tests_per_trial': trials[0]['tested'],
        'marginal_false_rate': false_total / max(tested, 1),
        'trial_level_false_alarm_probability': (
            sum(1 for t in trials if t['any_false_alarm']) / count),
        'mean_discoveries_per_trial': disc_total / count,
        'mean_false_discoveries_per_trial': false_total / count,
        'trials_with_degradation': sum(1 for t in trials if t['degraded_machines']),
        'bonferroni_share_of_tests': (
            math.fsum(t['bonferroni_tests'] for t in trials)
            / max(math.fsum(t['tested'] for t in trials), 1)),
        'fixed_3sigma_alarm_rate': (
            math.fsum(t['fixed_3sigma_alarms'] for t in trials)
            / max(math.fsum(t['fixed_3sigma_tests'] for t in trials), 1)),
        'fixed_3sigma_mean_alarms_per_trial': (
            math.fsum(t['fixed_3sigma_alarms'] for t in trials) / count),
    }
    alt_tests = math.fsum(t['alternative_tests'] for t in trials)
    if alt_tests:
        result['alternative_power'] = (
            math.fsum(t['true_discoveries'] for t in trials) / alt_tests)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description='在线 FDR 监控器长期模拟')
    parser.add_argument('--trials', type=int, default=30)
    parser.add_argument('--ticks', type=int, default=400)
    parser.add_argument('--global-alpha', type=float, default=0.05)
    parser.add_argument('--seed-base', type=int, default=20261113)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args(argv)
    all_trials = {}
    for mode in ('iid_null', 'ar_null', 'early_burst', 'alternating'):
        trials = [run_trial(args.seed_base + index, mode, args.ticks,
                            args.global_alpha)
                  for index in range(args.trials)]
        all_trials[mode] = {'aggregate': aggregate(trials, mode),
                            'trials': trials}
    payload = {'config': {'trials': args.trials, 'ticks': args.ticks,
                          'global_alpha': args.global_alpha,
                          'machines': MACHINES, 'reference_size': REFERENCE_SIZE,
                          'ar_cutoff': AR_CUTOFF, 'ar_min_samples': AR_MIN_SAMPLES,
                          'seed_base': args.seed_base},
               'scenarios': all_trials}
    text = json.dumps(payload, ensure_ascii=False, sort_keys=True, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + '\n', encoding='utf-8')
    print(text)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
