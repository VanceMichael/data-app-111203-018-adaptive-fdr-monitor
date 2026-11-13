import argparse
import json
import os
import sys
import tempfile
from pathlib import Path

from .io import load_records
from .monitor import Monitor, load_state
from .summary import summarize


def save_json(path, result):
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix='.analysis-', dir=destination.parent)
    try:
        with os.fdopen(descriptor, 'w', encoding='utf-8') as handle:
            json.dump(result, handle, ensure_ascii=False, sort_keys=True, allow_nan=False, indent=2)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def build_parser():
    parser = argparse.ArgumentParser(description='本地观测数据分析')
    commands = parser.add_subparsers(dest='command', required=True)
    summary = commands.add_parser('summary', help='按业务维度核对观测口径')
    summary.add_argument('--input', type=Path, default=Path('data/observations.csv'))
    summary.add_argument('--output', type=Path)

    monitor = commands.add_parser('monitor', help='持续筛查：在线 FDR 预算管理与异常发现')
    monitor.add_argument('--input', type=Path, default=Path('data/observations.csv'))
    monitor.add_argument('--output', type=Path, help='本次运行报告（预算流水、发现、处置与状态）的 JSON 输出路径')
    monitor.add_argument('--state', type=Path, default=Path('output/monitor-state.json'),
                         help='可恢复状态文件；存在则续写，不存在则初始化')
    monitor.add_argument('--global-alpha', type=float, default=0.05)
    monitor.add_argument('--reference-size', type=int, default=4)
    monitor.add_argument('--lord-share', type=float, default=0.8)
    monitor.add_argument('--ar-cutoff', type=float, default=0.5)
    monitor.add_argument('--ar-min-samples', type=int, default=30)
    monitor.add_argument('--family-weight', action='append', default=[],
                         metavar='FAMILY=WEIGHT', help='固定 family 预算权重，可重复；缺省按首次出现等分')
    monitor.add_argument('--force-conservative', action='store_true',
                         help='全部检验直接走在线 Bonferroni 保守流')
    return parser


def parse_family_weights(tokens):
    weights = {}
    for token in tokens:
        if '=' not in token:
            raise ValueError('family 权重格式应为 FAMILY=WEIGHT：' + token)
        family, raw = token.split('=', 1)
        weights[family] = float(raw)
    return weights


def run_monitor(args):
    from .monitor import MonitorConfig
    config = MonitorConfig(
        global_alpha=args.global_alpha, reference_size=args.reference_size,
        lord_share=args.lord_share, ar_cutoff=args.ar_cutoff,
        ar_min_samples=args.ar_min_samples,
        family_weights=parse_family_weights(args.family_weight),
        force_conservative=args.force_conservative)
    config.validate()
    records = load_records(args.input)
    if args.state.exists():
        monitor = Monitor.restore(load_state(args.state))
        if monitor.config.to_dict() != config.to_dict():
            raise ValueError('配置与既有状态不一致；错误预算规则初始化后不可更改')
    else:
        monitor = Monitor(config=config)
    from .io import fingerprint
    finger = fingerprint(sorted(records, key=lambda row: row['observation_id']))
    monitor.ingest(records, input_fingerprint=finger)
    monitor.finalize()
    report = monitor.report(input_fingerprint=finger, rows_total=len(records))
    report['summary'] = summarize(records)
    if args.output:
        save_json(args.output, report)
    monitor.save_state(args.state)
    print(json.dumps(report, ensure_ascii=False, sort_keys=True, allow_nan=False))
    return 0


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == 'summary':
            result = summarize(load_records(args.input))
            if args.output:
                save_json(args.output, result)
            print(json.dumps(result, ensure_ascii=False, sort_keys=True, allow_nan=False))
            return 0
        return run_monitor(args)
    except (ValueError, OSError) as error:
        print(str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
