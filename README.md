# adaptive-fdr-monitor

制造分析整理设备残差、检修标记和可见时间，按产线计算绝对残差与观测覆盖，并提供基于在线 FDR 预算的持续异常筛查。

## 使用

运行需要 Python 3.11 或更新版本，没有第三方运行依赖。数据在本地处理，不请求外部接口。

```bash
# 原有残差汇总（保持不变）
python3 -m analysis_app summary
python3 -m analysis_app summary --input data/observations.csv --output output/summary.json

# 持续筛查：在线 FDR 预算管理，状态可中断恢复
python3 -m analysis_app monitor \
  --family-weight F0=0.5 --family-weight F1=0.5 \
  --state output/monitor-state.json \
  --output output/monitor-report.json
```

标准输出是完整 JSON，输出文件使用临时文件原子替换；失败退出码为 2，不覆盖既有结果。空数据集保留零计数与空分组，缺失数值不进入数值分母。数据指纹由按 observation_id 排序后的解析记录计算。

### monitor 子命令参数

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `--global-alpha` | 0.05 | 全局 FDR 预算，按固定 family 权重切分 |
| `--reference-size` | 4 | 每个检修 epoch 建立冻结参考所需的稳定窗口观测数 |
| `--lord-share` | 0.8 | family 内 LORD3 流的预算比例，其余给在线 Bonferroni 降级流 |
| `--ar-cutoff` | 0.5 | 滞后 1 阶自相关绝对值的降级阈值 |
| `--ar-min-samples` | 30 | 依赖诊断滚动窗口长度（样本数） |
| `--family-weight` | 缺省等分 | 固定 family 预算权重，可重复；未登记 family 的记录隔离 |
| `--force-conservative` | 关 | 全部检验直接走在线 Bonferroni 流 |
| `--state` | output/monitor-state.json | 可恢复状态文件；存在则续写，不存在则初始化 |

报告包含五部分：`budget_ledger`（每次假设的 p 值、阈值、流、零假设表述与成立条件）、`discoveries`（发现记录）、`dispositions`（每条输入的处置：tested / cold_start / reference_member / maintenance / duplicate / late_arrival / quarantined）、`current_state`（水位线、预算位置、各流财富、参考、降级状态）、`summary`（原有残差汇总，口径不变）。

## 方法学要点

* **主规则 LORD3**：在线 FDR 规则，gamma_t = 6/(pi^2 t^2) 可求和、W0 = alpha/2；适用假设为 p 值独立或 PRDS 正向依赖、真零假设下超均匀。
* **保守降级在线 Bonferroni**：beta_t = alpha * gamma_t，任意依赖下按并集界控制 FWER；自相关诊断连续 3 次越界即切换，连续 3 次恢复可切回。诊断只负责触发降级，FDR 成立性来自规则本身的数学保证，不来自算法名称。
* **层级预算**：全局 alpha 按固定 family 权重、family 内按固定流比例切分，层级并集界给出全局 FDR <= global_alpha；预算位置只增不改，迟到记录隔离不插入历史位置。
* **参考隔离**：按检修切 epoch，参考仅来自本 epoch 稳定窗口，窗口不足给冷启动；检修重建参考但不清空全局错误账本。

## 长期模拟验证

```bash
python3 -m scripts.simulate_monitor --trials 30 --ticks 400 \
  --output docs/simulation-results.json
```

四个确定性场景（固定随机种子，仅标准库）：

1. `iid_null`：独立零假设，验证长期误报被 global_alpha 控制；
2. `ar_null`：AR(1) 相关零假设，验证自相关诊断触发 Bonferroni 降级后误报仍有界；
3. `early_burst`：参考冻结后的强突发信号后回到零假设，验证 LORD 财富补充与长期误报；
4. `alternating`：检修后持续均值漂移与检修切换，验证参考重建、发现能力与全局账本保留。

每个场景同时给出固定 3 倍标准差规则的对照（无错误预算管理，长序列下累计大量误报）。最近一次 30 试验 × 400 tick 的结果保存在 `docs/simulation-results.json`。

## 验证与容器

```bash
python3 -m unittest discover -s tests -v
python3 -m compileall -q analysis_app scripts
docker build -t adaptive-fdr-monitor .
docker run --rm adaptive-fdr-monitor
```

容器内也可执行 `python -m unittest discover -s tests -v`。镜像构建阶段取得 Python 镜像，产品启动、测试和演示均无需下载。小规模样例的完整单元测试和 summary 可以连续执行，不需要等待定时任务。

## 代码与数据

- `analysis_app/io.py`：有类型的 CSV 读取、完整性校验及内容指纹。
- `analysis_app/schema.py`：字段类型、主体与聚合维度。
- `analysis_app/domain.py`：当前业务核算函数。
- `analysis_app/summary.py`：确定性的分组与数值口径。
- `analysis_app/stats.py`：Student-t 双侧 p 值与滞后 1 阶自相关（仅标准库）。
- `analysis_app/fdr.py`：LORD3 与在线 Bonferroni 包络的顺序控制器。
- `analysis_app/references.py`：按检修 epoch 切分的冻结参考窗口管理。
- `analysis_app/monitor.py`：监控器配置、检验排序、层级预算、降级诊断、状态快照与恢复。
- `analysis_app/cli.py`：本地命令及原子输出。
- `scripts/simulate_monitor.py`：长期误报、降级与恢复的模拟验证。
- `docs/data.md`：数据契约与持续筛查语义。

演示数据全部为本地合成，不对应现实机构、设备或人员；不能用来推断真实业务效果。
