# 数据契约

制造分析整理设备残差、检修标记和可见时间，按产线计算绝对残差与观测覆盖。

主测量 `residual` 的单位为标准化残差，主体字段是 `machine`，核算维度是 `line`。`observation_id` 是文件内唯一观测编号，主体允许重复观测。数值只接受有限值；空单元格解析为 None，不能直接当成零。额外列不进入当前核算。

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `observation_id` | str | 必填 |
| `machine` | str | 必填 |
| `line` | str | 必填 |
| `tick` | int | 必填 |
| `residual` | float | 必填 |
| `maintenance` | int | 必填 |
| `visible_at` | int | 必填 |
| `family` | str | 必填 |

时间字段 day、week、tick、hour、minute、cycle、start_day、reported_day、known_at、visible_at、opened、observed 均以各字段名所述的相对时间单位表示整数；不同字段并不必然代表同一种时间。缺少的登记码与浓度等信息保持为空。示例只表达产品数据结构，不能代替实际调查资料。

## 业务字段语义

residual 为标准化有符号残差，maintenance 为检修事件标志。tick 是测量时间，visible_at 是可见时间，均以小时表示。family 为设备族；line 为产线，检修后参考分布可能变化。

## tick 与 visible_at 的区别

* `tick` 是残差实际产生的测量小时，用于在设备时间轴上界定参考窗口与检修 epoch；
* `visible_at` 是该观测进入监控器、可以参与检验的可见小时。同一 tick 的观测可能因上报延迟而具有更大的 `visible_at`；
* 持续筛查的检验顺序只由 `(visible_at, machine, family, tick, observation_id)` 决定，与文件中的行顺序无关；
* 已经推进的可见时间水位线（watermark）之后才到达、且 `visible_at` 低于水位线的记录属于迟到记录，只能隔离（`late_arrival`），不能插入已经花费错误预算的历史位置，也不会回写任何已冻结参考。

## 持续筛查（monitor）契约

### 参考分布与检修 epoch

* 每台设备按 `maintenance=1` 事件把时间轴切成 epoch；检修记录本身不做检验，下一 tick 开启新 epoch；
* 每个 epoch 的参考分布只能由该 epoch 开头 `reference_size` 条 `maintenance=0` 观测（稳定窗口）建立，窗口收齐即冻结，检修前后窗口永不混用；
* 参考窗口不足、窗口内残差无波动（标准差为 0）或 epoch 在窗口收齐前关闭时，相关观测进入 `cold_start` 状态：不产生 p 值、不消费错误预算；
* 检验对冻结参考做 Student-t 双侧检验（自由度 n-1，含参考估计不确定性），参考的均值、标准差、成员 observation_id 与 tick 区间随每次假设写入预算流水。

### 错误预算规则与适用假设

* 主规则为 LORD3 在线 FDR（可求和序列 gamma_t = 6/(pi^2 t^2)，初始财富 W0 = alpha/2），适用假设：流内 p 值独立或满足 PRDS 正向依赖，且真零假设下 p 值超均匀；
* 降级规则为在线 Bonferroni 包络（beta_t = alpha * gamma_t），仅要求超均匀性，任意序列依赖下按并集界控制 FWER，因而也控制 FDR；
* 降级触发：按 (设备, epoch) 维护滚动窗口的滞后 1 阶自相关，去均值后绝对值连续 3 次超过 `ar_cutoff` 即降级，连续 3 次回到阈值内可恢复；`force_conservative=true` 时全部检验直接走 Bonferroni 流；
* 全局预算 `global_alpha` 按初始化时冻结的固定 family 权重切分，family 内再按 `lord_share` 切 LORD/Bonferroni 两个流；层级并集界给出全局 FDR <= global_alpha，未登记 family 的记录隔离不补切预算；
* 预算位置（family 内流位置与全局位置）只增不改，状态中断恢复后从快照续写。

### 状态与恢复

* 状态文件包含水位线、已处理 observation_id、待处理缓冲、各 epoch 参考、各流控制器财富与拒绝位置、依赖诊断、预算流水、发现与处置记录；
* 重复 `observation_id` 只追加 `duplicate` 审计处置，不重复消费预算；
* 检修只新建参考 epoch，不清空 family 或全局错误账本；
* 配置在监控器初始化后冻结，恢复时配置不一致直接报错，不允许事后修改错误预算规则。
