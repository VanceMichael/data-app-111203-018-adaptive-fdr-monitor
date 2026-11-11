# 数据契约

制造分析整理设备残差、检修标记和可见时间，按产线计算绝对残差与观测覆盖。

主测量 `residual` 的单位为标准化残差，主体字段是 `machine`，核算维度是 `line`。`observation_id` 是文件内唯一观测编号，主体允许重复观测。数值只接受有限值；空单元格解析为 None，不能直接当成零。额外列不进入当前核算。

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `observation_id` | str | 必填 |
| `machine` | str | 必填 |
| `line` | str | 必填 |
| `tick` | int | 允许空值，空值不是零 |
| `residual` | float | 必填 |
| `maintenance` | int | 允许空值，空值不是零 |
| `visible_at` | int | 允许空值，空值不是零 |
| `family` | str | 允许空值，空值不是零 |

时间字段 day、week、tick、hour、minute、cycle、start_day、reported_day、known_at、visible_at、opened、observed 均以各字段名所述的相对时间单位表示整数；不同字段并不必然代表同一种时间。缺少的登记码与浓度等信息保持为空。示例只表达产品数据结构，不能代替实际调查资料。
