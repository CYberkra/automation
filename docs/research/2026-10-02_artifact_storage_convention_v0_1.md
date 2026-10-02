# 仿真产物长期存储规范 v0.1

日期：2026-10-02
状态：**v0.1 生效**（用户"有必要创建一个长期的存储的规范和方案吧"授权）
首个执行案例：`artifacts/research_checks/2026-10-02_halfspace_standard_hs/`（HS1/HS2/HS3 半空间标准件批）

---

## 1. 问题

HS 批之前，重要仿真产出的存放是即兴的：求解输出落在仓库外的 `verify_runs/`（git 不可见），验收分析用 heredoc 现场跑（脚本不落盘、无法复现），图散在工作区根目录。一旦本机文件清理或人员交接，**结论无法复现、证据链断裂**。

## 2. 三层存储模型

| 层 | 位置 | git | 内容 | 准则 |
|---|---|---|---|---|
| **T1 追踪层** | `artifacts/research_checks/<date>_<name>/` | ✅ 追踪 | .in 输入、小型 H5 迹线（单道/少道，<10 MB/件）、验收/分析脚本引用、npz 派生数据、验收图、manifest.json | 一切进入结论/台账的证据必须在 T1 有锚点 |
| **T2 本地层** | `artifacts/simulations/`、`artifacts/local_checks/`、`verify_runs/`（仓库外） | ❌ gitignore | 批量原始 H5（33 道×N 模型）、vtkhdf 几何转储（>100 MB）、中间缓存 | 可丢弃；必须能由 T1 的 .in + 脚本重建 |
| **T3 远端层** | GitHub `CYberkra/automation` main | — | T1 的全部内容随提交推送 | 推送前必须过"复现性自检"（§4） |

**单件体积阈值**：≤10 MB 进 T1；>10 MB 留 T2，T1 只存其 SHA256 + 生成路径。阈值以外的例外须记入决策日志。

## 3. T1 证据目录的内部构成（强制）

一个证据目录 = 一个自包含的"实验胶囊"：

```
<date>_<name>/
├── *.in                    # 求解器输入（逐字，含冻结数值设置）
├── *.h5                    # 原始迹线（小件时；大件则在 manifest 记 SHA + T2 路径）
├── manifest.json           # 目录内全部文件的 SHA256（自身除外）
├── *_metrics.json          # 验收数字（机器可读，台账 key_numbers 的来源）
└── *.png                   # 验收图（已目验版本）
```

配套要求：

- **验收脚本进 `scripts/`**（如 `run_hs_acceptance_v0_1.py`），参数化输入目录，禁止 heredoc 现场分析进结论；
- **台账条目** `research_artifact_registry.json` 的 `key_numbers` 逐字引用 metrics JSON 的数字，`directory` 指向证据目录；
- **决策日志**新条目插文件顶部，写清：谁授权的、跑了什么、验收数字、下一步。

## 4. 复现性自检（提交前必做）

1. 用证据目录里的副本（而非工作原件）重跑验收脚本，metrics 数字与首次运行一致；
2. manifest.json 重建后无漂移；
3. 图已目验（ReadMediaFile），无字形缺失/坐标轴颠倒；
4. `git status` 确认 T2 大文件未被误加。

## 5. 既有约定的衔接

- `research_checks/` 历史证据维持追踪（本规范追认其地位）；
- 冻结契约内嵌 SHA 链的文件**原地归档、永不改写**（见 2026-10-02 旧载频归档条目），新证据另开目录；
- 批量仿真（如 T1/T3 396 例）原始 H5 天然属于 T2，T1 只收其分析产物与冻结清单。

## 6. 待办（规范自身的演进）

- T2 层的异地备份策略未定（本机单点，3.3 GB 批量 H5 无副本）；
- manifest 自动化校验脚本未写（目前手动 python 一行）；
- 图的字形规范（中文用 msyh FontProperties，禁用上标字符）可抽成公共绘图模块。

## 7. 首次存量清理（2026-10-02 执行，用户"存量老 npz 你自己清理好"+"重要的保留、过时的清理"）

按"冻结契约 SHA 内嵌 → 必须保留；未内嵌且 >10 MB → 降级 T2"的口径核查全部 >5 MB 已追踪文件：**降级 33 件共 459.1 MB**（09-26 损伤阶梯 r1/r2 × 16、09-27 dev 阶梯 r1/r2 × 16 的 records.json，以及 22.5 MB 的统一适配器 npz——均为旧载频时代证据或可再生转储，本地文件原地保留、SHA256 全部入 [降级清单](2026-10-02_large_legacy_demotion_manifest.json)）；**保留**：损伤阶梯 test 族 records.json（SHA 内嵌于 g4_s5_test_confirmation 等冻结契约）、全部 MT33 母模型 H5（重跑需真实 GPU 成本）、r1/r2 成对 arrays.npz（复现性证据，<10 MB）。注意：降级只影响未来克隆的追踪集，**git 历史包体积不变**（改写历史需强推，未做）。
