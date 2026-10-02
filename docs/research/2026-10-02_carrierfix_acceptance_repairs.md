# 载频复审三项修复与验收（2026-10-02）

用户在复审后授权自主修复。本单元基于 main `5af6826` 和独立复审 `6ae9c12`，在 `codex/carrierfix-acceptance-20261002` 完成。应用 `receiving-code-review` 技能，逐项按归档证据确认后修改；未委派草稿。

## 修复结果

1. **阈值结论**：报告、START_HERE、decision_log、registry 均更正为 2 样本 D=0.955661/0.998495/1.010644 全部超过 0.95，最小裕量 0.005661；没有依据认定阈值必须修改。冻结容差保持原值。
2. **排序与比值结论**：仅本次四族 top-1 与可行性零翻转。完整排序的 S2X B7/B5 和 S2TZX B8 相对 B4/B3 变化已登记；加噪 Nb 和到时诊断的变化如实保留。不能将全局幅度缩放抵消扩大到载频变化、全尺度或截幅。增益/截幅及增广的旧链未由此单元升级为正式可用证据。
3. **复算验收**：增加 `scripts/sfcw_carrierfix_acceptance.py`，三组脚本在写文件前验证旧表示指标/峰值/选型与归档完全一致。缺行、重复行或失配即抛错。S1/S3 还核对冻结费马表、峰窗和 NC 净空，使用冻结的逐道 H5 哈希和恢复目录；不再优先选择任意可找到的 attempt。t3/B2 要求独立预先建立的逐道输入清单，不以正在验收的文件自生成历史身份。

官方 `processing.py` 在每次加载前核验完整 SHA-256；逐道 H5 在加载前后核验，拒绝并发变动。产物记录实际输入与恢复目录、脚本、加载器、相关契约/归档、Python/NumPy 和官方源码版本/哈希。新验收使用显式异常，在 `python -O` 下也生效。

加强后的脚本输出采用新的 `carrierfix_v0_3_{r1,r2}.json` 和 schema `/2`，仅在两遍全部验收及字节一致后创建新文件；既有输出一律拒绝覆盖。旧 v0.2 六份结果保留，不能补写并不存在的历史来源哈希。

## 验证与证据边界

- **17 项验收反例通过**：归档正常记录通过；改变旧指标/峰值/选型、缺行/重复行、遗漏字段、错误排名标记、H5 字节或目录不符、缺少/越界清单、官方代码变化、峰窗/NC 越界、两遍不一致及覆盖历史文件均被拒绝。另直接执行三组脚本的错误旧值路径，确认没有输出文件。普通模式和 `-O` 模式均通过。
- **独立数值检查通过**：仍有 264 项 t3 旧指标精确一致、B2 已报告字段与归档一致、132 个旧峰登记一致；80 MHz 构造输入和本机控制道的独立直接 DFT/逆求和重建通过。前述验证不是 FDTD 收敛或仪器标定。
- **182 项既有纯数组回归通过**：`python scripts/verify_workspace.py`，结果在忽略的 `artifacts/local_checks/20261002T080150Z-c1cadc7f/`。
- 新证据：`artifacts/research_checks/2026-10-02_carrierfix_acceptance_checks.json`，包含验收计数、新代码哈希、旧新归档核对、实际官方加载器及单控制道检查。

本机没有完整 t3/S1S3/B2 原始 H5，也没有 t3/B2 的完整历史逐道 H5 哈希清单，因此**没有生成全量 v0.3 科学结果**。验收代码及错误拒绝行为已验证；全量原始输入重跑仍待原始归档所在主机完成，不能冒充本机已经端到端复现。G4、冻结契约、测试族与求解/训练门未改，未运行求解器或训练。

## 复跑入口

```powershell
& 'artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe' scripts/check_sfcw_carrierfix_acceptance.py
```

默认验证产物写入忽略的 local_checks，拒绝覆盖同名文件。完整开发侧复算仅在拥有原始归档和可追溯输入清单时进行：

```powershell
& 'artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe' scripts/study_t3_ladder_carrierfix_v0_2.py --input-manifest <t3清单全路径>
& 'artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe' scripts/study_b2_reward_carrierfix_v0_2.py --input-manifest <B2清单全路径>
& 'artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe' scripts/study_s1s3_refwindow_carrierfix_v0_2.py
```

清单 JSON 格式为 `{"inputs":[{"run_id":"...","source_dir":"...","h5_sha256":"64位小写十六进制"}]}`，精确覆盖对应脚本要求的 99/132 个开发 run_id，无额外案例；source_dir 是 `artifacts/simulations/` 下归档目录的单个名称。哈希须来自预先核验的归档身份依据，不从本次待验收输出倒推。S1/S3 的清单已经在冻结契约中。若新主机缺少这些来源，脚本失败退出；不启动求解器补造输入。

后续先在原始归档所在主机完成 v0.3 后处理验收，再复核官方表示的敏感性与下游增益/截幅证据。此修复单元不关闭上轮几何、PML、2D/3D 配对、材料及学习设计问题。
