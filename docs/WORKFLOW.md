# 版本管理与跨电脑接续

仓库：<https://github.com/CYberkra/automation>。2026-09-24 接入时为私有空仓库。主分支 `main` 保存可接续的研究状态；这里的稳定版本指文档和检查可复现，不代表雷达性能已验证。

## 新电脑开始

**先读 [START_HERE.md](../START_HERE.md)。** 当前设备为SFCW 20–170MHz，项目暂按地下约20m以内理解；gprMax V4仿真须先与用户敲定，当前未批准。该入口列出已完成工作与剩余事项，优先于历史交接。

需要 Git 和 Python 3.10；当前核验环境为 Python 3.10.9、NumPy 2.2.6。先确保新电脑的 GitHub 账号有该私有仓库访问权限，使用它自己的登录凭据。

在 PowerShell 中运行：

```powershell
git clone https://github.com/CYberkra/automation.git
Set-Location automation
py -3.10 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe scripts/verify_workspace.py
```

直接调用虚拟环境的 Python，无需修改 PowerShell 执行策略。若系统没有 `py` 启动器，用已安装的 Python 3.10 路径执行 `-m venv .venv`。Linux/macOS 对应为 `python3.10 -m venv .venv`，后续解释器路径为 `.venv/bin/python`；跨系统数值结果仍需在目标机器实际验证。

当前预期为 **136 项通过**：代数10、算子及参考23、加权评分反例5、评价与标签26、损伤试验72。损伤试验默认不绘图，只依赖NumPy；不读取实测、不联网、不运行gprMax或训练。加权评分仅为反例，用户已撤回“人为权重先行”的临时表述。通过仅证明构造数组行为。

核心依赖仍按requirements.txt固定NumPy2.2.6。新增损伤例另在Python3.12.14/NumPy2.3.5和Python3.10.21/NumPy1.26.4验证。可选绘图依赖见requirements-figures.txt。**gprMax V4源码要求Python3.11–3.13，须另建环境，上述数组环境不是求解器安装步骤。**

2026-09-24 已在独立本地克隆、全新虚拟环境中通过全部 64 项，55 个已跟踪文件与原目录逐字节一致，见[版本接续验证记录](../artifacts/research_checks/2026-09-24_repository_portability.json)。该检查仍在原 Windows 电脑执行，不能替代目标电脑验证。此次 pip 访问索引发生 TLS 中断，后改由 curl 获取官方 PyPI 元数据和安装包，匹配官方 SHA-256 后本地安装；未关闭证书验证。不要把这种网络错误解释成 NumPy 版本不存在。

新结果写入忽略的 `artifacts/local_checks/<UTC时间>-<随机后缀>/`，包含各组 JSON、环境、脚本及结果哈希。历史 `artifacts/research_checks/` 不被覆盖。失败时命令返回非零值，不输出成功汇总；请保留失败输出用于定位。

然后按顺序阅读：

1. 根目录 [START_HERE.md](../START_HERE.md)：最新状态、用户决定和避免重复工作的入口。
2. [最新研究进度](research/continuous_research_status.md)：从最上面的最新记录开始。
3. [评价与标签 v0.2](research/2026-09-24_evaluation_and_labels_v0.2.md)及[机器契约](../configs/research/evaluation_label_contract_v0.2.json)。
4. [损伤试验](research/2026-09-24_damage_pilot_findings.md)与[V4报告](research/2026-09-24_gprmax_v4_review.md)。旧cycle05中的损伤阶梯和D误差包络已完成首版。

新会话可直接使用START_HERE.md内的接续文本。本地对话历史不随Git迁移，决定、失败原因及下一步必须落档；历史报告中的旧待办不自动代表当前状态。

## 日常同步

开始工作先检查：

```powershell
git status --short --branch
git fetch origin
git log --oneline --decorate -5
```

工作区干净、当前在 `main` 且没有分叉时：

```powershell
git pull --ff-only
```

不要覆盖另一台电脑尚未提交的改动。若两台电脑可能同时工作，从更新后的 `main` 建立各自短分支，例如 `git switch -c codex/reference-calibration`；完成后合并并保留历史。`pull --ff-only` 遇到分叉会停止，此时先比较两边改动、解决冲突并重跑相关检查，不使用强制推送或 `reset --hard` 跳过问题。

每个完整研究单元完成后，更新最新进度，运行与改动相关的检查，显式暂存本次文件，检查 `git diff --cached --stat` 和 `git diff --cached`，再提交与推送。提交信息例：`docs(research): record reference uncertainty assumptions`、`test(evaluation): add weak-event damage controls`。首次推送新分支使用 `git push -u origin <分支名>`，已建立上游后使用 `git push`。Git 提交保存在当前电脑；只有推送成功，其他电脑才能取得它。推送失败要明确记录，不报告“已同步”。

当前不自动创建新定时任务，也不迁移用户维护的触发器。切换电脑后由用户决定在哪个会话继续定时工作，避免两处同时推进同一分支。

## 仓库包含与不包含

| 内容 | 管理方式 |
|---|---|
| 工作约定、研究报告、来源台账、JSON 契约、脚本 | Git 跟踪 |
| 自编综述 Word/PDF、少量历史数组结果 JSON | Git 跟踪；原始数值不来自实测测线 |
| `探地雷达背景资料/` | 留在本机，约 1.57 GB；不上传 Git |
| `artifacts/initial_audit/` | 原始资料派生审计、文本和预览，只在本机保留 |
| `artifacts/research_sources/` | 公开文献及数据缓存，不进 Git；来源台账记录 URL、阅读范围和哈希 |
| `.venv/`、`__pycache__/`、`artifacts/local_checks/` | 本机重建 |
| 旧备用调度配置及脚本 | 忽略，不能用克隆操作迁移调度 |

纯研究和136项检查不需要原始资料。确需后期使用时，由用户通过合适渠道另行迁移到同名目录，核对历史清单和实际存在性，仍遵守开发与后期验证分离。历史资料审计不应作为新电脑初始化命令自动执行。

`requirements-audit.txt` 只记录旧资料审计的可选依赖，不代表当前授权重读测线。综述源文档是 Markdown，现有 Word/PDF 已入库；`render_literature_review.cjs` 依赖 Node 和外部 `docx` 包，不属于核心验证入口，尚未提供锁定的渲染环境。

历史报告中的 `D:\自动处理`、`E:\python\python.exe` 或缓存路径是当时的来源记录，新机器无需复建相同盘符。资料审计链接、公开缓存链接和已忽略的调度文件链接在纯克隆中可能不存在；这不表示代码缺失。优先依据台账里的原始 URL 获取所需公开材料，核对哈希；来源更新导致哈希不同需另建版本记录，不能覆盖旧来源记录。当前没有自动下载所有论文的脚本。

## 哈希与证据版本

`.gitattributes` 使用 `* -text`，让 Git 原样保存文件字节，覆盖机器上 `core.autocrlf` 的换行转换行为。已有证据按原始字节计算 SHA-256；不得为了统一换行而批量重写旧文件。新编辑也会产生新哈希，应追加相应版本和检查记录。

历史结果对应历史代码，不能要求旧版 25 项评价输出与后来 26 项脚本哈希相同。当前有效检查结果为 `2026-09-24_cycle04_evaluation_labels_r2.json`，当前协议为 v0.2。将来修改脚本时保留旧证据，新输出另存。不要把检查通过、Git 提交或标签名写成正演或实测验证通过。
