# 八站密采样前缀与四组波场控制执行接续

## 新增证据与研究边界

持续目标仍是解释模型、B-scan与实测差异。本单元没有把目标缩成完成5m段或制作快照脚本，也没有认定低损耗材料真实。上一目标回合核查了确实存活的46项任务；本单元继续同一进程树，取回16份新增native及2份复用，独立复算后交付六张中文图，并将后续四组输入和有界执行/输出审核代码放到ROG。没有新GPU求解、快照、正式材料替换、三维或训练。

公共前缀为`artifacts/research_checks/2026-10-09_v401_dense_loss_partial_r4/`。对应私有快照`artifacts/local_checks/2026-10-08_v401_dense_loss_results_r1/snapshots/1791478761367623200/`，取回ZIP SHA256 `be5f4172fa7970ae51e452d672f2317ee32cfff1e4c054228283a0244f5dcd7e`。共有190→188.25m八站两种损耗H1，只有190m具备双材料H0；未算站位仍掩码，不插值。16/46是这个不可变快照的完成数，不是实时状态。最新另行SSH核查已完成20/46，low_x18750_H0运行，两个owned Python进程存在；GPU1%/841MiB为该项启动阶段瞬时读数，不能据此认为任务停止。接续时必须重新poll。

八站低损耗Hann/Blackman总场峰均随局部底砂模板变浅，每站峰位比各自预测提前0.83167ns。高损耗Hann的预测残差范围−7.48503至+5.82169ns，Blackman为−6.65336至+5.82169ns，峰位继续跳变。这加强了“当前高损耗假设压弱底砂、其他响应叠加改变最强峰选取”的解释，但其余七站没有H0，不能逐站认证高损耗峰的物理身份，更不能称低损耗配方为场地标定。

独立18份native/元数据/完整501点DFT/逆求和/缺列掩码通过，最大DFT相对差约1.06e−13。另190m共同窗的24项标量独立核对通过：底砂差场低/高仍724.029倍（Hann）/746.791倍（Blackman），H0绝对范数低/高仍0.97647/0.93329。剩余响应没有消失，这些范数比也不是能量占比或实测SNR。图件已人工查看中文标签、色标与缺列显示。

图件：

- `artifacts/research_checks/2026-10-09_v401_dense_loss_partial_r4/dense_loss_hann_main_gray.png`
- `artifacts/research_checks/2026-10-09_v401_dense_loss_partial_r4/dense_loss_blackman_main_gray.png`
- `artifacts/research_checks/2026-10-09_v401_dense_loss_partial_r4/dense_loss_hann_late_gray.png`
- `artifacts/research_checks/2026-10-09_v401_dense_loss_partial_r4/dense_loss_blackman_late_gray.png`
- `artifacts/research_checks/2026-10-09_v401_dense_loss_partial_r4/basal_zoom/dense_loss_hann_basal_zoom.png`
- `artifacts/research_checks/2026-10-09_v401_dense_loss_partial_r4/basal_zoom/dense_loss_blackman_basal_zoom.png`

同样的合同/清单/批前审核通过`identity_references.json`按原路径与SHA引用，未重写旧证据。`delivery_manifest.json`逐文件记录本次交付SHA。原生、完整复谱和原始实测资料仍留忽略目录。

## 波场执行代码与实际两机检查

下一批仍为[非局部覆盖层候选与四项控制](2026-10-09_nonlocal_cover_candidate_and_controls.md)中的原几何/全平层快照，以及远/近覆盖层底反差移除。输入没有改变。新增代码：

- `scripts/line9_v401_cover_wavefield_run.py`：freeze/run/verify/collect；沿用既有有界监督器，只能4项一次attempt。
- `scripts/check_line9_v401_cover_wavefield_execution.py`：不启动求解器的内存选择、元数据、终止门禁和复数差分检查。
- `scripts/analyze_line9_v401_cover_wavefield.py`：精确501点、同源尺度、Hann/Blackman，六张单站配置图与固定窗指标。
- `scripts/audit_line9_v401_cover_analysis.py`：独立native DFT及直接逆求和重算指标。
- `scripts/render_line9_v401_cover_wavefield.py`：原/平层同一固定灰度、差场单独固定灰度的Ricker GIF，按原生细网格提取真实地表/覆盖层底曲线。

当前安装的4.0.1 `Snapshot.nbytes`在场数组初始化前仍为0，内存检查发生在此前。已有`gprmax_snapshot_cuda_entry.py`只在官方内存预检期间补六分量输出字节数、随后恢复，官方函数选择GPU→主机存储；没有改场更新/快照内核。本次实际调用两机安装的官方内存选择函数，16GiB容量夹具、250帧六分量历史1,323,000,000字节，能启用官方存储；非快照模型本身超过显存时拒绝并恢复计数。超限夹具日志中的容量警告是预期反例，不是本批实际显存不足。

两机各17项CPU检查通过，代码与官方内存选择源SHA一致。检查还覆盖快照版本、六分量、float64/有限值、原点与电/磁半步时间；46个独立收据/原生哈希一致要求；重复/缺失收据、未终止或仍存活的前批进程拒绝；四个固定窗的复数幅度与反相差分。这些使用小HDF5及明确标为fixture的执行证书，不代表实际GPU分配、观察器不扰动接收记录、四组求解完成或传播机制认证。

首份ROG准备检查因`unittest.mock`间接导入asyncio、Windows应用控制策略拦截`_overlapped` DLL而失败，发生在任何新契约/求解之前。仅将检查脚本的替身替换改为作用域context manager，未更改系统策略、解释器或求解环境。失败准备目录`E:\line9_basal_pair_20261008_r1\v401_cover_wavefield_r1`保留；可接续的是全新r2目录。开发中间检查保留于本机忽略目录，不以过时源码哈希冒充当前代码验证。

实际部署：`E:\line9_basal_pair_20261008_r1\v401_cover_wavefield_r2`，32份payload文件逐字节核验，ZIP 334609字节、SHA256 `b9e1ae2c3e2cce29e60dfdeaebd89195237c9a3b70830be868a7737b91f3d92f`。manifest仍`e3e4b8297c325fb4816e2d2a6a1c7505d713983acb829415a18417003e6955b2`。ROG输入审核通过；实际当前46项处于运行状态时，前批门禁拒绝，且不存在新execution目录。

证据：`artifacts/research_checks/2026-10-09_v401_cover_execution_cpu_r3/checks.json`与`artifacts/research_checks/2026-10-09_v401_cover_execution_target_r1/`。这里`execution_ready=false`仍然正确，未冻结/未启动四项；分析和动画代码仅已做语法检查，其真实完整输出审核待执行后才能成立。

## 执行顺序和复现入口

先用已有`artifacts/local_checks/2026-10-08_v401_dense_loss_deploy_r1/poll.ps1`续600s租约并核验进程。完成后用同目录`retrieve.ps1`取完整46份，运行既有dense分析/独立审核/共同窗审核，特别检查187.5m与185m新H0锚点。不得把观察超时当作停止或重启旧attempt。

四项的新freeze强制要求旧46项最终COMPLETED事件、46份唯一收据、完整native审计与证书哈希一致，以及前批owned进程已退出；同时核对前批manifest与190m原H0哈希。之后在目标机重新核验空闲RAM/VRAM/磁盘和实际运行时/编译器。预算为快照历史加6GiB设备预留约7.232GiB，主机空闲20GiB加一份历史约21.232GiB，磁盘两份历史加2GiB。都是预检预算，实际峰值仍待运行。每项1800s/全批9000s、600s租约/USER_STOP/只停owned进程/无重试。

目标机r2已有`freeze.cmd`和`run.cmd`（CUDA13.3/MSVC及已核验4.0.1解释器）。须等旧46项完成审核、退出后依次运行，不能现在执行：

```powershell
& 'E:\line9_basal_pair_20261008_r1\v401_cover_wavefield_r2\freeze.cmd'
# 检查execution_contract/preflight，记录SHA，再启动；运行期间仍须续session_heartbeat租约。
& 'E:\line9_basal_pair_20261008_r1\v401_cover_wavefield_r2\run.cmd'
```

运行后在同目录用其解释器及`repo/scripts/line9_v401_cover_wavefield_run.py verify --out execution`核验全部500帧六分量、电/磁时间、采样与原接收/source字节一致，再用`collect --out execution --export-out <全新目录>`导出小native及收据；完整约2.464GiB场史留目标机。`render_line9_v401_cover_wavefield.py --execution execution --out <全新图目录>`在目标机渲染；图、GIF与审核/哈希一起取回，不能将远端收据说成本地重新读取全场。

`analyze_line9_v401_cover_wavefield.py --source <取回目录> --out <全新公开分析目录> --numerical <全新私有H5>`后接`audit_line9_v401_cover_analysis.py --source ... --public ... --numerical ... --out <全新审核JSON>`。预声明300–450ns、原底砂345.618097–374.608117ns、0–120ns及450–1100ns窗口保持；四项都是去深部底砂的H0，仍保留覆盖层真实界面，不能把所有H0叫噪声。区域替换增加边缘、改变下部延续与全部交互；平层改变整个横向几何，差场都不是纯多次波。

当前约330ns强波仍只具有非局部覆盖层候选证据。后续按四组控制和真实波场核查，不用射线接近到时替代物理归因；真实天线、有限三维及实测响应差异仍需继续，研究目标保持active。
