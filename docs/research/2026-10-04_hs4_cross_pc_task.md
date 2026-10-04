# HS4跨电脑完整接续任务

用户已授权自主推进并要求“若不行，就把完整任务git上去，我另一个电脑去做”。本任务无需聊天上下文。入口是本页、[机器任务及20份输入](../../artifacts/research_checks/2026-10-04_hs4_cross_pc_design_r1/task.json)、[12份已完成细网格控制目录](../../artifacts/research_checks/2026-10-04_hs4_cross_pc_design_r1/completed_cases.json)。用户的自主GPU研究授权持续有效；旧门禁里的历史未批准描述不撤销后续授权，但每个新批次仍先冻结具体契约。

## 目标与已经完成的工作

解释15m航高时，为什么HS4T2D的B-scan上部地下条带没有直接复制起伏界面外形；分清几何/处理错误、数值边界或网格影响，以及真实传播与多位置相干叠加。正确性优先于把图画成模型形状。先完成数值稳健性，再核查三维对应模型。此任务不扩大处理产品、启动训练或使用留作后期验证的营山/雅安数据。

此前[47道连续场及局部扰动](2026-10-04_hs4_height_wavefield_findings.md)完成。高航高的共同较浅拱顶解释取得区域因果证据；固定几何PML对照不支持当前带限上部条带主要由侧边反弹造成。中心三个网格的区域关系稳定，但精确波形没有绝对误差认证。

本次再完成左右端1.25cm八道，原生V4.0.0/CUDA double/Ey float64。输入、实际收发格点、CFL、有限值、八份全材料图以及与既有中心同变体几何一致均独立核查PASS。每道约154–169秒，实际监督所有子进程，峰值RSS约2.04–2.05GB；与旧35道只量启动器的历史限制分开。第一轮实时容量预检拒绝发生在创建执行事件/求解之前，不消耗attempt；资源恢复后按原冻结契约执行一次。

| 固定诊断 | 左端 | 中心既有结果 | 右端 |
|---|---:|---:|---:|
| 1.25cm早期160–180ns拱顶/坡面变化复包络L2比 | 58.9035 | 7.556（既有） | 4.1751 |
| 1.25cm后期180–220ns拱顶/坡面变化L2比 | 0.4317 | 见既有结果 | 0.3807 |
| 2.5cm→1.25cm基准地下160–220ns带符号相对差 | 27.1598% | 不是此相邻对，中心已有中间档 | 32.0962% |

上部波列在左右端仍更敏感于拱顶，后续波列更敏感于坡面；这些是等面积扰动引起的波形变化比，不能解释为各区域回波能量百分比。端点相位仍明显依赖网格，不能用中心的7.586%相邻差代表全扫描。新[数组、结果及对比图](../../artifacts/research_checks/2026-10-04_hs4_station_grid_results/summary.json)从原始H5重新构建，数组/summary/图三文件均字节一致。旧47与新增8共55道已消耗，勿重跑历史契约。

## 尚未运行的准确任务

以下20份输入已生成并通过CPU几何/资源门禁检查，**不是20份正演结果**。先按阶段执行，不一次排全扫描。

| 阶段参数 | 道数 | 用途 | 启动前可用RAM / 空闲VRAM门槛 |
|---|---:|---|---|
| `replay` | 4 | 在目标机器另冻中心1.25cm四控制，记录跨机器/构建差异 | 2.25 / 3.5 GiB |
| `centre1cm` | 4 | 原界面、全覆盖层、拱顶、坡面；补第四档原生网格 | 3.75 / 4.3 GiB |
| `3d-centre` | 3 | COL6实际中心站位起伏0.8m剖面、平面、全覆盖层 | 8 / 5 GiB |
| `3d-ends` | 6 | 左右端同三种模型，中心检查后才扩展 | 8 / 5 GiB |
| `3d-xwide-centre` | 3 | 同中心位置、全原内部几何平移后，检验近X边界 | 10 / 6.5 GiB |

门槛是保守启动预算，不是已实测三维峰值或物理标准。三维240×240×660格的字段/ID/单极Debye ADE主数组设备下界约3.432GiB，X扩宽320×240×660约4.572GiB，**尚未含PML、上下文、编译/初始化临时空间**。主机另有solid/rigid、PML与初始化峰值。当前6GiB显卡/16GiB主机的可用RAM不满足三维门槛；建议另一台至少32GiB主机、8GiB显卡，并按实际空闲量预检。不要通过减少域、改精度、改变原界面或停用户程序绕过预算。

三维注意事项：

- 旧三维实际COL6 X=1.60m，不是输入名义1.625m；本任务将输入显式落在实际格点。中心Tx/Rx为[1.6,5.6,27]/[1.6,6.9,27]m，地表Z=12m，航高15m。
- 对旧48×48界面统一用 `Z_new=9.05+2*(Z_old−9.05)`。第6个X bin对应二维已执行剖面，48点逐项相同、该剖面起伏0.8m。**整个三维面的起伏为1.6m，Z范围8.15–9.75m**；不能把剖面0.8m称为全三维面的幅度。
- 新平面底界Z=9.05m（覆层2.95m），匹配当前放大后的参考面；旧HS1 Z=9.0m不代替本参考。全覆盖层差分包含想保留的界面，不自动成为clean训练标签。
- 近X边界控制将域12m扩至16m，整个旧内部箱体与收发X均+2m，新增左右各2m区域延续原最外行界面。原内部几何逐箱体不变，PML仍1m，Tx到内侧PML边缘由0.6m增至2.6m。不要把扩域同时重塑内部地形的结果当纯边界检验。
- 这里三维是理想X电偶极/Ex接收，基线沿扫描Y，属于继承的inline机制模型。项目仪器已确认cross-track基线，故当前模型**尚不认证设备几何**。二维是旋转后的Ey线源机制控制，也不等于真实有限天线。需单列后续姿态/有限天线试验，不能混进本批后声称纯维度差异。

## 另一台电脑可执行步骤

先同步干净main，再开本机独立`codex/`分支；若本地有其他改动，不覆盖。Python数组环境与V4环境分开。

```powershell
git pull --ff-only
git switch -c codex/hs4-cross-pc-controls
```

使用用户原先提供的同一V4源码/本地wheel建立Python3.11–3.13独立GPU环境。安装、编译方法见[V4构建记录](2026-09-24_v4_build_success.md)、[GPU环境记录](2026-09-24_gpu_environment.md)和两个依赖清单；本次[实际依赖版本](../../artifacts/research_checks/2026-10-04_hs4_cross_pc_design_r1/observed_dependencies.json)另存。前述构建文档的旧“未授权/下一步”仅为当时状态。Git不包含V4发行包、venv、wheel或CUDA/MSVC；在另一台已有的V4上先预检，不因为目录叫gprmax就使用V3.1.6。若尚无发行包，取得用户已提供且源码指纹相同的副本，不随意以同版本号的另一源码替换。不同编译二进制可以新冻，但必须做下面的目标机replay，跨构建相位差没有误差预算时明确未决，不自动继承旧精确波形认证。

下列变量只需设成本机路径。所有设计输入使用相对路径，历史冻结契约保持原字节，**不手改旧契约里的D盘路径**。

```powershell
$env:HS4_PYTHON = 'C:\your-v4-env\Scripts\python.exe'
$env:HS4_VCVARS = 'C:\your-Visual-Studio\VC\Auxiliary\Build\vcvars64.bat'
$env:HS4_CUDA_BIN = 'C:\your-CUDA\bin'
$hs4Design = 'artifacts/research_checks/2026-10-04_hs4_cross_pc_design_r1'
& $env:HS4_PYTHON scripts/verify_hs4_cross_pc_archive.py --catalog "$hs4Design/completed_cases.json" --out artifacts/local_checks/hs4_archive_target.json
& $env:HS4_PYTHON scripts/check_hs4_cross_pc.py --design $hs4Design --out artifacts/local_checks/hs4_cpu_target.json
.\scripts\run_hs4_cross_pc.cmd preflight --out artifacts/local_checks/hs4_runtime_target.json
```

`preflight`只导入、核对源码/编译库身份并检查资源/编译器，不跑官方示例或求解。只写新的输出路径，目录需先存在。`check_hs4_cross_pc.py`当前10项均已通过；其四道原始资料的CPU复核不叫新电脑replay。

先另建目标机器replay；这个阶段用于跨机器比较，是一个新批准的四道契约，不是重用历史attempt：

```powershell
$hs4Replay = 'artifacts/research_checks/target_hs4_replay_r1'
.\scripts\run_hs4_cross_pc.cmd freeze --design $hs4Design --stage replay --out $hs4Replay
.\scripts\run_hs4_cross_pc.cmd run --out $hs4Replay --execute
& $env:HS4_PYTHON scripts/analyze_hs4_cross_pc.py --capsule $hs4Replay --out artifacts/research_checks/target_hs4_replay_results_r1
```

之后一次冻结并运行一个阶段。先`centre1cm`，再`3d-centre`；中心核查后才扩展`3d-ends`及`3d-xwide-centre`。

```powershell
$hs4Stage = 'centre1cm'
$hs4Capsule = "artifacts/research_checks/target_hs4_${hs4Stage}_r1"
.\scripts\run_hs4_cross_pc.cmd freeze --design $hs4Design --stage $hs4Stage --replay-report artifacts/research_checks/target_hs4_replay_results_r1/summary.json --out $hs4Capsule
.\scripts\run_hs4_cross_pc.cmd run --out $hs4Capsule --execute
& $env:HS4_PYTHON scripts/analyze_hs4_cross_pc.py --capsule $hs4Capsule --out "artifacts/research_checks/target_hs4_${hs4Stage}_results_r1"
```

换阶段时修改`$hs4Stage`，使用新的目录。三维扩X结果分析额外传 `--reference-capsule artifacts/research_checks/target_hs4_3d-centre_r1`，报告相同位置粗糙/平面/全覆层原始复谱变化；比较双方必须核对同一构建。所有命令成功/失败如实留档。

监督器禁止覆盖已有胶囊、重复已创建execution日志的attempt、自动重试失败求解和CPU回退。容量拒绝在首次创建执行日志前不消耗attempt，可以资源恢复后使用同一未启动胶囊；**只要任何组已STARTED/有execution日志，就保留原胶囊，分析原因，另冻只需补做的具体任务，不重跑已完成组**。当前脚本阶段固定分组；若部分失败需新的更小分组，由接续agent按原冻结输入和已完成目录创建并审核新的有界契约，不直接删除日志/原始H5使入口放行。仅终止本任务自己的子进程树；不关其他用户软件。单道/整批时限7200秒/7200×道数，实时RSS及最低系统余量有监督。

所有运行明确 `-gpu 0 -gpu_precision double`，检查GPU0，不要通过另一个`CUDA_VISIBLE_DEVICES`映射把容量预检与求解卡分开。实际H5必须原生float64，不能用转dtype代替double。Windows运行包装脚本只根据本机三个变量配置工具链，历史硬编码包装器不用在另一台执行。

## 评价与交付标准

保持实际501个20–170MHz含端点、0.3MHz频点；源归一化、原生200ns尾渐消、单位均值Hann/8倍补零、固定160–180/180–220/160–220ns。差分复谱后才取模。输出复谱、带符号波形、复包络和可比较尺度图，记录源/接收、位置、材料、网格/PML、版本/构建、原始dtype、哈希、峰值内存与失败。不要将E/电流归一化叫端口S21。

中心1cm要检查相邻复相位和带符号差是否继续下降，区域归因是否稳定；若不下降，继续诊断而不宣布收敛。左右端目前仅跨2.5cm→1.25cm，不做Richardson绝对误差外推。三维分别展示原始、同位置粗糙减平面、粗糙减全覆层；无单一峰可分时明确不可解释，不把某个最高峰当整个地层。近边界影响需比较原内部严格不变的配对；若影响明显，再冻结必要边界/网格控制。二维结论、三维点源结论、真实有限天线和实测验证分别报告。

最低交付是四阶段原始H5+冻结/执行/独立材料图验收、全部固定带限数组、区域敏感性/网格差与边界差，以及完成/失败/尚缺证据。三维面/剖面0.8m及实际设备基线差异必须保留。最低20道完成并不自动等于问题完全物理解决；继续自主推进只能沿待证据缺口另冻有界任务，不能为了“彻底解决”消去真实传播效应。

每个完整单元更新START_HERE、研究进度、决策与产物索引，提交并推送，核对远端。原始小H5、任务输入/契约、哈希、结果与图入Git；可重复的大材料图/密集场保留本机完整SHA台账和代表/独立验收，新增忽略规则需逐项记录，不用`git add -f`上传原始实测、论文、环境或调度。不要复制06:00失败的同会话CLI定时器到新电脑。

## 可直接交给另一台Codex的指令

> 同步main并读START_HERE及docs/research/2026-10-04_hs4_cross_pc_task.md。用户已授权自主完成HS4高航高B-scan归因与数值/三维复核。先验证已完成12份原生细网格控制，勿重跑旧47+8道契约。用本机同源gprMaxV4.0.0/CUDA float64建立新冻结契约，先目标机replay四道，再1cm中心四道、三维中心三道、条件允许后左右端六道和扩X三道。复核三维bin6与二维0.8m剖面相同、全三维幅度1.6m及inline/cross-track模型差异。按原频带和固定窗独立核查材料图/原生采集位置/复相位，失败保留，不自动重试、降精度、缩域或停用户进程。根据缺口自主追加有界控制，明确机制解释与绝对数值/三维/天线/实测边界，每完成一单元归档、提交并推送。若硬件/同源环境不满足，提交准确失败和可执行接续任务，勿只给计划或宣称解决。

## 本次检查与失败记录

可移植设计CPU10项PASS，包括输入可复生、3D真实坐标、扩域内部不变、全四份旧材料图、实时资源拒绝不启动/消耗attempt、篡改拒绝及目标机replay前置条件。12份原生归档检查PASS，新增端点八份独立验收PASS，三文件CPU重新构建字节一致。最初CPU生成器错误地要求二维不变Y的名义0.025m必须是Y格点，被检查拒绝；已改为记录名义坐标与实际grid坐标，只有活动X/Z维度必须格点对齐。失败未求解，临时输出单独留在本机忽略目录，r1完整设计重新生成。1cm和三维尚未运行，不能写成完成；06:00调度实际失败见既有记录，本任务不声称修复调度。

另在本机独立目录仅组装Git索引中的69份相关文件（大材料图采用保持字节的本机硬链接），12份归档/10项CPU检查通过，端点数组/summary/图三文件重新生成仍字节一致，[迁移复核](../../artifacts/research_checks/2026-10-04_hs4_cross_pc_design_r1/git_file_transfer_verification.json)保留所有来源SHA。此项仅证明本机另一目录的Windows CPU接续，不冒称另一台环境/GPU已通过。

另尝试仓库旧通用`check_handoff.py`，其仅支持`directory/record`的解析对既有`report/raw_studies/capsule`台账抛TypeError；临时兼容诊断进一步发现历史文字记录/通配占位路径和旧单道门禁与当前批次格式不符。该通用检查**未通过**，不是本次专用核验的PASS项；未改写历史记录或冻结门禁让它通过，临时兼容改动撤回。新机按上面的专用入口核验准确12份科学文件及20份设计，不依赖该旧检查。通用检查迁移另属仓库维护问题；main上未发现可用GitHub Actions运行，不能声称CI通过。
