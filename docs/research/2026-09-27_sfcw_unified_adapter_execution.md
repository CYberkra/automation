# C3开发数据官方SFCW统一适配执行记录

本次仅对已有开发侧C3二维H5做CPU后处理，没有运行FDTD、没有读取C5/C8、没有训练，也未改评价/阈值文件。主入口为[research_sfcw_unified_adapter.py](../../scripts/research_sfcw_unified_adapter.py)；输出在[统一适配产物目录](../../artifacts/research_checks/2026-09-27_sfcw_unified_dev_adapter_mt33_rerun/)。

脚本对CO BG/TGT各11个单道H5，以及MT33 BG/TGT两个H5中的mt01–mt33，调用已安装gprMax V4官方`load_source`、`load_receiver`、`direct_frequency_response`和`reconstruct_time_response`。共24个H5，输出48组（每个输入无taper及200 ns tail-taper敏感性各一组）；MT频响数组为501×33，CO单道存储为一维501点（11道可按接收位置堆叠为501×11）。结果NPZ包含复频响、源/接收谱、source-valid、官方窗权重、复包络/带通信号、real-bandpass及物理时间轴；逐条input SHA、源与接收采样偏移/坐标、源类型、单位与空间尺度、taper、重建选项和官方源码SHA在manifest中。

频率为20–170 MHz等间隔501点。重建固定使用rectangular窗、zero-pad=4、shift=0；填零仅把显示采样率插值至601.2 MHz，以表达170 MHz载波，不增加带宽或分辨率。complex envelope周期为1/df≈3.333 μs；保存它作为首选表示。real-bandpass是调制后的参考量，端点周期连续性不作假设。两种记录处理均保留全长原始接收历史，不做时间截短；官方尾部诊断统计记录末5%样本峰值相对全记录峰值，范围为−46.47至−32.82 dB。**记录衰减/收敛未认证**；200 ns taper只是敏感性对照，不能隐藏或证明尾部收敛。

合成自检只调用官方变换：源半采样偏移相位最大误差3.72×10⁻¹⁶；构造延迟峰落在预测索引92；双道输出shape为501×2，线性配对最大误差1.78×10⁻¹⁵。它验证API对齐与变换链，不验证物理精度。主控独立验收从原始H5直接计算了1,584个复数DFT值，最大误差/所选峰值为8.87×10⁻¹²；对每道9个时点直接求逆变换，最大误差为2.16×10⁻¹³；CO t05与MT rx17四个BG/TGT锚点通过。200ns taper相对无窗谱L2变化约MT 3.376%、CO 3.474%。详见产物目录中的`root_acceptance.json`。输入是二维x-invariant机制模型；H5中x=0为占位坐标，不能视作实机位置。本批只作为开发侧SFCW处理代理，不构成clean实机数据或物理验收。

正确产物SHA：manifest `00ecdf7c62fabd62b45af7c43aaed125161f9456ea95c20d077c933bb3c6125a`；NPZ `52f544f8516dceb9a8573e4f1f35efbb78f6b971d7c315526fb95df4f6a30984`。官方`processing.py` SHA为`adad556f09140956f0ee19d3038430e06a9ae8be6a826dcad723096d99624a3b`；适配脚本SHA为`4a8a6a4ed319f8ca6bd57f7a99d9c8d60b95fec64f2a0ca273e317f7ad6f2088`。复现命令：`artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe scripts/research_sfcw_unified_adapter.py`。脚本拒绝覆盖已存在输出目录；主agent实测拒绝覆盖通过。现有产物无需重导：只验收可运行 `python artifacts/research_checks/2026-09-27_sfcw_unified_dev_adapter_mt33_rerun/root_verify.py`。如需独立复现，给适配脚本传 `--output-dir artifacts/local_checks/sfcw_reproduction_new`（目录须不存在）。另一台机器使用已安装gprMax V4、NumPy、SciPy、h5py的Python；版本和源码SHA见manifest，环境变化不自动视为等价。

**输入选择失误留档：**首次执行误用了不带`-MT33`的两个母模型单道H5，故其MT结果无效、不纳入上述计数及解释。原目录保留在`artifacts/local_checks/2026-09-27_sfcw_adapter_wrong_mt_input/`；旧manifest SHA `24a1c75682986ad547545176993102adb9717c7f662a34503f5bf74e7072892e`，旧NPZ SHA `25cc28657ab2987b50bedd7bf679209836e5f59045c2ca442385efe7f50b9bdf`。更正后的脚本默认指向独立新目录，没有覆盖首跑结果。
