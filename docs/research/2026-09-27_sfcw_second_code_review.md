# SFCW统一适配第二次只读代码审查

**范围：**审查[适配脚本](../../scripts/research_sfcw_unified_adapter.py)、修正后manifest/NPZ、[root独立验收代码与结果](../../artifacts/research_checks/2026-09-27_sfcw_unified_dev_adapter_mt33_rerun/)、[采集与SFCW数据契约](2026-09-27_acquisition_and_sfcw_data_contract.md)及本机官方V4 `processing.py`。未读测试组，未运行整批、FDTD或任何求解器，未改代码与产物。

## 结论

修正后的正式产物没有发现会推翻当前数值结果的代码或数据形状错误。MT数据确为指定C3 BG/TGT `-MT33` H5，每组按H5 `rx1`至`rx33`数值顺序调用官方`name:mt01`至`name:mt33`后堆叠；只读核对其Name和Position属性与manifest次序一致。root_verify独立从原始H5按`rx1..rx33`重建，并验证频响shape、坐标、CO t05/MT rx17锚点。CO为22个单道H5，正确使用`name:measurement`。首轮误选母模型单道H5的结果在独立归档中，并未混入正式rerun。

| 检查 | 证据与审查结果 |
|---|---|
| 输入边界、道序与shape | 适配器输入表固定为2个C3 MT33 H5与22个C3 CO11 H5；MT自然数字段序，堆为`(nt,33)`，CO存储为一维`(nt,)`。正式manifest列24个源H5 SHA、receiver names/positions；root_verify逐H5 SHA回读、核对H5属性位置，并在9个频点上比较MT全部33道和CO单道，结果shape符合501×33/一维501点。未发现C5/C8路径。 |
| 时间轴与时间原点 | V4 `load_source/load_receiver`读取H5采样间隔与`TimeSampleOffset`；MT stack前断言33个receiver的dt/time offset相同，源dt与接收dt一致。源样本原点为`dt/2`、Ex接收原点为0；NPZ含两条物理样本时间轴，root_verify逐元素核对`offset+n·dt`，并以原点相位做显式复指数和。未见半采样偏移丢失。 |
| 频响与单位 | `direct_frequency_response`使用官方engineer DFT/CZT计算接收谱除以源谱；H5源历史的时间偏移计入频谱。记录了source `electric_current`、A单位及`SpatialScale=0.025`，接收Ex为V/m；脚本不重复乘/除SpatialScale，响应单位按接收样本/源样本。它是源历史归一化的场响应，不是端口S21、校准天线电压或实机接收链。root_verify独立用显式复指数求和核验1584个复数值（9/501个频点覆盖全部trace/records），最大误差相对所选峰8.87e−12。 |
| 重建、零填充与边界 | 官方重建为矩形窗、pad=4、shift=0；保存source-valid、权重、复包络、complex bandpass及real bandpass。2004点采样间隔约1.663 ns、等效采样率601.2 MHz，高于170 MHz载波的Nyquist要求；3.333 μs包络周期由df决定。填零只插值，不增带宽/分辨率；超过1.2 μs原始记录的包络部分是有限频带的周期重建，不是新传播观测。root_verify按显式DFT和9个时点的逆和核验，最大相对误差2.16e−13。 |
| 尾部与窗 | 主组未加尾窗且保存全长度；200 ns组按契约`(round(200ns/dt)-0.25)/N`，与官方raised-cosine实现一致。末5%样本峰/全记录峰范围−46.47至−32.82 dB；200 ns窗相对无窗谱L2变化MT约3.376%、CO约3.474%。`record_decay_certified=false`正确，taper没有证明收敛。当前MT tail值是跨33道的聚合最大值，不是逐道衰减认证；原始H5可支持后续逐道核对。 |
| 求解器边界 | 适配脚本只导入`gprMax.toolboxes.SFCW.processing`及NumPy/H5读取，实际调用load、频响和重建函数；无solver、输入生成、子进程或GPU调用。root_verify只读H5并做显式DFT和矩阵求和。当前处理不是隐式FDTD。 |

## 可复现边界与改进项

1. 当前验证已足以支持这批开发数据的处理链，但root_verify的独立DFT只抽查9/501个频点，逆变换只抽查9/2004个时间点；它不是全数组独立重算。若以后要扩大数值验收，可对一组MT33及CO代表输入做完整501频点、2004时点显式求和，或明确把当前结果称为抽样独立核验。
2. manifest没有显式记录`direct_frequency_response`的`source_floor_db`实参；脚本使用本机冻结源码的默认−100 dB，且保存了501点source-valid掩码、源码SHA和gprMax版本。当前H5在本次均通过。未来复现若要脱离默认值，建议调用时显式传−100并写入manifest；这属于元数据完备性，不是已发现的数组错误。
3. MT每个记录只保存跨道聚合的尾部dB标量。当前未窗原始H5完整可读，故没有隐藏尾部；如需逐道讨论尾部，直接对这些同一输入按每道计算，勿将现有聚合量解释为每道收敛证据。

**旧首跑处理：**不带`-MT33`的MT母模型误选产物已隔离为`artifacts/local_checks/2026-09-27_sfcw_adapter_wrong_mt_input/`，保留但不属于正式产品。当前正式manifest SHA `00ecdf7c62fabd62b45af7c43aaed125161f9456ea95c20d077c933bb3c6125a`，NPZ SHA `52f544f8516dceb9a8573e4f1f35efbb78f6b971d7c315526fb95df4f6a30984`。本地源码审查不替代主控正在做的官方文档/空间分辨率与材料物理核查；这些限制决定解释边界，不推翻本批代码数值验收。
