# 三包后处理专项排查与本地修正基线

日期：2026-10-08；源码基线158ba85。用户要求先排查本地后处理，避免新增正演成本。本轮只读取既有428道native/SFCW，全部运算CPU完成；没有启动FDTD、远程任务、训练，亦未改变原包、材料、几何或执行契约。

## 结论

**确有后处理和展示层面的问题，不能将三包看不清脉络全部归为物理不可行。** 已提供独立修正基线：精确501频点、实际源相位与电流元归一化、保留复数、正确逆变换幅度、正确实带通、实际体素参考曲线、全完成道数及同一绝对色标。没有把有风险的增强操作强制串入默认结果。

原包复谱的基础计算在前轮已经独立复算通过；本轮进一步用已知延时/倾斜回波验证逆变换。**数学计算可以自洽，频窗和非线性展示仍可能遮住真实结构或使旁瓣被误认为地层。** 对实际三包，规范重绘后覆盖层底附近有明显响应，深部砂岩仍未被独立确认；本轮不宣称已经恢复完整模型。

## 查到的问题及影响

|项目|实际代码/结果|修正与判断|
|---|---|---|
|原始场预览与SFCW结果有混用风险|`plot_bscan_bipolar.py`直接取native Ez并隔8样本绘图，没有20–170MHz SFCW变换|作为原始宽带FDTD预览可以保留，但不能代表SFCW效果。新图使用2Re[复基带×起始频率载波]，不把复基带实部直接当带通信号|
|原始双极图色条标签反了|上述脚本画`-mat`配`gray`，正场为黑；色条却将−1端的黑色标为“负相位(白)”|新图明确标注零=灰、正=黑、负=白，并保留显示裁剪与数值结果的区别|
|低频子带逆变换漏除频点数|`plot_bscan_agl1.py`、`plot_bscan_compact.py`、`plot_bscan_deep_enhance.py`用IFFT×补零后长度，而官方定义是IFFT×补零倍数|20–39.8MHz实际67点，旧幅度大67倍，即36.52dB；新实现按频点数归一化。**该常数会被后续AGC抵消，不单独解释图形变化或同频带局部对比度**，但禁止据旧值跨频带比较绝对幅度|
|所谓背景抑制改变了物理量|`plot_bscan_stacked_bs.py`取包络、时间中值平滑后，除以空间中值|这是局部幅度比图，不是带符号/复数背景残差。空间共同的真实地层也变为0dB背景；不能用它判断地层是否不存在。保留原复数总场，单列风险对照|
|包络域“去多次波”不保相位|`plot_bscan_deep_enhance.py`用幅度自相关估系数，减延迟幅度后截到非负|不满足复数回波叠加关系，正负相位、多次反射与主波可能被错误处理；不能称已认证的预测反褶积。新基线不执行此操作|
|AGC、包络中值、相干加权被当作恢复证据|旧增强图对每幅结果自行归一化；有脚本硬编码“深部信号已恢复”和“对比度+4.7dB”|强旁瓣同样能被AGC抬亮，阈值截零/中值/权重会改变弱事件。新图无这些处理，采用共享绝对参考；删除由标题代替验收的解释|
|实际频窗与标题冲突|批处理默认Gaussian，但其头部写Hann；`plot_bscan_sfcw.py`同一图分别标Hann/Gaussian|以H5元数据为准。新图按实际处理标签，Gaussian与Hann并列，不能把换窗等同于地下回波变强|
|接收记录末200ns渐消|800ns包从约600ns开始，1200ns包从约1000ns开始，不是地层吸收|新基线不渐消；有限记录尾端仍有截断风险，取消渐消不等于获得更长时窗。真实底砂参考窗的总场复数变化仅约0.30%/1.04%/0.076%，不足解释该窗全部形态缺失；前轮600–800ns差异大，尾窗仍须限制|
|层位和坐标标注错误|固定7m+7.3m水平线、pkg2裁剪坐标+80而非−80、固定328ns/375ns目标窗|前轮已按实际体素重绘。本轮继续使用真实层位参考，明确95MHz局部相位路径只是近似，不是完整色散到时真值|
|未满5道尾组被丢弃|旧叠加脚本整除分组|pkg1丢4道、pkg2丢2道；新基线逐道呈现194/142/70道，不插值未完成区|

源码逐项行号和SHA-256见[源证据](../../artifacts/research_checks/2026-10-08_line9_postprocessing_audit_r1/source_findings.json)。三个包均附有这些脚本的副本，但“代码存在”不证明每张旧图都执行了所有步骤；没有伪造原批运行轨迹，也没有执行/覆写其原脚本。

## 纯数组反例：正确地层可以被后处理遮住

构造一个确定存在的倾斜界面，沿70站由290ns变化到258ns，复幅度1e−4；加入4.3ns幅度1和8ns幅度0.3的两条早时回波。仅是声明的解析信号，不是新FDTD结果，也不是实际三包的clean真值。

- 仅目标时，Gaussian和Hann均正确显示倾斜脉络；逐站最大峰位误差0.411ns，小于当前采样间隔0.832ns的一半。
- 加入早时强回波后，Gaussian在真实目标±10ns窗的早时泄漏/目标复数L2比为7.67，目标轮廓被规则横纹掩盖；Hann同一比值降到0.0529，轮廓可见。
- 已知背景的复数差分恢复该解析目标，相对L2≤4.7e−12；这里只证明信号域和运算机制，不暗示实测有这样的已知背景。
- 一个强度任意但各道相同的水平真地层，经“幅度/空间共同背景”得到0dB；均值道相减也会删除它。强制背景抑制不能作为地层存在性的检验。
- 相位反例：主波1与反相多次波−0.5叠加后的幅度为0.5，减去多次波幅度0.5再截零得到0，真实主波幅度却是1。取模后去多次波无法在不同相位下保持正确结果。

图：[已知倾斜回波反例](../../artifacts/research_checks/2026-10-08_line9_postprocessing_audit_r1/known_echo/known_inclined_echo_counterexample.png)。四面板共同幅度参考、无AGC。这直接验证“看不到模型脉络”可能发生于后处理阶段，而不是必须由错误几何或高航高造成。

## 实际三包的六组对照

逐道重建原包Gaussian、无尾渐消Gaussian、无尾渐消Hann、两端各5%余弦权重、先移除早时场再Gaussian、复频谱均值道扣除。后两项仅为有目标损伤风险的诊断，不设为生产流程。5%边缘权重为本轮探索，不称设备窗口或已验证最优配置。

所有面板共享本包无尾渐消Gaussian全记录峰值参考，禁止逐面板归一化使弱场看起来变强。真实底砂参考到时±10ns内的复幅度RMS：

|包|Gaussian无尾渐消/dB|Hann无尾渐消/dB|先平滑移除早时场、Gaussian/dB|
|---|---:|---:|---:|
|pkg1|−67.76|−103.04|−102.64|
|pkg2|−67.03|−81.13|−81.21|
|pkg3|−65.60|−83.54|−83.51|

这些是**该参考时间窗内全部成分**的幅度，不是已分离的砂岩回波。Hann和早时移除两种不同操作都使原深部背景大幅减弱；不能把“数值降低”当作目标被删掉或恢复的独立证明。保留的较浅响应更易观察，但深部仍存在归因缺口。

对5道相干叠加，底砂参考窗内“叠后幅度/叠前平均幅度”的中位数损失约0.00009/0.00476/0.00252dB，说明当前该窗**总场**不是被5道叠加整体抵消。它可能主要由相干旁瓣占据，故不能由此宣称弱真实目标已经免于叠加损伤。新基线仍取消叠加，避免改变空间采样。

## 已排除的基础实现错误与尚未排除事项

精确20–170MHz、0.3MHz、501点，源半时间步及接收时间原点、复谱除实际源谱、I/Q实虚部、载波回补、IFFT符号与补零幅度，均有原始审计及本轮已知延时验证支持。逆变换采样为0.83167ns，周期1/Δf=3333.33ns；不是原始800/1200ns记录被延长，也不是补零提高物理带宽。不能把native0.05897ns、重建0.83167ns与旧实测501样本混用。

本轮已知回波直接复数求和与IFFT相对L2≤9.0e−12；补零前后共同样点相对L2≤7.4e−16。结果支持基础公式正确，但不构成模型收敛或实测仪器链验证。材料、FP32弱场底、二维/天线差别、边界等问题本轮不重跑、不宣布排除。

## 本地后续使用规则

先看规范的全频带总场：Gaussian/Hann并列、正确双极与复包络并列、同一幅度参考、无默认AGC/叠道/包络背景除法/包络去多次波；保留频谱与复数结果，再决定是否添加背景抑制或增益。未经目标保真检查的增强图只作诊断。

本轮私有复数缓存留在忽略目录 `artifacts/local_checks/2026-10-08_line9_postprocessing_cache_r1/`，可继续CPU处理，不需要重新FDTD。与实际三包相关的目标身份尚未确认；接下来的本地研究仍可比较分段时间场、频谱边缘和保相位背景估计，不能据一张更漂亮的图宣称找到砂岩。

## 图和复现入口

- pkg1：[六组处理对照](../../artifacts/research_checks/2026-10-08_line9_postprocessing_audit_r1/line9_pkg1_v5_d02m_194st_processing_ablation.png)、[正确双极SFCW](../../artifacts/research_checks/2026-10-08_line9_postprocessing_audit_r1/line9_pkg1_v5_d02m_194st_signed_sfcw.png)
- pkg2：[六组处理对照](../../artifacts/research_checks/2026-10-08_line9_postprocessing_audit_r1/line9_pkg2_x160-110_agl35_142st_processing_ablation.png)、[正确双极SFCW](../../artifacts/research_checks/2026-10-08_line9_postprocessing_audit_r1/line9_pkg2_x160-110_agl35_142st_signed_sfcw.png)
- pkg3：[六组处理对照](../../artifacts/research_checks/2026-10-08_line9_postprocessing_audit_r1/line9_pkg3_x160-120_agl1_70st_processing_ablation.png)、[正确双极SFCW](../../artifacts/research_checks/2026-10-08_line9_postprocessing_audit_r1/line9_pkg3_x160-120_agl1_70st_signed_sfcw.png)

本批图件位于仓库相对路径 `artifacts/research_checks/2026-10-08_line9_postprocessing_audit_r1/`，JSON含脚本、原审计、native/SFCW来源核验及CPU检查。原包逐文件身份仍由前序审计绑定，未上传原始资料。

```powershell
D:\gprmax_v4_gpu_env\Scripts\python.exe scripts/audit_line9_postprocessing.py --root artifacts/local_checks/2026-10-07_line9_three_packages_review_r1 --review artifacts/research_checks/2026-10-07_line9_three_packages_review_r3 --out artifacts/local_checks/postprocessing_reproduce_report --cache artifacts/local_checks/postprocessing_reproduce_complex
D:\gprmax_v4_gpu_env\Scripts\python.exe scripts/check_line9_postprocessing_known_echo.py --out artifacts/local_checks/postprocessing_known_echo_reproduce
```

依赖本机V4 SFCW processing、NumPy、h5py、Matplotlib及中文字体，输出路径必须不存在。没有求解器调用。本轮只做针对处理变更的CPU验证，不重新跑无关的研究数组或仿真批次。
