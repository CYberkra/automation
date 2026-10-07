# 三包后处理修复交付与地层对应复核

日期：2026-10-08；源码基线5f27370。用户要求继续解决问题，优先本地排查。本轮全为CPU离线处理，没有FDTD、远程任务、训练、材料/几何修改。

**总数更正：194＋142＋70＝406道。此前报告和入口写成428，是合计笔误，不是另有22道、漏处理或本轮新增道数。逐包原始清单、旧数值JSON及其哈希不改；只更正文字合计和索引元数据。**

## 现在解决了什么

已提供直接从native重建的修正SFCW入口和三份可读取数值HDF5，不再依赖原包中值比、包络相减、AGC或硬编码层位的绘图流程。最明确的新结果是：**1m包原始数据中存在随实际底部砂岩走时变化的弱回波；原Gaussian背景遮挡和旧参考线使其难以辨认。** 它在Hann和Blackman中均与完整材料谱的一次反射近似高度对应，不能再说这个包完全没有模型脉络。

这不是三个包已经全部物理验收。8m/v5和固定y35m/v3两包的深部对应明显较弱；不能用1m结果替它们签认，也不能把三包同时改变的几何、域宽和航高当成纯航高实验。

## 修复后的可执行处理入口

[reprocess_line9_result_packages.py](../../scripts/reprocess_line9_result_packages.py)读取全部406份native，以现有逐文件审计作身份约束：

- 实际源波形、源/接收器Yee时间偏移、0.025m源空间归一化；严格20–170MHz/0.3MHz/501点，保留复谱。
- Hann作为当前诊断基线，另存Gaussian/Blackman对照。不是已核验的设备频窗，也不是给ML新增频域算子。
- 无尾渐消、时间分割、背景扣除、AGC、叠道、包络相减、逐道归一化或额外零时调整。
- 保存复基带、复带通和`2Re[复带通]`有符号数据；灰度图仅对复幅度取对数，全部面板共享各包的同一Gaussian峰值参考。
- 只保存完成道；另存计划站位及完成掩码，不插值补道。图只画实际完成区，不能当作完整195m测线。
- pkg2清单x坐标与native裁剪域x相差80m；输入卡和native逐项核对，偏移单独记录，不用全局坐标索引局部H5。二维TMz中名义z=0.0125m被存为唯一平面z=0，不能误当几何位移。

三份数值文件位于本机忽略目录`D:\自动处理\artifacts\local_checks\2026-10-08_line9_corrected_products_r3\`：

|文件|已完成道/计划道|大小约|原始精度|
|---|---:|---:|---|
|line9_pkg1_v5_d02m_194st_corrected_sfcw.h5|194/966|95.5MB|FP32|
|line9_pkg2_x160-110_agl35_142st_corrected_sfcw.h5|142/251|70.1MB|FP32|
|line9_pkg3_x160-120_agl1_70st_corrected_sfcw.h5|70/201|34.6MB|FP32|

这些是派生数据，不覆盖旧包。complex128输出不提升原求解精度，场/电流矩代理也不是端口S21。HDF5完整逆变换周期为3333.33ns；原记录分别仅约1200/800/800ns，周期长度不等于已求解时间或探测深度。

常用键：`frequency`、`response`、`station_id`、`chainage_m`；Hann为`time_response/time`、`time_response/complex_envelope`、`time_response/complex_bandpass`、`time_response/real_bandpass`；另两窗为`diagnostic_gaussian`、`diagnostic_blackman`。第一维时间/频率，第二维完成站位。绘图不参与数值变换。

## 1m包为什么现在能对应上

前序标线是95MHz相位射线时间。对于Debye色散、频率相关反射与损耗，该时间不等于带限包络峰位。本轮用逐道实际层厚、完整复折射率、正入射Fresnel系数及上覆双程传播构造局部一次反射近似：

`H_primary(f) = (overlying down/up transmissions) × R_target(f) × exp[-i·4πf Σ(n_j(f)d_j)/c]`

空气段加入旧相位射线的非垂直路径小修正。该模型不含二维扩散、坡面/侧向散射、多次反射、PML和实体天线，绝不当成新的FDTD。回波位置不是通过实际数据调延时拟合出来的。

1m包底砂Hann近似峰位沿完成测段为286.09→255.32ns，比旧相位线提前约2.61–4.17ns；实际峰在这个近似的±12ns参考窗内，最大偏差约1.66ns。上一界面近似峰与底砂间隔约44.1–49.9ns，不是把紧邻的上一界面直接改名为底砂。

|底砂参考窗复数波形相关|pkg1 / 8m / v5|pkg2 / 固定y35m / v3|pkg3 / 1m / v3|
|---|---:|---:|---:|
|Gaussian|0.028|0.099|0.120|
|Hann|0.558|0.496|**0.980**|
|Blackman独立复算|0.577|0.506|**0.985**|

1m包首个砂岩顶的Hann相关约0.825。靠近测段末端它与覆盖层底趋近，不能因此声称这两条已完全分辨。

这里的相关保留复相位，定义为`|Σconj(model)·actual|/(||model||·||actual||)`，不是能量占比。所有结果在**由模型给出的±12ns窗**计算，属于模型辅助的对应诊断，不能叫盲检。

### 反例与独立检查

- 把1m底砂模板改为水平走时，相关降至0.00450；改为反向走时，降至0.01475。两者均无逐道拟合。
- 用前半站位估一个全局复比例，再检查后半站位，Hann/Blackman相关仍为0.980/0.985；幅相拟合相对L2误差约0.214/0.192。这是站位系数检验，不是地质母模型独立测试。
- 1m右侧PML单次空气镜像路径沿线增加约92.1ns，而底砂近似峰减少约30.8ns；该简单路径本身不能解释这条反向变化的回波。不排除其他边界路径或数值误差，也不签认PML反射幅度。
- 9个固定站位（每包首/中/末）直接求和复算全部501频点，最大相对L2约1.97e−12；三包每道另有7频点独立检查。弱底砂窗显式逆求和相对L2最大约1.54e−10，均通过预设1e−9检查。
- 无色散二层的解析相位、Fresnel透射/反射和峰位反例通过；HDF5写后回读、复数/有符号数据及完成掩码检查通过。
- 四项拒绝检查通过：输入目录用作输出、输出目录别名、native哈希变化、材料哈希变化。处理后406份native哈希全部不变。

**以上支持1m包存在与底砂近似一致的真实晚时记录趋势；缺少同域同站位底砂消融及FP64收敛，仍不能作独立物理归因证书或实测性能结论。** 本轮没有人为注入回波来生成此结果，也没有依据模型曲线修改任何实际数据。

## 较高航高结果为什么仍不签认

两包的Hann底砂相关仅约0.50–0.56，前半系数在后半的误差也很大。局部水平层近似失配、横向传播、多次反射和弱场数值可信度仍未独立拆分。材料损耗预算仍沿用前序，不能用换频窗消除介质本身的衰减。

还检验了覆盖层内部的第二次往返回波：一次覆盖层底反射再乘`R_cover→air·R_cover→below·exp(-i4πfn_cover d_cover/c)`。8m包近似峰378.41→292.75ns，确实接近部分深部斜带，但三包Hann/Blackman复数相关仅约0.41–0.49。**到时接近不足以认定它主导，更不能照旧做包络相减去多次波。** 本轮只记录候选，未据此删信号。

现阶段不继续扩道、不把“高飞必然不行”当结论。继续物理归因时，最小下一步仍是一个代表站位、同域同站位FP64底砂保留/移除配对；它不是本轮已执行事项。现有数据中的后处理问题已有可执行修正，尚不具备本地无求解手段让两份较高航高结果获得独立目标归因。

## 图件、数值和复现

下列均为同一实际数据的处理对照，不是新增正演或训练样本。第一组中间列无参考线，右列才叠近似层位：

- [pkg1无标线/全界面对应](../../artifacts/research_checks/2026-10-08_line9_layer_kinematics_r1/line9_pkg1_v5_d02m_194st_unmarked_and_contacts.png)
- [pkg2无标线/全界面对应](../../artifacts/research_checks/2026-10-08_line9_layer_kinematics_r1/line9_pkg2_x160-110_agl35_142st_unmarked_and_contacts.png)
- [pkg3无标线/全界面对应](../../artifacts/research_checks/2026-10-08_line9_layer_kinematics_r1/line9_pkg3_x160-120_agl1_70st_unmarked_and_contacts.png)
- [1m三个固定站位：实际波形与近似，两个频窗](../../artifacts/research_checks/2026-10-08_line9_corrected_independent_r2/pkg3_basal_actual_and_primary.png)
- [pkg1修正交付/三个频窗](../../artifacts/research_checks/2026-10-08_line9_corrected_products_r3/line9_pkg1_v5_d02m_194st_corrected_sfcw.png)
- [pkg2修正交付/三个频窗](../../artifacts/research_checks/2026-10-08_line9_corrected_products_r3/line9_pkg2_x160-110_agl35_142st_corrected_sfcw.png)
- [pkg3修正交付/三个频窗](../../artifacts/research_checks/2026-10-08_line9_corrected_products_r3/line9_pkg3_x160-120_agl1_70st_corrected_sfcw.png)

数值：[层位与反例](../../artifacts/research_checks/2026-10-08_line9_layer_kinematics_r1/layer_kinematics.json)、[修正交付及文件哈希](../../artifacts/research_checks/2026-10-08_line9_corrected_products_r3/corrected_processing_report.json)、[独立复算](../../artifacts/research_checks/2026-10-08_line9_corrected_independent_r2/independent_checks.json)、[拒绝/406道哈希复核](../../artifacts/research_checks/2026-10-08_line9_reprocessing_guards_r1/guard_checks.json)、[覆盖层多次路径候选](../../artifacts/research_checks/2026-10-08_line9_cover_reverberation_r1/cover_reverberation.json)。

原ZIP、地质H5和派生数值H5均留私有忽略目录；Git包含代码、聚合报告和图。原始native/SFCW逐文件清单及其科学数值不修改。

```powershell
D:\gprmax_v4_gpu_env\Scripts\python.exe scripts/reprocess_line9_result_packages.py --root artifacts/local_checks/2026-10-07_line9_three_packages_review_r1 --review artifacts/research_checks/2026-10-07_line9_three_packages_review_r3 --out artifacts/local_checks/line9_corrected_reproduce --report-out artifacts/local_checks/line9_corrected_reproduce_report
D:\gprmax_v4_gpu_env\Scripts\python.exe scripts/diagnose_line9_layer_kinematics.py --root artifacts/local_checks/2026-10-07_line9_three_packages_review_r1 --review artifacts/research_checks/2026-10-07_line9_three_packages_review_r3 --cache artifacts/local_checks/2026-10-08_line9_postprocessing_cache_r1 --out artifacts/local_checks/line9_kinematics_reproduce
D:\gprmax_v4_gpu_env\Scripts\python.exe scripts/check_line9_corrected_products.py --root artifacts/local_checks/2026-10-07_line9_three_packages_review_r1 --review artifacts/research_checks/2026-10-07_line9_three_packages_review_r3 --products-report artifacts/research_checks/2026-10-08_line9_corrected_products_r3/corrected_processing_report.json --kinematics artifacts/research_checks/2026-10-08_line9_layer_kinematics_r1/layer_kinematics.json --out artifacts/local_checks/line9_independent_reproduce
D:\gprmax_v4_gpu_env\Scripts\python.exe scripts/diagnose_line9_cover_reverberation.py --root artifacts/local_checks/2026-10-07_line9_three_packages_review_r1 --review artifacts/research_checks/2026-10-07_line9_three_packages_review_r3 --cache artifacts/local_checks/2026-10-08_line9_postprocessing_cache_r1 --out artifacts/local_checks/line9_reverberation_reproduce
D:\gprmax_v4_gpu_env\Scripts\python.exe scripts/check_line9_reprocessing_guards.py --root artifacts/local_checks/2026-10-07_line9_three_packages_review_r1 --review artifacts/research_checks/2026-10-07_line9_three_packages_review_r3 --out artifacts/local_checks/line9_guards_reproduce
```

各输出路径必须不存在。前两条分别独立生成修正数值和近似诊断；第三条示例核对本次实际交付，核对另一次重建时改为其报告路径。跨机需另带私有原包；仅Git不足以重建。依赖已核验V4 SFCW工具箱/NumPy/h5py/Matplotlib，未另装环境。

开发过程：修正处理r1因把二维z唯一平面当成名义半网格位置而拒绝；r2在pkg2全局/局部x清单差异处拒绝。两次均未修改原包，r2仅有一个包的中间产物，没有当作整批成功。修正后的r3完成三包。独立图r1的统一纵轴裁掉最后一站峰顶，r2改为由六面板共同峰值确定纵轴，数值未变；这些中间产物归档于私有`artifacts/local_checks/2026-10-08_line9_processing_attempts/`。
