# 独立全波半空间参考与离散误差预算

## Material Passport
2026-09-24；AI辅助数值研究；确定性频域校验。未运行新的FDTD，未读取实测测线，未修改gprMax环境或系统安全策略。独立参考数值检查及复算通过，尚无FDTD网格收敛认证。

## 主要结论

对已有M01_fixedPML的501点20–170MHz反射频谱，独立无限介质半空间全波积分给出：最大复数相对差9.6172%、中位1.7978%，最大幅度差0.22308dB，最大相位差5.38594°，主要高频端偏差在170MHz。固定PML已排除前轮提前差分，此处比较的是保留完整近场项的反射频谱，不是直达场主导的总场。

| 频率 | 复数相对差 | 幅度差 FDTD/参考 | 相位差 |
|---|---:|---:|---:|
| 20MHz | 0.04040% | −0.003273dB | −0.008358° |
| 95MHz | 1.79780% | −0.066841dB | −0.935328° |
| 170MHz | 9.61717% | −0.223082dB | −5.385942° |

原始FDTD频谱没有经验校正、拟合相位或时间对齐。物理准确度阈值仍未设定，不能把“数值检查通过”写成“达到工程精度”或20m可探测。

## 参考的来源和方法

使用官方empymod1.10.6的未修改全波wavenumber内核，外加本项目的分支点分段Gauss积分包装。它是独立于gprMax FDTD的半解析频域积分参考，**不是闭式半空间解析解，也不是empymod默认Hankel变换，更不是对官方SFCW流程的替换**。

几何以地表为0、向下为正：源(0,0,−15)、接收(0,1.3,−15)m，两者X向；源电流元长度0.05m。上层εr1、下层εr9、μr1、各向同性，导电率严格为零。直接构造频域eta/zeta并仅计算反射部分；没有用有限大电阻率近似，也未采用空气扩散近似。

官方内核返回PJ0/PJ1/PJ0b；横向角pi/2的被积函数为(PJ0−PJ0b)J0−PJ1 J1/r。传播段用lambda=k0 sin(theta)，倏逝段用lambda=sqrt(k0²+u²)，在下介质分支点进一步分段。Gauss节点避开端点；对称和Bessel系数沿用官方内核约定。[官方1.10.6内核](https://github.com/emsig/empymod/blob/v1.10.6/empymod/kernel.py)明确允许直接调用，但由调用者保证输入正确。代码与[预先设计](2026-09-24_halfspace_reference_design.md)均留档。

当前官方[dipole文档](https://empymod.emsig.xyz/en/stable/api/empymod.model.dipole.html)用于初始接口调查；实际1.10.6接口以安装包源码核对，旧版本使用fht而不是dlf，错误尝试没有冒充有效结果。

## 参考自身的核验

- Gauss128/256相对512阶的全频带最大变化分别1.8593e−8、1.5265e−8；倏逝积分上限U=2改3变化1.0776e−8。这是数值敏感性，不是严格误差上界或完全独立的多软件验证。
- 官方全空间函数对完整Hertzian公式最大相对差5.4440e−10，支持单位、源长度、相位约定和横向几何正确。
- 两层参数相同时反射逐点为零。
- 下层导电率从1e6增至1e10S/m（只用于PEC极限测试），与严格PEC镜像参考的5频点最大相对差从1.9432e−4降为1.9447e−6，符合趋近PEC。介质主参考仍σ=0。
- 参考有效性检查阈值在脚本中执行前固定：积分变化<1e−5、全空间差<1e−7、PEC极限差<1e−5且随导电率增大减小；这些阈值不用于给FDTD判合格。
- 501频点积分和核验部分实测约0.59秒（不含依赖配置、导入、绘图）；CPU足够，无GPU迁移理由。
- 独立重复分析后，除耗时外JSON一致，CSV和图hash一致，所有参考数组逐位一致。

默认积分的失败也保留：5频点fht结果严重不一致，qwe给出全零，quad发出未收敛警告且结果不一致。全部拒绝用作真值，而不是择优隐藏。完整日志位于halfspace_reference目录。

## 幅相误差的机制解释

此前读取PEC结果前保存的轴向Yee相位预算，在170MHz给出−5.41755°；本轮实际−5.38594°，全频带最大残差0.03160°。这支持长路径网格色散解释，但不是细化收敛实验。

读取本轮结果后，追加无拟合正入射界面离散反射预算：界面电场节点εr=(1+9)/2=5，qj=2asin(sqrt(εrj)sin(ωdt/2)/(c dt/dz))，Γnum=(sin(q1)−sin(q2))/(sin(q1)+sin(q2))。该式等价于[Schneider第7章](https://eecs.wsu.edu/~schneidj/ufdtd/chap7.pdf)7.96的平均节点反射式。

170MHz预测幅度相对连续正入射Γ=−0.5低0.21969dB，实际低0.22308dB；全频带幅度预算残差最大0.003391dB。该检查是事后提出的机制解释，无拟合参数，未用于校正数据；正入射平面波不等于完整三维点源参考，仍需细化验证。

## 下一项明确工作

先评估并执行一个方向性网格细化的空气/介质配对，候选dx=dy=5cm、dz=4cm，保持24m域与1m物理PML厚度（侧面20格、上下25格），双方显式相同PML。单元从110592000增至138240000（+25%），静态轴向相位预算由5.418°降至2.916°，正入射界面幅度预算由−0.21969降至−0.13747dB。**这些是预测，不是已运行的改善结果。**

旧空气结果不能直接复用：细化后CFL dt改变，需要新的空气/介质配对。显存可按既有实测作初估，但必须核对具体PML/临时数组和主机Job上限后再形成执行契约。dz=2.5cm单元翻倍，当前16GiB GPU很可能不足，不直接启动。此轮没有新增执行契约，现gate保持上一项已完成/attempt已消耗。

## 环境、失败记录和复现

独立环境`artifacts/local_checks/halfspace_reference_env`，Python3.12；使用empymod1.10.6、NumPy1.26.4、SciPy1.13.1、Matplotlib3.9.4。[实际依赖锁定](../../configs/research/halfspace_reference_dependencies.txt)列出所需包；environment_freeze.txt另保留所有安装包，含未使用的2.6尝试残留依赖。pip check通过。

最初empymod2.6.0导入Numba原生组件时被Windows应用控制拒绝，未获得参考输出。切换1.10.6后首次NumPy导入也曾报同类错误；随后同一未改动环境的积分、全空间与复算导入成功，中间没有安全策略更改，首次失败原因未进一步归因。失败日志保留。这证明本机已完成核验，不保证新机器自动通过原生模块策略。

原始gprMax/SFCW环境未被修改；FDTD仍必须GPU。这里只有独立CPU频域积分校验。

```powershell
python -m venv artifacts/local_checks/halfspace_reference_env
& artifacts/local_checks/halfspace_reference_env/Scripts/python.exe -m pip install -r configs/research/halfspace_reference_dependencies.txt
& artifacts/local_checks/halfspace_reference_env/Scripts/python.exe scripts/compare_halfspace_reference.py --input artifacts/research_checks/2026-09-24_fixedPML_analysis/frequency_comparison.csv --output artifacts/local_checks/halfspace_replay_new
```

已存在且验证好的环境不要重建；输出目录须不存在。`scripts/budget_interface_discretization.py`重算后验机制/细化预算。完整结果在`artifacts/research_checks/2026-09-24_halfspace_comparison/`；失败/默认算法诊断在`2026-09-24_halfspace_reference/`；record.json含hash。没有训练集或实测性能结论。
