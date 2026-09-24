# 固定PML因果对照 V4-FIXED-PML-01

## Material Passport
2026-09-24；AI辅助数值机制验证；用户对固定PML对照明确“做吧”“继续”；一个新GPU算例，确定性实验，不读取实测或作统计推断。

问题：77ns原始提前差分是否主要来自空气/介质侧PML自动参数不同？

保持M01_T800所有几何、介质εr9/σ0、源/接收、5cm网格、24m立方域、800ns记录、20格HORIPML，仅增加官方命令：
```
#pml_cfs: constant forward 0 0 constant forward 1 1 quartic forward 0 0.21235349838321013
```
数值为本地官方CFS.calculate_sigmamax(.05,1,1)所得空气值。默认alpha=0,kappa=1,sigma四次正向；明确数值避免按介质截面平均自动调整。仅一阶，不叠加默认项。官方命令实现与文档均核对：[PML commands](https://docs.gprmax.com/en/latest/input_hash_cmds.html#pml-cfs)。

用已安装官方CFS函数在小数组上验证固定字面值与空气自动值的alpha/kappa/sigma共6个E/H参数数组逐位相等。该检查未构建网格、未调用FDTD，记录在fixedPML_sources/equivalence.json。因此复用T800原空气输出，不重复空气求解。当前介质自动侧sigma为空气0.654654倍，底面1/3，顶面相同；此全局固定命令同时改变侧面与底面，不能声称只改侧面。对0–80ns提前响应的时间隔离是本轮重点；晚期结果须考虑底面也变化。

预声明：新M01_fixedPML减相同空气T800，与旧M01_T800减空气比较。报告0–80ns最大绝对差、RMS、抑制比；80–100/100–140/140–400/400–800ns保留窗口诊断。0–80ns是前轮已用窗口，不根据新结果选取。官方direct处理501点20–170MHz，对0/200/400ns尾渐消均保存复数差分；比较新旧频谱的每频点相对变化及以旧频谱峰值归一的最大变化，不把变化当解析准确度。若早到显著减少，支持PML自动参数机制；如果不减少，不强行归因。无物理合格阈值拟合，无裁剪修饰原始场。

CUDA0/double，30分钟、32GiB主机Job committed memory、1GiB输出、启动可用RAM32GiB/VRAM12GiB，1次0重试，不回退CPU。环境、输入、launcher、supervisor hash由gate锁定。保留原始HDF5/日志/执行记录及分析脚本。后处理官方CPU。

边界：固定参数为机制隔离，不是介质PML最优设计或网格收敛；已知170MHz长空气路径约5.4°数值相位偏差仍在。不能据此认证真实粘土砂岩或20m探测能力。
