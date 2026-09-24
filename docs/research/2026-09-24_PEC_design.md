# PEC平面独立参考对照 V4-PEC-01

## Material Passport
研究验证设计；2026-09-24；AI辅助执行；用户授权自主设计/GPU仿真；非现场材料参数。一个确定性对照，不涉及统计样本或显著性推断。

问题：当前24m立方域、5cm网格及有限记录，能否在20–170MHz再现已知平面反射？

只将M01_T800的εr9半空间换为内置pec，删除自定义材料。位置、800ns记录、20格HORIPML、X电流元/Ex、1.3m横基线及15m高度保持。输入hash与环境hash见gate。CUDA0/double，30分钟/32GiB host committed/1GiB输出，启动空闲RAM32GiB、VRAM12GiB，单次0重试，无CPU回退。原空气T800复用，不重跑。

独立参考：地面z=4m，源z=19m的镜像z=-11m，水平电流反向；镜像到Rx距离sqrt(30²+1.3²)m，Rx仍在X电流元赤道面。精确无限PEC反射Ex/I = -transverse_dipole(f,r,0.05)，含近场项、不用远场近似、不拟合高度/时移/相位/幅度。实际有限网格PML域与无限平面解差异正是被检查对象。有限电流元与点元差异、边界和时间窗仍可能影响结果。

方法来源：[Michigan EECS530 Image Theory](https://www.eecs.umich.edu/courses/eecs530/lec%206.pdf)；[Cornell完整电流元场](https://courses.cit.cornell.edu/ece303/Lectures/lecture28.pdf)。官方SFCW沿用实际源历史/时间偏移和direct方法。

预声明分析：保持官方800ns截取及0/200/400ns尾渐消，PEC与空气配对相减，比较501点复数响应与解析反射；报告最大/中位复数相对误差、幅度dB、相位deg，无事后物理合格阈值。保留原始提前响应及频谱诊断。不能用总场误差掩盖弱反射误差；不能用PEC结论替代介质半空间结论。不扩大为网格/PML扫描。执行、分析、日志、hash、接续说明均归档。
