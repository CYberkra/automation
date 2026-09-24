# gprMax V4 GPU环境与复建

## Material Passport

- 用户要求：安装配置环境，后续使用GPU运行。
- 状态：PyCUDA安装、真实CUDA双精度核函数、V4导入与pip依赖检查已通过；首个CUDA M00校准已正常完成。
- 位置：`artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe`。
- 后端策略：FDTD必须CUDA，缺少GPU/编译器/显存则失败，不自动回退CPU。
- 审计：`artifacts/research_checks/2026-09-24_gpu_setup/`，包含失败与成功日志、版本、身份、哈希。

## 固定配置

Windows、Python3.12.14、gprMax4.0.0（与前次相同的本地wheel）、PyCUDA2026.1、CUDA Toolkit13.3.73、MSVC14.44.35207、RTX4090 Laptop 16GiB、驱动610.88。Python依赖精确版本见`configs/research/gprmax_v4_gpu_win_py312_dependencies.txt`。原V4源码和CPU环境未修改，GPU环境另建。

[V4官方加速器文档](https://docs.gprmax.com/en/latest/accelerators.html)区分CPU精度与GPU精度，CUDA默认single。本项目本次显式指定`-gpu 0 -gpu_precision double`，保持物理输入和dtype不变。GPU负责FDTD时间步运算；建模、编译、文件读写仍需要主机CPU。官方SFCW目前以NumPy/SciPy实现，不能把其后处理称为CUDA计算。

## 安装证据

PyCUDA由[PyPI官方2026.1源码发行包](https://pypi.org/project/pycuda/2026.1/)在本机编译；源码SHA256为`759516160628ba06f32ce7e563e3f5b9214691dc9528a03ea99ea1073f4e14ba`。本机wheel SHA256为`49f9b291e88d76b4a1a20a602833f8e34b7f8e44d011223109a6e835af635189`。重建wheel不保证二进制逐字节相同，应重新归档身份。

首次尝试加载官方`activate_cuda.ps1`被本机脚本签名策略拒绝，接续核函数编译因找不到cl.exe失败。未修改安全策略；改用官方文档同样支持的MSVC x64 Native Tools命令环境（`vcvars64.bat`），CUDA核函数编译和执行通过。C4819字符编码警告保留在日志中，不隐瞒；无数值校验错误。

真实GPU核函数处理4096个float64整数样本，平方加1，与可精确表示的参考逐项一致。首次CUDA事件计时包含启动/JIT等影响，不作为FDTD性能基准。

## 复建与启动

1. 在新机器建立Python3.12独立venv，安装上述锁定依赖；Cython设置`NO_CYTHON_COMPILE=true`并使用`--no-binary=Cython`。PyCUDA需要CUDA工具包和兼容MSVC。
2. 按前次[V4构建记录](2026-09-24_v4_build_success.md)取得/构建经过核对的V4 wheel，再安装到该环境。不要默认下载同名新版替代本次源码。
3. 在MSVC x64环境中运行`inspect_cuda_runtime.py`及`inspect_v4_runtime.py`，输出到新路径；核对GPU和版本，并更新契约的本机路径/哈希。
4. `scripts/run_airborne_gpu.cmd`初始化编译器并调用当前批准的启动器；默认仅前检，带`--execute`才求解。已使用的attempt禁止重复运行；另一次必须创建并留档新运行契约。

本次CUDA M00预算为30分钟、32GiB主机Job提交内存、1GiB输出、8线程建模准备、至少12GiB空闲显存。显存阈值是启动前检查，**不是硬性显存配额**。观察时GPU总使用约10GiB、利用率98%；这些是整卡瞬时指标，完整遥测随结果归档。

## 首个GPU结果

4156步全部完成，退出码0，监督器总耗时262.719秒（约4分23秒）。主机Job提交峰值20.883GiB；运行中38次整卡遥测的显存使用最高10259MiB、利用率最高99%，采样从启动后开始，不能视为精确进程峰值。原始HDF5已通过10项元数据/数值基本检查：V4、坐标、float64、有限非零、单样点源、长度与时刻、偶极长度均符合契约。

完整证据位于`artifacts/research_checks/2026-09-24_M00_x_3d_gpu_r1/`。启动入口复查已正确拒绝重复attempt。官方SFCW仍未执行，等待用户对GPU要求是否包含后处理的澄清；原始记录末5%相对峰值为−41.414dB，后续必须分析有限时窗与静电尾项，不能直接当作合格频域校准。

## 尚不能声称的事项

两次CPU试跑均主动中止，没有CPU完整输出，不能据此声称已通过完整CPU/GPU逐样本一致性。当前CUDA算例仍要进行输出元数据和独立物理核查，也没有证明20m目标可探测、实际天线匹配或网格/PML收敛。后续如改变单精度、材料或域大小，要记录独立误差与资源依据。
