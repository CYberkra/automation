# V4独立构建与导入验证完成

2026-09-24。用户明确“许可你构建”。此次授权用于构建/安装和只导入验证，**不改变仿真需另行敲定的约定**。

## 结果与可信范围

已从用户提供的V4源码副本构建并安装 `gprMax 4.0.0`，独立环境位于 `artifacts/local_checks/gprmax_v4_env`，Python3.12.14。原V3环境未改变；原始V4源码2236个文件重新核对SHA-256均未改变。

验证通过：安装成功；`pip check`无损坏依赖；实际导入路径位于独立venv；版本4.0.0；四个本地编译模块（普通场更新、几何primitive、电/磁HORIPML更新）均可导入；官方SFCW命令行解析接受20–170MHz、501点、direct。没有读取模型输入、分配FDTD模型或执行时间步，不能据此宣称数值/性能验证通过。

[导入身份](../../artifacts/research_checks/2026-09-24_v4_build_r2/runtime_identity_r2.json)、[本轮记录](../../artifacts/research_checks/2026-09-24_v4_build_r2/record.json)、[只导入检查脚本](../../scripts/inspect_v4_runtime.py)。安装包为 `gprmax-4.0.0-cp312-cp312-win_amd64.whl`，大小12,632,417字节，SHA-256 `ca10b53674e211cd6f85efdc7a428023ae8ce435b34817cbb639630fa2394ebd`。wheel留本地忽略目录，Git保留源码指纹、日志、wheel哈希和依赖版本，不上传二进制依赖。

## 怎样解决前次Cython阻断

前次失败是pip隔离构建环境加载Cython预编译扩展时被Windows代码完整性策略拒绝。本次采用[Cython官方支持的非编译安装方式](https://docs.cython.org/en/latest/src/quickstart/install.html)：`NO_CYTHON_COMPILE=true`，并对Cython使用`--no-binary=Cython`。实际生成 `cython-3.3.0-py3-none-any.whl`，导入的是 `Cython/Compiler/Code.py`。

随后在已配齐构建依赖的独立环境中，用`--no-build-isolation`构建V4，防止构建工具重新创建另一套依赖；这不是关闭Windows隔离或安全策略。V4自身的数值扩展仍经MSVC编译为原生`.pyd`，没有改写为Python慢速求解器。系统安全配置未修改，未对被阻止的原生文件做重命名/搬迁尝试。

构建依赖为Cython3.3.0、NumPy2.5.3、setuptools84.0.0、wheel0.48.0、Jinja2 3.1.6；MSVC2022工具目录14.44.35207，构建并行度2。全部26项依赖版本见[Windows/Python3.12清单](../../configs/research/gprmax_v4_win_py312_dependencies.txt)。这是一份本次观察到的版本固定清单，不是跨平台通用锁文件，也不保证不同机器编出的wheel逐字节一致。

## 首次导入未返回的记录

安装后的首次只导入检查数分钟没有输出，已确认命令行属于本任务后停止该检查进程。在所查时间段内未找到针对此venv的新Code Integrity拒绝事件；这不足以确定停顿原因，不能称已定位。

随后用同一脚本加 `faulthandler.dump_traceback_later(45, exit=True)` 再试，在45秒期限内成功返回，记录四个原生导入及官方SFCW解析结果。没有新增修改求解器的修复，也不将一次成功说成长期稳定性保证。以后启动仍需限时与日志。

## 新地点复建次序

1. 准备Python3.12及Windows C++工具链；取得用户认可的V4源码，核对已有源码清单。创建独立venv，复制源码到独立构建目录，避免在唯一原始目录生成文件。
2. 在该venv安装[依赖清单](../../configs/research/gprmax_v4_win_py312_dependencies.txt)，PowerShell当前进程设置 `$env:NO_CYTHON_COMPILE='true'`，使用 `python -m pip install --no-binary=Cython -r <依赖清单路径>`；检查 `Cython.Compiler.Code.__file__` 指向`.py`。
3. 设置 `$env:GPRMAX_BUILD_JOBS='2'`，执行 `python -m pip wheel <源码副本> --no-build-isolation --no-deps --wheel-dir <本地wheel目录>`。安装生成的wheel，然后`python -m pip check`。
4. 用新环境Python运行 `scripts/inspect_v4_runtime.py --output <新的JSON路径>`，建议套45秒外部/诊断超时。保留本机版本、模块路径/哈希、日志和实际依赖版本，不能直接复制本机“通过”的状态。

以上只含构建与导入步骤，不包含求解示例、pytest或FDTD。另一台机器的原生库仍受其本地策略约束；若被拒绝，保留失败原因并使用管理员认可的环境，不能保证此方法普遍解除所有阻断。

## 下一步

已完成构建任务。剩余是把已验证的Job执行限制程序与批准记录、输入哈希、此环境身份及元数据验收连接起来，再给用户提交单次 `M00_x_3d` 的启动范围。当前批准仍为false，其他输入不自动排队。既有设备PNG工作区删除状态仍保留，其原图在Git历史；不是本轮构建造成。
