# 独立V4环境尝试与Windows执行限制验证

2026-09-24。用户同意继续环境/执行程序准备；无仿真授权。**未运行FDTD、geometry-only或任何gprMax测试，未导入求解器。**

## 环境实际完成到哪里

已用现有Python3.12.14建立独立venv：`artifacts/local_checks/gprmax_v4_env`。把用户所供V4源码复制到忽略目录`artifacts/local_checks/gprmax_v4_build/source`，记录2236份文件的SHA-256；在副本上构建wheel，原始源码目录和V3环境未修改。构建并行度设为2。

构建失败于pip隔离构建环境的Cython导入，尚未生成V4 wheel或安装V4。Windows `Microsoft-Windows-CodeIntegrity/Operational` 在本地20:01:21记录事件3077/3033：`Cython/Compiler/Code.cp312-win_amd64.pyd`不满足签名/代码完整性策略。错误不是FDTD输入错误，也没有证据归因于V4数值代码。

没有禁用、修改或绕过策略，没有把被阻止文件搬到别处尝试规避。后续需要管理员批准的原生扩展环境/签名策略处置，或另一个明确允许该依赖的计算环境；不把未经批准的“继续”视为系统安全策略修改授权。这个阻断解决后仍须完成V4构建、模块路径/版本/依赖身份与导入检查。

[记录、原始构建日志和源码清单](../../artifacts/research_checks/2026-09-24_v4_environment_preparation/record.json)已入Git。pip隔离依赖没有成功完整固化，不能把build-system的宽松依赖要求称为可复现的已验证锁文件。原始日志中本地化错误字符有编码损失，Code Integrity事件用于独立确认原因。

## 执行限制程序与实际检查

[bounded_windows_process.py](../../scripts/bounded_windows_process.py)使用Windows Job Object。子进程先以隐藏、挂起方式启动，加入Job后再恢复，避免启动到归入监管之间自行创建未受管子进程；退出时关闭Job清理后代。内存由JOB_MEMORY实施为**进程组提交内存上限，不是RSS**；超时和文件大小轮询检查，文件上限不是磁盘硬配额，可能超调。

[check_windows_supervisor.py](../../scripts/check_windows_supervisor.py)只运行普通Python样例，6项通过：正常完成、非零退出、超时停止、后代进程确已退出、超过96MiB提交预算的256MiB申请被拒绝、输出超限停止。没有跑gprMax。测试结果保留命令、退出码、墙钟时间、观察到的提交内存峰值和文件大小；见[结果](../../artifacts/research_checks/2026-09-24_v4_environment_preparation/supervisor_results.json)。

内存超限可能让Python抛MemoryError或让原生进程失败，并不保证可标成特定退出码；监督器保留非零退出与日志，不伪造精确失败原因。程序本身是通用执行函数，没有公开的自动求解入口；未来求解包装器仍须绑定全局批准记录、包/输入哈希、实际V4环境身份、可用内存前检及CPU线程设置。**目前没有把批准检查和求解命令接起来，因此不能宣布整个仿真执行链已经就绪。**

首跑原先“24GiB RSS”拟值须在最终提交时改为“24GiB Job提交内存”；该语义修改不自动获得用户批准。普通测试使用小资源限制，不是24GiB/30分钟真实压力测试。

复核只需标准库：

```text
python scripts/check_windows_supervisor.py --output artifacts/local_checks/<新的检查目录>
```

输出目录必须不存在。源码或限制逻辑有变化才重跑；无需重复无变化的6项检查。已知跟踪PNG缺失仍留作工作区状态，未擅自恢复；本轮新增文件和哈希单独检查，不声称整体handoff已通过。

## 当前待办顺序

1. 解决Windows原生扩展签名/策略阻断，或选择用户批准的可用环境。
2. 完成V4独立构建与只导入身份验证、环境锁定。
3. 完成批准记录与Job监督器的连接、实际源/接收元数据验收代码，再把单次M00_x_3d、参数哈希和资源上限交用户敲定。
4. 只有明确批准后才启动该次FDTD，不因环境修好自动运行后续模型。
