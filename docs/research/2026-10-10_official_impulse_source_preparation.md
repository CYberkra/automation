# 按用户授权将新输入改为官方内置 impulse

2026-10-10。用户在明确询问激励尚未修改后指示“改为官方的”。本轮完成源替换输入、对应源验收代码及 CPU 检查。**输入已改好，尚未运行新激励正演，既有 B-scan 仍来自 Ricker。**

## 官方依据与改动

已读取 [gprMax 4.0.1 官方 SFCW 文档 Preparing the model / Information](https://docs.gprmax.com/en/latest/inc_SFCW.html#preparing-the-model)、本机官方 `cylinder_sfcw_2D.in` 和 `waveforms.py` 的 impulse 分支。采用官方示例的内置单样本单位冲激：

```text
#waveform: impulse 1 1 impulse
#hertzian_dipole: z 169.35 38.95 0.0125 impulse
```

第一个1是单位幅度，作为 Hertzian 电流源为1 A；第二个1是必填频率语法占位，impulse不使用中心频率，**不是1 Hz发射波**。不自写 `t==0` 函数；内置实现能在电源首个 `dt/2` 样本正确激发。源起始时间保持0，源坐标、方向和0.025m单元长度保持原值。单位源不代表实际设备发射功率已经标定。

SFCW继续使用官方推荐的 `direct`：准确20–170MHz、0.3MHz、501点，保留实际源/接收时间原点与复相位。Hann/Blackman、补零8、尾窗0、显示移位0以及场/电流矩参考保持既有约定。用户本次指令针对激励；本轮未额外切换到homodyne或改变其他处理参数。

现行源约定为[配置 v0.2](../../configs/research/line9_source_v0_2.json)。[准备脚本](../../scripts/prepare_line9_official_impulse.py)读取已有模型包，检查父输入/几何/材料哈希，另存新包，只改波形定义和偶极的波形ID，拒绝多波形、错绑定或延迟源。该脚本没有求解入口，不生成 execution contract 或复用旧求解输出。

## 已准备的输入

复用已核验的190m高损耗站位 `high_x19000_H0` / `high_x19000_H1`，只用于展示并验收这次源替换；没有启动两项计算。域210×42.5m、dl0.025m、时窗1200ns、FP64约定、HORIPML80格、材料/几何、约8m离地高度、收发位置保持；不新增探针或快照，不改变横向基线未表示的二维限制。

本机完整包：`artifacts/local_checks/2026-10-10_line9_official_impulse_package_r1/`。源准备 manifest、两份输入文本和验证结果入库：`artifacts/research_checks/2026-10-10_line9_official_impulse_r1/`。该公共目录是输入及证据摘录，几何/材料仍在完整模型包内；公共输入文本不能脱离其几何/材料单独运行。

可在有原模型包的电脑上复现准备（新 out 目录必须不存在）：

```powershell
& artifacts/local_checks/gprmax_v401_gpu_env/Scripts/python.exe scripts/prepare_line9_official_impulse.py --parent artifacts/local_checks/2026-10-08_v401_dense_loss_prepared_r1 --out artifacts/local_checks/line9_impulse_new --groups high_x19000_H0 high_x19000_H1
```

它也可对同一类已有包显式选择其他站位，逐项保留原输入的几何与数值配置；本轮没有把旧批次整体转换或恢复。

## 验收与限制

官方输入解析器在两份新卡中均读到 `wave_type=impulse, amp=1, freq=1, id=impulse`；偶极波形ID一致，start/stop未加延迟。独立逐行对照只变两行；几何和材料字节哈希不变。

五项 CPU 回归：[check_line9_official_impulse.py](../../scripts/check_line9_official_impulse.py)。检查源替换的LF/CRLF保持及其他输入不变，拒绝歧义/延迟/错绑定/方向输入，实际调用官方内置波形证明首个半步为1、其余255样本为0，验证501点源谱为 `dt·exp(−i2πf·dt/2)`，拒绝假多样本冲激/时间偏移/波形元数据，以及保留历史Ricker验收。

`line9_v401_version_controls.audit_source`已按批次声明核验波形类型、幅度和频率字段，新增impulse首样本/起始时间/半步检查。旧Ricker批次仍按原参数检验；本机已有真实Ricker H5通过。未改任何旧执行契约、原生输出或历史哈希；旧冻结批次按其原源码提交复现，新批次必须冻结更新后的验收代码身份。

准备检查不等于新正演通过，更不代表条纹/实测差距改善。换源后必须重新正演才能取得该源的接收历史；不能把旧Ricker H5元数据改成impulse。单冲激还激发目标频段外的数值响应，官方明确其宽谱不自动认证所有频率；正式结论仍限定于经网格、时窗和源/接收原点检查的20–170MHz。

后续新的求解批次、homodyne对照或其他模型改变，仍按用户要求先列明具体范围与成本并取得许可。当前 `execution_authorized=false` 只表明本轮是源输入准备，不沿用历史自主求解授权。
