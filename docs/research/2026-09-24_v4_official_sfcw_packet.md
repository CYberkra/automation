# V4-CAL-02-SFCW：采用官方集成流程

2026-09-24。当前有效草案，替代Ricker优先的V4-CAL-01。用户明确“SFCW在gprMaxV4有官方支持，优先官方集成好的”。**本轮仅准备，未运行仿真或官方后处理。**

## 主线已经确定

使用官方 `gprMax.toolboxes.SFCW`，不另写频域除法、混频、I/Q或逆变换实现。官方推荐的内置impulse激励提供离散系统响应，`direct`计算请求频点；`homodyne`仅作一条M00_base记录的官方独立交叉检查。数值impulse不等于设备发射信号，输出也暂不称硬件S21。

本地官方 README、CLI、processing、示例输入及版本文件共5份SHA-256保存在[新包清单](../../configs/research/v4_official_sfcw_packet_v2/manifest.json)。对应[官方SFCW文档](https://docs.gprmax.com/en/latest/inc_SFCW.html)。已静态核查CLI用`linspace(start,stop,steps)`，因此此处`--steps 501`是501个频点，而非501个间隔。

8个新输入使用 `#waveform: impulse 1 1 impulse`、零启动延时，官方说明中频率参数1不参与impulse波形。接收器命名measurement，输出Ez。旧Ricker输入保留为搁置历史，不运行，也不自动加入额外对照。

## 明确参数，不隐藏采用默认值

| 项目 | 本包提案 |
|---|---|
| 频率 | 20e6至170e6 Hz，501点，间隔300000Hz |
| 方法 | 官方direct |
| 接收 | `name:measurement`，Ez |
| 源 | 官方inspect确认唯一源后填写实际路径，不写死src1 |
| 频率窗 | rectangular，用于初始诊断；这是显式选择，官方默认是gaussian |
| 零填充/展示移位 | 1 / 0秒，先保存原始复响应，暂不美化时域显示 |
| 尾部taper | 0；不截除物理尾波来消除告警 |
| 源谱floor | 官方默认−100dB，显式登记，仅为数值除法保护，不认证精度 |
| 独立检查 | M00_base另用官方homodyne，cycles=8；比较同一频率/源/接收器定义 |

待批准、生成真实输出后使用的命令模板：

```text
python -m gprMax.toolboxes.SFCW inspect <actual_raw.h5>
python -m gprMax.toolboxes.SFCW process <actual_raw.h5> --source <inspected_source_path> --receiver name:measurement --component Ez --f-start 20e6 --f-stop 170e6 --steps 501 --method direct --window rectangular --zero-pad 1 --time-shift 0 --tail-taper 0 --source-floor-db -100 --output <new_processed.h5>
```

尖括号是审查占位符，不是可直接执行命令。独立homodyne使用相同记录另存文件，方法改为homodyne并显式`--homodyne-cycles 8`。不运行官方示例来替代本项目结果。

CLI在存在无效源频点时抛出错误，可能尚未写出处理结果；应保存失败日志/原始文件，不能假设一定已有mask文件，也不自动放宽floor。成功时保存官方复响应、I/Q、源和接收频谱、源有效性及各类重建量。设备的真实导出定义未核对前，不擅自挑其中一种当作CSV等价量。

## 小域与资源提案

两个场景：M00空气、M01理想无损εr=9半空间。一个收发位置，高0.5m、距0.4m。每场景base/fine/boundary/long四变体，共最多8次FDTD；具体几何、坐标和哈希全在清单。

- base：16×10m，网格0.025m，240ns；地表y=7m。
- fine：网格0.0125m，其他物理条件不变；PML保持0.3m厚度。
- boundary：20×14m，源/接收器/地表平移(+2,+2)m，其他条件不变。
- long：基准几何，480ns。时窗不由1/步长代替。

CPU double、8线程、串行；最多每次30分钟、累计求解120分钟，官方后处理累计另限30分钟；进程树RSS 8GiB、输出1GiB、零自动重试。超限或元数据异常停止并留档。以上为**拟上限，未实施也不是速度承诺**。CFL工作量估算沿用相同几何的旧草案，但旧Ricker频谱保护理由不沿用。

独立V4环境、Windows构建身份和执行监督器尚未准备完；需要先补齐再请用户敲定实际启动。生成器只写输入和清单，不是求解/后处理脚本，不会因运行生成器而触发FDTD。全局批准记录保持false。

## 评价必须随impulse修订

**不比较细网格与粗网格的单位样本impulse原始幅度来判断物理误差。** dt变化改变其脉冲面积；以官方实际源归一化的复频响为主，保留正确Yee时间偏移。两域、两网格和两记录长度分别比较复响应相对L2差异及逐频绝对误差；近零响应不以逐频相对误差代替。

暂提频带整体相对L2差异≤1%作为继续准备下一阶段的工程筛查值，仍待审查；不是现场精度、弱目标误差界或物理标签阈值。报告逐频误差可揭示整体指标掩盖的差异。M01−M00在相同变体、相同归一化的频域量上计算，并独立检查差分稳定性；强直达波稳定不能替代弱差分稳定。

170MHz在εr=9、0.025m网格上约23.5单元/波长，仅支持作为起始网格，不能据此宣称收敛。impulse激发更宽数值频带，采用范围仍限于经验证的目标频带。记录末尾官方−60dB告警、base/long频响差异都必须报告；触发告警不能自动taper或加跑。两级细化不提供渐近收敛证明。

rectangular重建可能有振铃和周期边界环绕，属于窗/显示效应，不能按视觉形态判成新地层。3.333μs周期不等于240/480ns FDTD记录，也不补回晚到信号。

该包只校准官方数值链。SFCW设备驻留时间、运动、端口/天线、接收机与标定仍未模拟；20m粉质粘土—砂岩场景将在此后单列，不属于这8次运行。
