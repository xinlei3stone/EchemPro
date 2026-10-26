# EchemPro 产品需求文档 (PRD)

| 项目 | 内容 |
|---|---|
| 产品名称 | EchemPro |
| 版本 | v1.0 (打包基线: 2026-10-11; 源码正本: 2026-10-26) |
| 文档日期 | 2026-10-26 |
| 平台 | macOS (Apple Silicon, arm64) |
| 分发形态 | 独立 .app (PyInstaller `--windowed` 单目录打包, 无需用户安装 Python) |
| 源码正本 | `EchemPro-src/` (2026-10-26 起为活代码仓库, 不再是发布快照) |
| 源码入口 | `EchemPro.py` (PyQt5 GUI, 构建 app 的入口) |
| 计算内核 | `echem_process.py` 的 `process()` (CLI 与 GUI 共用) |

---

## 1. 背景与问题

碱性 HER 电催化研究每组实验会产生多个样品、每种样品多个技术文件 (EIS / LSV / CV) 的 CS Studio 导出数据。研究人员目前需要:

1. 手工把电位换算到 RHE 标度并做 iR 补偿;
2. 对 EIS 数据做等效电路拟合取 Rs / Rct;
3. 把多个样品的结果画成对比图 (Nyquist、LSV、Tafel)。

这些步骤重复、易错 (单位、补偿系数、Tafel 窗口), 且散落在 Jupyter/脚本之间, 无法给非编程背景的合作者直接使用。

## 2. 产品目标

- **一键批处理**: 选定日期文件夹 → 勾选样品 → 填参数 → Run, 完成 EIS 拟合 + LSV 校正 + Tafel 对比全流程。
- **零环境门槛**: 双击即用, 不依赖 Python/venv/命令行。
- **结果即时可见**: 处理日志实时滚动, 生成的图直接在应用内预览。
- **可复现**: 应用输出与命令行 `echem_process.py` 逐字节一致 (同一计算内核), 满足论文数据可追溯要求。

## 3. 目标用户

| 用户 | 场景 |
|---|---|
| 电化学方向研究生 (主要) | 日常实验数据处理, 已熟悉 CS Studio 导出流程 |
| 课题组成员 (次要) | 复现/检查他人数据, 不需要看懂代码 |
| 未来的自己 (打包维护者) | 脚本更新后重建 app |

## 4. 功能需求

### 4.1 数据选择与参数配置 (P0)

- F-1 通过文本框 + "Browse…" 选择**日期文件夹** (如 `DATA/EC/20260921`); 支持直接粘贴路径, 回车刷新样品列表。
- F-2 自动列出日期文件夹下的**样品子文件夹** (仅目录, 排序), 复选框控制参与处理; 提供 Select All / Clear。
- F-3 参数表单 (默认值按本组实验惯例预设):
  - pH (默认 14.0, 范围 0–14, 步进 0.5)
  - Eref vs. SHE / V (默认 0.098, Hg/HgO)
  - iR 补偿比例 (默认 0.9, 即 90% 补偿)
  - EIS 模型: `auto` / `single` / `double` (默认 auto: 优先双时间常数 `Rs+L+(R1‖CPE1)+(R2‖CPE2)`, GoF 改善 ≤20% 时回退单 Randles)
  - Tafel Log|J| 窗口 (默认 0.8–1.6)
  - Tag (输出文件命名; 留空自动 = 去公共前缀后 `_vs_` 连接)
  - Labels (图例名, 逗号分隔; 留空用文件夹名)
  - 手动 Rs (空格分隔 Ω·cm², 数量须与样品一致; 留空用 EIS 拟合值)
- F-4 参数校验在运行前完成, 非法输入弹出 QMessageBox 警告且不启动任务: 未选文件夹、未勾样品、Labels/Rs 数量不匹配、Tafel 下限 ≥ 上限。

### 4.2 处理流程 (P0, 后台线程执行)

对每个勾选样品:

1. **文件识别**: 按文件名后缀 `-EIS` / `-LSV` 识别技术类型; `.txt` (CS Studio 导出) 优先, 其次 `.csv` (galvani 转出) 与 `.mpr` (BioLogic, 运行时转换)。
2. **EIS 拟合**: 比例加权 least-squares 等效电路拟合; auto 模式先双时间常数、不满足改善阈值回退单时间常数。输出单样品 Nyquist/Bode (实测点 + 拟合线)。
3. **LSV 校正**: 自动取回程段 (扫描顶点后); 电位换算 RHE (`E + 0.0592·pH + Eref`); 按补偿比例做 iR 补偿 (Rs 取同样品 EIS 拟合值, Ω·cm² × A/cm² 直接相乘); 保存 `*_RHE.csv`。
4. **Tafel 分析**: 在 Log|J| 窗口内线性拟合, 输出斜率与对比图。

跨样品输出: Nyquist 对比图、LSV 对比图、Tafel 对比图、EIS 拟合汇总 CSV。

### 4.3 运行反馈与结果预览 (P0)

- F-5 后台 QThread 执行, 界面不阻塞; Run 按钮处理中禁用。
- F-6 日志实时追加显示 (QPlainTextEdit, 上限 5000 行); 失败时日志末尾附 `ERROR:` 并弹窗。
- F-7 完成后右侧 Plots 标签页按序加载生成的 PNG (对比图在前、单样品图在后, 标签页显示文件名, tooltip 显示完整路径), 缩放至 960×720 预览。
- F-8 "Open Output Folder" 按钮在 Finder 中打开日期文件夹。
- F-9 状态栏显示 就绪 / 处理中… / 失败 / 完成 (tag: …, N 张图)。

### 4.4 输出物规格 (P0)

- 图: PNG, matplotlib 300 DPI, 项目配色 (最多 8 样品循环取色)。
- 数据: `*_RHE.csv` (逐样品校正后 LSV), EIS 拟合汇总 CSV (含 Rs / Rct / GoF 等)。
- 命名: 带 tag 前缀, 写入日期文件夹内, 与命令行脚本输出完全一致。

## 5. 非目标 (Non-goals, v1.0 明确不做)

- 不做 .mpr 原生解析 (galvani 未打入 app; 遇到 .mpr 仅提示转换方法)。
- 不做 CV 分析、扫速 CV 提取 (保留在独立脚本 `extract_scanrate_cv.py`)。
- 不做 DRT 反演 GUI (保留在 `drt_analysis.py` / `peis_separate_plots.py`)。
- 不做数据管理/数据库、项目管理、多实验对比归档。
- 不打包 Windows/Linux 版本 (当前仅 macOS arm64)。
- 不做自动更新; 脚本更新后由维护者手动重建。

## 6. 技术架构

源码布局 (`EchemPro-src/`, 正本):

```
EchemPro-src/
├── EchemPro.py               # PyQt5 GUI: MainWindow + PipelineThread<QThread>
├── echem_process.py          # 计算内核 process() + CLI
├── drt_analysis.py           # DRT (Tikhonov + NNLS)
├── peis_separate_plots.py    # 单 EIS 三图 (复用内核 + DRT)
├── extract_scanrate_cv.py    # 扫速 CV 最后一圈 + Ep-扫速
├── EchemPro.spec             # PyInstaller 配置
├── build.sh                  # 一键重建
├── requirements.txt          # 运行时依赖
└── PRD.md
```

`RuHeterojunctionHER/scripts/` 保留同名薄封装, 旧命令仍可用, 但实现全部在本目录。

打包产物:

```
EchemPro.app (PyInstaller, --windowed)
├── PyQt5 GUI  (EchemPro.py: MainWindow + PipelineThread<QThread>)
│     └── 信号: logged / failed / succeeded
├── 计算内核    (echem_process.py: process())
│     └── numpy / pandas / scipy.optimize / matplotlib(Agg)
└── 运行时: 内嵌 Python 3.14.3 解释器, 不读系统环境
```

- 绘图在后台线程以 Agg 后端完成 PNG, GUI 仅做 QPixmap 预览, 避免 Qt 跨线程绘图问题。
- app 是冻结产物: 源码修改后须用 `build.sh` 重建, 不会自动生效。

## 7. 数据规格 (输入契约)

| 文件 | 格式 | 说明 |
|---|---|---|
| `<样品>-CV.txt` / `<样品>-LSV.txt` | CS Studio 文本导出, 首行 `CSStudioFile` 元数据, tab 分隔, 列 `E(V) / i(A/cm²) / T(s)` | i 已是电流密度 (几何面积 1 cm²); 非数值行 (如 "N CYCLE") 自动跳过 |
| `<样品>-EIS.txt` | 同上, 表头 `Freq(Hz) / Z'(Ohm.cm²) / Z''(Ohm.cm²) / …` | Z'' 电容性为负; 列按表头名匹配, 不依赖列序; 感抗点 (Z''>0) 保留由拟合模型电感项吸收 |
| `*.mpr` | BioLogic EC-Lab | 仅提示, 不解析 (v1.0) |

约定: 跳过 macOS AppleDouble 垃圾文件 (`._` 前缀); 多样品文件时取排序后第一个。

## 8. 验收标准

1. 双击 `EchemPro.app` 在无 Python 环境的 Mac 上启动, 窗口标题为 "EchemPro"。
2. 选择含 3 个样品的日期文件夹, 默认参数 Run, 日志完整滚动, 完成后产出 ≥3 张对比图 + 每样品 Nyquist/Bode, Plots 标签页可预览。
3. 与 `../RuHeterojunctionHER/.venv/bin/python echem_process.py <同参数>` 运行, 输出文件内容一致。
4. 参数非法 (如 Tafel 下限 ≥ 上限) 时弹窗且不启动后台任务。
5. 处理期间界面可正常拖动/切换, Run 按钮禁用, 不崩溃。

(第 1 条已于 2026-10-11 构建时实测通过; 2–5 基于计算内核同源推定, 建议用真实数据回归一次。)

## 9. 已知限制与风险

| 限制 | 影响 | 缓解 |
|---|---|---|
| galvani 未打包 | 独立版无法转换 .mpr | 先用 EC-Lab/CS Studio 导出 txt |
| ad-hoc 签名 | Gatekeeper 首次需右键 → Open | 如需分发可考虑 Apple 开发者签名/公证 |
| app 体积 220 MB | 分发不便 | 可接受 (内嵌完整 Python + 科学计算栈) |
| app 是冻结产物 | 源码改进不会自动进 app | 改完运行 `./build.sh` 重建 (旧 app 自动备份); 后续可考虑 CI |
| 后台线程异常经 `failed` 信号兜底 | 极端错误可能丢栈 | 日志含 ERROR 消息, 可回命令行复现 |
| pip 脚本 shebang 指向旧路径 (`~/VibeCoding/...`) | 重建时 `pip` 命令失效 | 统一用 `.venv/bin/python -m pip` |

## 10. 路线图 (候选, 未承诺)

- v1.1: 打包 galvani (解决 .mpr); 结果表格内嵌预览 (替代仅 PNG 预览); 参数预设保存/载入。
- v1.2: 集成 DRT (复用 `drt_analysis.py`, Nyquist/Bode/DRT 三图 GUI 化); CV/扫速分析模块。
- v2.0: 实验项目管理 (多日期文件夹归档、历史结果对比); 自动更新或至少版本提示。

## 11. 运行与重建

开发运行 (复用工作一的 venv):

```bash
cd EchemPro-src
../RuHeterojunctionHER/.venv/bin/python EchemPro.py
```

重建 app:

```bash
cd EchemPro-src
./build.sh          # 或 PYTHON=/path/to/python ./build.sh
```

`build.sh` 用 PyInstaller 重新打包, 并把旧 `../EchemPro.app` 改名为带时间戳的 `.bak` 而非删除。PyInstaller 安装在 `../RuHeterojunctionHER/.venv` (6.22.3, 支持 Python 3.14)。如需独立环境, 按 `requirements.txt` 新建 venv 后 `PYTHON=... ./build.sh`。

数据默认位置不变: `../RuHeterojunctionHER/DATA/EC/YYYYMMDD/`。`extract_scanrate_cv.py` 已改为 `--root/--date/--sample/--files` 参数, 不再依赖 `__file__` 推导数据路径。
