#!/usr/bin/env python3
"""EchemPro — 电化学数据处理 PyQt5 图形界面 (DRTtools 风格卡片式布局)。

顶部导航栏 (应用名 / 状态 / Run / Open Output / About), 左侧卡片式
边栏 (Data Input / Analysis Parameters / Log), 右侧结果区 (LSV / Tafel /
Nyquist 对比图卡片网格 + 单样品图 Tabs + 空状态提示)。
计算逻辑完全复用同目录 echem_process.process, 后台线程执行, 日志实时显示。

本目录 (EchemPro-src/) 是 EchemPro 的**正本源码仓库**, 不再是打包快照。
开发运行 (复用 RuHeterojunctionHER 的 venv):
  ../RuHeterojunctionHER/.venv/bin/python EchemPro.py
重新打包 (生成 ../EchemPro.app):
  ./build.sh
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # 绘图在后台线程完成, GUI 只用 QPixmap 预览 PNG

from PyQt5.QtCore import Qt, QThread, QUrl, pyqtSignal
from PyQt5.QtGui import QDesktopServices, QPixmap
from PyQt5.QtWidgets import (
    QApplication, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout,
    QFrame, QGridLayout, QHBoxLayout, QLabel, QLineEdit, QListWidget,
    QListWidgetItem, QMainWindow, QMessageBox, QPushButton, QPlainTextEdit,
    QScrollArea, QSizePolicy, QTabWidget, QVBoxLayout, QWidget,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from echem_process import process


STYLE = """
QMainWindow, QWidget#central, QWidget#sidebarContent, QWidget#resultsContent {
    background: #f5f6f7;
}
QWidget {
    color: #1f2328;
    font-size: 13px;
}
QLabel { background: transparent; }

QFrame#navbar {
    background: #ffffff;
    border-bottom: 1px solid #e5e7eb;
}
QLabel#appName { font-size: 18px; font-weight: 700; }
QLabel#appSubtitle { color: #6b7280; font-size: 12px; }
QLabel#statusLabel { color: #6b7280; }

QFrame#card {
    background: #ffffff;
    border: 1px solid #e5e7eb;
    border-radius: 10px;
}
QLabel#cardTitle { font-size: 14px; font-weight: 600; }

QPushButton {
    background: #ffffff;
    border: 1px solid #d1d5db;
    border-radius: 6px;
    padding: 6px 14px;
}
QPushButton:hover { background: #f3f4f6; }
QPushButton:disabled { color: #9ca3af; background: #f3f4f6; }
QPushButton#primaryButton {
    background: #111827;
    color: #ffffff;
    border: none;
    border-radius: 6px;
    font-weight: 600;
    padding: 8px 24px;
}
QPushButton#primaryButton:hover { background: #1f2937; }
QPushButton#primaryButton:disabled { background: #9ca3af; }

QLineEdit, QComboBox, QDoubleSpinBox {
    background: #ffffff;
    border: 1px solid #d1d5db;
    border-radius: 6px;
    padding: 5px 8px;
    min-height: 20px;
}
QLineEdit:focus, QComboBox:focus, QDoubleSpinBox:focus {
    border: 1px solid #111827;
}
QListWidget, QPlainTextEdit {
    background: #ffffff;
    border: 1px solid #e5e7eb;
    border-radius: 6px;
}
QPlainTextEdit { font-family: Menlo, monospace; font-size: 12px; }

QTabWidget::pane {
    border: 1px solid #e5e7eb;
    border-radius: 6px;
    background: #ffffff;
}
QTabBar::tab { padding: 6px 14px; color: #6b7280; }
QTabBar::tab:selected { color: #111827; font-weight: 600; }

QLabel#emptyState { color: #9ca3af; font-size: 16px; }
QLabel#emptyHint { color: #c3c8cf; font-size: 12px; }

QScrollArea { border: none; background: transparent; }
"""


class PipelineThread(QThread):
    """后台执行 process(), 通过信号回传日志 / 结果 / 错误"""

    logged = pyqtSignal(str)
    failed = pyqtSignal(str)
    succeeded = pyqtSignal(dict)

    def __init__(self, params, parent=None):
        super().__init__(parent)
        self.params = params

    def run(self):
        try:
            result = process(log=self.logged.emit, **self.params)
        except Exception as exc:  # ValueError (参数) 与运行期错误都在这里
            self.failed.emit(str(exc))
        else:
            self.succeeded.emit(result)


class PlotLabel(QLabel):
    """随窗口尺寸平滑缩放的 PNG 预览标签"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._pix = None
        self.setAlignment(Qt.AlignCenter)
        self.setMinimumSize(240, 180)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

    def set_plot(self, path):
        self._pix = QPixmap(str(path))
        self.setToolTip(str(path))
        self._rescale()

    def clear_plot(self):
        self._pix = None
        self.clear()

    def resizeEvent(self, event):
        self._rescale()
        super().resizeEvent(event)

    def _rescale(self):
        if self._pix is not None and not self._pix.isNull():
            self.setPixmap(self._pix.scaled(
                self.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation))


def make_card(title):
    """白底圆角卡片: 标题行 + 内容布局"""
    card = QFrame()
    card.setObjectName("card")
    layout = QVBoxLayout(card)
    layout.setContentsMargins(16, 14, 16, 16)
    layout.setSpacing(10)
    header = QLabel(title)
    header.setObjectName("cardTitle")
    layout.addWidget(header)
    return card, layout


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("EchemPro — EIS / LSV / Tafel Batch Processor")
        self.resize(1360, 880)
        self.thread = None

        central = QWidget()
        central.setObjectName("central")
        root_layout = QVBoxLayout(central)
        root_layout.setContentsMargins(0, 0, 0, 0)
        root_layout.setSpacing(0)
        root_layout.addWidget(self._build_navbar())

        body = QHBoxLayout()
        body.setContentsMargins(16, 16, 16, 16)
        body.setSpacing(16)
        body.addWidget(self._build_sidebar())
        body.addWidget(self._build_results(), 1)
        root_layout.addLayout(body, 1)

        self.setCentralWidget(central)

    # ---------- 顶部导航栏 ----------

    def _build_navbar(self):
        bar = QFrame()
        bar.setObjectName("navbar")
        bar.setFixedHeight(64)
        layout = QHBoxLayout(bar)
        layout.setContentsMargins(20, 0, 20, 0)
        layout.setSpacing(10)

        title_box = QVBoxLayout()
        title_box.setSpacing(0)
        name = QLabel("EchemPro")
        name.setObjectName("appName")
        subtitle = QLabel("EIS · LSV · Tafel batch processor")
        subtitle.setObjectName("appSubtitle")
        title_box.addWidget(name)
        title_box.addWidget(subtitle)
        layout.addLayout(title_box)
        layout.addStretch(1)

        self.status_label = QLabel("Ready")
        self.status_label.setObjectName("statusLabel")
        layout.addWidget(self.status_label)

        open_btn = QPushButton("Open Output")
        open_btn.clicked.connect(self.open_output_folder)
        layout.addWidget(open_btn)

        about_btn = QPushButton("About")
        about_btn.clicked.connect(self.show_about)
        layout.addWidget(about_btn)

        self.run_button = QPushButton("Run")
        self.run_button.setObjectName("primaryButton")
        self.run_button.setMinimumHeight(38)
        self.run_button.clicked.connect(self.on_run)
        layout.addWidget(self.run_button)
        return bar

    # ---------- 左侧边栏 ----------

    def _build_sidebar(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFixedWidth(420)
        content = QWidget()
        content.setObjectName("sidebarContent")
        layout = QVBoxLayout(content)
        layout.setContentsMargins(0, 0, 4, 0)
        layout.setSpacing(16)

        # --- Data Input ---
        data_card, data_layout = make_card("Data Input")
        choose_row = QHBoxLayout()
        choose_btn = QPushButton("Choose Folder")
        choose_btn.clicked.connect(self.browse_root)
        choose_row.addWidget(choose_btn)
        self.root_edit = QLineEdit()
        self.root_edit.setPlaceholderText("e.g. DATA/EC/20260921")
        self.root_edit.returnPressed.connect(self.refresh_samples)
        choose_row.addWidget(self.root_edit, 1)
        data_layout.addLayout(choose_row)

        data_layout.addWidget(QLabel("Sample folders (check to include)"))
        self.sample_list = QListWidget()
        self.sample_list.setMaximumHeight(150)
        data_layout.addWidget(self.sample_list)
        sel_row = QHBoxLayout()
        select_all = QPushButton("Select All")
        select_all.clicked.connect(lambda: self.set_all_checked(True))
        select_none = QPushButton("Clear")
        select_none.clicked.connect(lambda: self.set_all_checked(False))
        sel_row.addWidget(select_all)
        sel_row.addWidget(select_none)
        sel_row.addStretch(1)
        data_layout.addLayout(sel_row)
        layout.addWidget(data_card)

        # --- Analysis Parameters ---
        param_card, param_layout = make_card("Analysis Parameters")
        form = QFormLayout()
        form.setSpacing(8)

        self.ph_spin = QDoubleSpinBox(minimum=0, maximum=14, value=14.0,
                                      decimals=1, singleStep=0.5)
        form.addRow("pH", self.ph_spin)

        self.eref_spin = QDoubleSpinBox(minimum=-1, maximum=1, decimals=3,
                                        value=0.098, singleStep=0.001)
        form.addRow("Eref vs. SHE (V)", self.eref_spin)

        self.ir_spin = QDoubleSpinBox(minimum=0, maximum=1, value=0.9,
                                      decimals=2, singleStep=0.05)
        form.addRow("iR compensation", self.ir_spin)

        self.model_combo = QComboBox()
        self.model_combo.addItems(["auto", "single", "double"])
        form.addRow("EIS model", self.model_combo)

        tafel_row = QHBoxLayout()
        self.tafel_lo = QDoubleSpinBox(minimum=-3, maximum=3, value=0.8,
                                       decimals=2, singleStep=0.1)
        self.tafel_hi = QDoubleSpinBox(minimum=-3, maximum=3, value=1.6,
                                       decimals=2, singleStep=0.1)
        tafel_row.addWidget(self.tafel_lo)
        tafel_row.addWidget(QLabel("to"))
        tafel_row.addWidget(self.tafel_hi)
        tafel_row.addStretch(1)
        form.addRow("Tafel Log|J| window", tafel_row)

        self.tag_edit = QLineEdit()
        self.tag_edit.setPlaceholderText(
            "Empty = auto (strip common prefix, join with _vs_)")
        form.addRow("Tag", self.tag_edit)

        self.labels_edit = QLineEdit()
        self.labels_edit.setPlaceholderText(
            "Comma-separated; count must match samples. Empty = folder names")
        form.addRow("Labels", self.labels_edit)

        self.rs_edit = QLineEdit()
        self.rs_edit.setPlaceholderText(
            "Space-separated Ω·cm². Empty = use EIS-fitted Rs")
        form.addRow("Manual Rs", self.rs_edit)

        param_layout.addLayout(form)
        layout.addWidget(param_card)

        # --- Log ---
        log_card, log_layout = make_card("Log")
        self.log_view = QPlainTextEdit(readOnly=True)
        self.log_view.setMaximumBlockCount(5000)
        self.log_view.setMinimumHeight(140)
        log_layout.addWidget(self.log_view)
        layout.addWidget(log_card, 1)

        scroll.setWidget(content)
        return scroll

    # ---------- 右侧结果区 ----------

    def _build_results(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        content = QWidget()
        content.setObjectName("resultsContent")
        layout = QVBoxLayout(content)
        layout.setContentsMargins(0, 0, 4, 0)
        layout.setSpacing(16)

        self.empty_container = QWidget()
        self.empty_container.setObjectName("resultsContent")
        empty_layout = QVBoxLayout(self.empty_container)
        empty_layout.setContentsMargins(0, 0, 0, 0)
        self.empty_state = QLabel("Run analysis to see results")
        self.empty_state.setObjectName("emptyState")
        self.empty_state.setAlignment(Qt.AlignCenter)
        self.empty_hint = QLabel("Choose a date folder, check samples, then press Run")
        self.empty_hint.setObjectName("emptyHint")
        self.empty_hint.setAlignment(Qt.AlignCenter)
        empty_layout.addStretch(1)
        empty_layout.addWidget(self.empty_state)
        empty_layout.addWidget(self.empty_hint)
        empty_layout.addStretch(1)
        layout.addWidget(self.empty_container, 1)

        self.results_widget = QWidget()
        self.results_widget.setObjectName("resultsContent")
        results_layout = QVBoxLayout(self.results_widget)
        results_layout.setContentsMargins(0, 0, 0, 0)
        results_layout.setSpacing(16)

        grid = QGridLayout()
        grid.setSpacing(16)
        self.lsv_card, lsv_layout = make_card("LSV Comparison")
        self.lsv_plot = PlotLabel()
        self.lsv_plot.setMinimumHeight(300)
        lsv_layout.addWidget(self.lsv_plot)
        grid.addWidget(self.lsv_card, 0, 0)

        self.tafel_card, tafel_layout = make_card("Tafel Comparison")
        self.tafel_plot = PlotLabel()
        self.tafel_plot.setMinimumHeight(300)
        tafel_layout.addWidget(self.tafel_plot)
        grid.addWidget(self.tafel_card, 0, 1)

        self.nyquist_card, nyquist_layout = make_card("Nyquist Comparison")
        self.nyquist_plot = PlotLabel()
        self.nyquist_plot.setMinimumHeight(320)
        nyquist_layout.addWidget(self.nyquist_plot)
        grid.addWidget(self.nyquist_card, 1, 0, 1, 2)
        results_layout.addLayout(grid)

        self.samples_card, samples_layout = make_card("Per-sample Plots")
        self.plot_tabs = QTabWidget()
        self.plot_tabs.setMinimumHeight(360)
        samples_layout.addWidget(self.plot_tabs)
        results_layout.addWidget(self.samples_card)

        self.results_widget.hide()
        layout.addWidget(self.results_widget, 1)

        scroll.setWidget(content)
        return scroll

    # ---------- 文件夹 / 样品列表 ----------

    def browse_root(self):
        start = self.root_edit.text().strip() or "DATA/EC"
        path = QFileDialog.getExistingDirectory(self, "Choose date folder", start)
        if path:
            self.root_edit.setText(path)
            self.refresh_samples()

    def refresh_samples(self):
        self.sample_list.clear()
        root = Path(self.root_edit.text().strip())
        if not root.is_dir():
            return
        for d in sorted(p for p in root.iterdir() if p.is_dir()):
            item = QListWidgetItem(d.name)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked)
            self.sample_list.addItem(item)

    def set_all_checked(self, checked):
        state = Qt.Checked if checked else Qt.Unchecked
        for i in range(self.sample_list.count()):
            self.sample_list.item(i).setCheckState(state)

    def checked_samples(self):
        return [self.sample_list.item(i).text()
                for i in range(self.sample_list.count())
                if self.sample_list.item(i).checkState() == Qt.Checked]

    # ---------- 运行 ----------

    def collect_params(self):
        """收集界面参数; 非法输入抛 ValueError"""
        root_text = self.root_edit.text().strip()
        if not root_text:
            raise ValueError("Please choose a date folder first")
        root = Path(root_text)
        if not root.is_dir():
            raise ValueError(f"not a directory: {root}")
        names = self.checked_samples()
        if not names:
            raise ValueError("Please check at least one sample folder")

        labels_text = self.labels_edit.text().strip()
        labels = [s.strip() for s in labels_text.split(",")] if labels_text else None

        rs_text = self.rs_edit.text().strip()
        if rs_text:
            try:
                rs_manual = [float(v) for v in rs_text.split()]
            except ValueError:
                raise ValueError(f"Manual Rs contains invalid numbers: {rs_text!r}")
        else:
            rs_manual = None

        tag = self.tag_edit.text().strip() or None
        lo, hi = self.tafel_lo.value(), self.tafel_hi.value()
        if lo >= hi:
            raise ValueError("Tafel window lower bound must be below upper bound")

        return dict(root=root, folders=[root / n for n in names],
                    labels=labels, tag=tag, ph=self.ph_spin.value(),
                    eref=self.eref_spin.value(), ir_comp=self.ir_spin.value(),
                    rs_manual=rs_manual, tafel_window=(lo, hi),
                    eis_model=self.model_combo.currentText())

    def on_run(self):
        try:
            params = self.collect_params()
        except ValueError as e:
            QMessageBox.warning(self, "Invalid parameters", str(e))
            return
        self.log_view.clear()
        self.run_button.setEnabled(False)
        self.status_label.setText("Processing…")
        self.thread = PipelineThread(params)
        self.thread.logged.connect(self.log_view.appendPlainText)
        self.thread.failed.connect(self.on_failed)
        self.thread.succeeded.connect(self.on_succeeded)
        self.thread.start()

    def on_failed(self, message):
        self.run_button.setEnabled(True)
        self.status_label.setText("Failed")
        self.log_view.appendPlainText(f"ERROR: {message}")
        QMessageBox.warning(self, "Processing failed", message)

    def on_succeeded(self, result):
        self.run_button.setEnabled(True)
        self.status_label.setText(f"Done (tag: {result['tag']})")
        self.show_outputs(result["outputs"])

    # ---------- 结果预览 ----------

    def show_outputs(self, outputs):
        pngs = [Path(p) for p in outputs
                if Path(p).suffix.lower() == ".png" and Path(p).exists()]

        mapping = [("LSV_comparison", self.lsv_plot, self.lsv_card),
                   ("Tafel_comparison", self.tafel_plot, self.tafel_card),
                   ("Nyquist_comparison", self.nyquist_plot, self.nyquist_card)]
        for _, plot, card in mapping:
            plot.clear_plot()
            card.hide()
        self.plot_tabs.clear()

        extra = []
        for p in pngs:
            for key, plot, card in mapping:
                if key in p.name:
                    plot.set_plot(p)
                    card.show()
                    break
            else:
                extra.append(p)
        for p in sorted(extra):
            label = PlotLabel()
            label.set_plot(p)
            index = self.plot_tabs.addTab(label, p.stem)
            self.plot_tabs.setTabToolTip(index, str(p))
        self.samples_card.setVisible(bool(extra))

        self.results_widget.setVisible(bool(pngs))
        self.empty_container.setVisible(not pngs)
        if pngs:
            self.status_label.setText(
                self.status_label.text() + f" — {len(pngs)} plots")

    def open_output_folder(self):
        root_text = self.root_edit.text().strip()
        if root_text and Path(root_text).is_dir():
            QDesktopServices.openUrl(QUrl.fromLocalFile(root_text))

    def show_about(self):
        QMessageBox.about(
            self, "About EchemPro",
            "<b>EchemPro</b> — EIS / LSV / Tafel batch processor<br><br>"
            "Input: CS Studio text exports in per-sample folders "
            "(<code>*-EIS.txt</code>, <code>*-LSV.txt</code>).<br>"
            "Pipeline: equivalent-circuit EIS fitting (auto single/double "
            "time constant), LSV RHE conversion with iR compensation, "
            "and Tafel analysis, with comparison plots across samples.<br><br>"
            "Computation: <code>echem_process.py</code>")


def main():
    app = QApplication(sys.argv)
    app.setApplicationName("EchemPro")
    app.setStyleSheet(STYLE)
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
