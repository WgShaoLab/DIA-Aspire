import sys
import re
import os
import subprocess
import tarfile
import requests
import json
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QDesktopWidget,QButtonGroup,
                             QLabel, QLineEdit, QPushButton, QComboBox, QGroupBox, QGridLayout,
                             QListWidget, QMessageBox, QTextEdit, QFileDialog, QCheckBox, QRadioButton,
                             QMenu,QCompleter, QSizePolicy,QProgressDialog)
from PyQt5.QtCore import (Qt, QStringListModel, QProcess, QEvent, QTimer, 
                          QThread, pyqtSignal)

# # 将src目录添加到导入路径
# src_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "src")
# sys.path.append(src_dir)

from dia_aspire.converters import fragpipe_api
from dia_aspire.converters import systemhc_api
from dia_aspire.converters import sptxt2tsv
from dia_aspire.converters import irt_alignment
from dia_aspire.converters import maxquant_api
from dia_aspire.converters import msms2tsv

# 在类定义之前添加
# import io
from pathlib import Path
from importlib import resources
from contextlib import redirect_stdout, redirect_stderr
from datetime import datetime
from alpharaw import register_all_readers

import warnings
warnings.filterwarnings("ignore")
import logging

logging.basicConfig(level=logging.INFO)
register_all_readers()

try:
    # Python 3.9+
    from importlib import resources as importlib_resources
except ImportError:
    # Python < 3.9
    import importlib_resources


class TeeLogger:
    def __init__(self, output_widget, log_file_path):
        self.output_widget = output_widget
        self.log_file = open(log_file_path, 'a', encoding='utf-8')

        self.original_stdout = sys.stdout
        self.original_stderr = sys.stderr

    def write(self, message):
        # 处理可能是 bytes 的情况
        if isinstance(message, bytes):
            try:
                message = message.decode('utf-8')
            except Exception:
                message = str(message)
        
        # 忽略空白
        if not message or not message.strip():
            return

        # 显示到 GUI（QTextEdit.append 接受 str）
        if self.output_widget and hasattr(self.output_widget, 'append'):
            lines = message.rstrip().split('\n')
            for line in lines:
                if line.strip():  # 只输出非空行
                    self.output_widget.append(line)
                    # 实时刷新 GUI
                    QApplication.processEvents()

        # 写入日志文件（确保未关闭）
        if self.log_file and not self.log_file.closed:
            self.log_file.write(message)
            self.log_file.flush()

    def flush(self):
        if self.log_file and not self.log_file.closed:
            self.log_file.flush()
        self.original_stdout.flush()
        self.original_stderr.flush()

    def close(self):
        if self.log_file and not self.log_file.closed:
            self.log_file.close()

    def fileno(self):
        return self.log_file.fileno()
    
    def isatty(self):
        return False

class QTextEditHandler(logging.Handler):
    """将 logging 输出重定向到 TeeLogger"""
    def __init__(self, tee_logger):
        super().__init__()
        self.tee_logger = tee_logger

    def emit(self, record):
        try:
            msg = self.format(record)
            if not msg.endswith('\n'):
                msg += '\n'
            self.tee_logger.write(msg)

        except Exception:
            self.handleError(record)

class SingleDownloadThread(QThread):
    download_done = pyqtSignal(bool, str, str)  # success, message, file_path

    def __init__(self, allele, selected_class, output_dir):
        super().__init__()
        self.allele = allele
        self.selected_class = selected_class
        self.output_dir = output_dir

    def run(self):
        file_name = f"HCD_cons_{self.allele}_top12_bynam_ptm.tsv"
        url = f"https://systemhc.sjtu.edu.cn/data/Systemhc_v2_2023/Data/230623_build/SysteMHC_Library/{self.selected_class}/Allele-specific/{file_name}"
        dest_path = os.path.join(self.output_dir, file_name)

        try:
            response = requests.get(url, stream=True)
            if response.status_code != 200:
                raise Exception(f"HTTP status code: {response.status_code}")
            
            with open(dest_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            
            self.download_done.emit(True, f"Downloaded: {file_name}", self.allele)
        except Exception as e:
            self.download_done.emit(False, f"Error: {str(e)}", self.allele)

class CommandLineGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.parameters = {
                    'sample_library_files': None,
                    'systemhc_library_files': None,  # New list for SysteMHC libraries
                    'extra_params': {}
                }
        self.rescore_enabled = True
        self.rescore_model = "DNN"

        self.default_params = {
            'threads': '16',
            'qvalue': '1',
            'matrix-qvalue': '0.01',
            'matrices':'true',
            'report-decoys': 'true',

            'mass-acc': '15',
            'mass-acc-ms1': '15',
            'export-quant': 'true',
            'xic': '30',
            'xic-theoretical-fr': 'false',

            'reanalyse':'false',
            'no-prot-inf': 'true',
            'rt-profiling': 'true',
            'pg-level': '1',
            'verbose':'4',
        }

        self.batch_download_attrs = {
            'alleles': [],
            'selected_class': '',
            'output_dir': '',
            'progress': None,
            'success_count': 0,
            'failed_alleles': [],
            'current_allele_idx': 0,
            'download_threads': [],
        }

        self.allele_list = {}
        self.load_allele_list()
        self.process = QProcess(self)
        self.main_layout = QVBoxLayout()
        self.extra_params_widget = None
        self.selected_pipeline = "SysteMHC-pipeline"  # Default pipeline
        self.library_filter = "SPTXT Files (*.sptxt)"
        self.log_file = None  # 日志文件句柄
        self.main_logger = None  # 初始化logger
        self.initUI()

        # Connect signals
        self.process.readyReadStandardOutput.connect(self.handle_output)
        self.process.readyReadStandardError.connect(self.handle_error)
        self.process.finished.connect(self.task_finished)

    def load_allele_list(self):
        """Load allele list from JSON file"""
        #json_path = "constants/allele_list.json"  # 或你的路径
        with importlib_resources.files('dia_aspire.constants').joinpath('allele_list.json') as json_path:
            with open(json_path, "r") as f:
                self.allele_list = json.load(f)

    def initUI(self):
        self.setWindowTitle('DIA-Aspire')
        self.setGeometry(300, 300, 800, 600)
        self.center_window()
        
        # 将logo路径和输入类型选择放在同一水平布局中
        # top_layout = QHBoxLayout()
        logo_io_layout = QHBoxLayout()
        logo_io_layout.setSpacing(10)  # Logo与输入项的间距
        logo_io_layout.setContentsMargins(10, 10, 10, 10)  # 整体边距

        # 1. Logo 部分（垂直拉伸匹配输入项高度）
        logo_widget = QWidget()
        logo_layout = QVBoxLayout(logo_widget)
        logo_layout.setContentsMargins(0, 0, 0, 0)

        # 设置应用程序logo
        # 创建logo标签并放在左侧
        logo_label = QLabel()
        #logo_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "constants", "DIA-Aspire_logo.jpg")
        with importlib_resources.files('dia_aspire.constants').joinpath('DIA-Aspire_logo.jpg') as logo_path:
            if logo_path.exists():
                from PyQt5.QtGui import QIcon, QPixmap
                self.setWindowIcon(QIcon(str(logo_path)))
                pixmap = QPixmap(str(logo_path))
                scaled_pixmap = pixmap.scaled(120, pixmap.height(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
                logo_label.setPixmap(scaled_pixmap)
                logo_label.setScaledContents(False)  # 关闭自适应容器大小（避免Logo被拉伸）
                logo_label.setAlignment(Qt.AlignCenter)  # Logo在容器内居中
        
        # 关键：让Logo垂直拉伸至与输入项总高度一致
        from PyQt5.QtWidgets import QSizePolicy 
        logo_label.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        logo_layout.setAlignment(Qt.AlignCenter)
        logo_layout.addWidget(logo_label)
        
        logo_io_layout.addWidget(logo_widget)  # Logo加入水平布局

        # 2. IO 输入项部分（4行，标签左对齐）
        io_widget = QWidget()
        self.create_io_layout(io_widget)  # 传入IO容器，重构IO布局
        logo_io_layout.addWidget(io_widget, stretch=1)  # 输入项占剩余空间

        # 将 Logo+IO 布局加入主布局
        self.main_layout.addLayout(logo_io_layout)
        
        # 以下代码保持不变
        # Create pipeline selector
        self.create_pipeline_selector()
        
        # Create library input sections
        self.create_parameter_inputs()
        
        # Create Allele Specific section
        self.create_allele_specific_layout()

        # Create additional parameters section
        toggle_extra_params_btn = QPushButton("DIANN Parameters")
        toggle_extra_params_btn.clicked.connect(self.toggle_extra_params_section)
        self.main_layout.addWidget(toggle_extra_params_btn)
        self.create_extra_params_section()

        # Action buttons
        button_layout = QHBoxLayout()
        button_layout.setContentsMargins(0, 0, 0, 0)  # 每行也不留边距
        # button_layout.setSpacing(20)

        
        # ==Rescore 配置 ==
        rescore_layout = QHBoxLayout()
        rescore_layout.setContentsMargins(0, 0, 0, 0)  # 每行也不留边距
        
        self.rescore_checkbox = QCheckBox("Post-processing (rescore)")
        self.rescore_checkbox.setChecked(False)

        # 模型标签
        model_label = QLabel("Model:")

        # 模型选择（初始禁用）
        self.rescore_svm_radio = QRadioButton("SVM")
        self.rescore_dnn_radio = QRadioButton("DNN")
        self.rescore_dnn_radio.setChecked(True)

        self.rescore_model_group = QButtonGroup(self)
        self.rescore_model_group.addButton(self.rescore_svm_radio)
        self.rescore_model_group.addButton(self.rescore_dnn_radio)

        # 初始禁用（包括按钮和标签）
        for w in [model_label, self.rescore_svm_radio, self.rescore_dnn_radio]:
            w.setEnabled(False)

        # 连接复选框状态，控制模型控件启用/禁用
        def toggle_rescore_widgets(state):
            enabled = (state == Qt.Checked)
            model_label.setEnabled(enabled)
            self.rescore_svm_radio.setEnabled(enabled)
            self.rescore_dnn_radio.setEnabled(enabled)

        self.rescore_checkbox.stateChanged.connect(toggle_rescore_widgets)

        # 添加到布局
        rescore_layout.addWidget(self.rescore_checkbox)
        rescore_layout.addSpacing(20)  # 可选：增加一点间距
        rescore_layout.addWidget(model_label)
        rescore_layout.addWidget(self.rescore_svm_radio)
        rescore_layout.addWidget(self.rescore_dnn_radio)
        rescore_layout.addStretch(1)  # 剩余空间右对齐（可选）

        # 将整行添加到主布局
        self.main_layout.addLayout(rescore_layout)


        # ----- Run & Stop -------
        self.run_btn = QPushButton('Run')
        self.run_btn.clicked.connect(self.execute_command)
        self.stop_btn = QPushButton('Stop')
        self.stop_btn.clicked.connect(self.stop_task)
        button_layout.addWidget(self.run_btn)
        button_layout.addWidget(self.stop_btn)

        # Output area
        self.output_area = QTextEdit()
        self.output_area.setReadOnly(True)
        self.output_area.setStyleSheet("""
        QTextEdit {
            background-color: #DCDCDC;
            color: #000000;
            font-family: Consolas;
            font-size: 12pt;
            #border: 2px solid #FFFFFF;
            padding: 5px;
        }
        """)

        # 设置垂直策略为扩展，以便在垂直方向上占据更多空间
        self.output_area.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        # 设置最小高度
        self.output_area.setMinimumHeight(200)
        self.main_layout.addLayout(button_layout)
        self.main_layout.addWidget(self.output_area)
        self.setLayout(self.main_layout)

    def center_window(self):
        screen_geometry = QDesktopWidget().screenGeometry()
        x = (screen_geometry.width() - self.width()) // 2
        y = (screen_geometry.height() - self.height()) // 2
        self.move(x, y)

    def create_io_layout(self,parent_widget):
        """Create input and output path layout"""

        io_layout = QVBoxLayout(parent_widget)
        io_layout.setContentsMargins(0, 0, 0, 0)  # 统一内边距
        io_layout.setSpacing(8)  # 行间距
        
        # ========== 1. DIA Input 行 ==========
        dia_input_layout = QHBoxLayout()
        dia_input_layout.setContentsMargins(0, 0, 0, 0)  # 每行也不留边距
        # dia_input_layout.setSpacing(20)
        dia_label = QLabel("DIA Input (raw, mzML, diaPASEF.d):")
        dia_label.setWordWrap(False)  # 关键：开启自动换行
        # dia_label.setFixedWidth(60)  # 宽度
        dia_label.setAlignment(Qt.AlignVCenter)  # 垂直居中，适配换行后的高度
        dia_input_layout.addWidget(dia_label)
        
        self.input_type_combo = QComboBox()
        self.input_type_combo.addItems(["Folder", "Files"])
        self.input_type_combo.currentTextChanged.connect(self.toggle_input_type)
        # self.input_type_combo.setMinimumWidth(150)  # 最小宽度，避免挤压
        dia_input_layout.addWidget(self.input_type_combo)
        io_layout.addLayout(dia_input_layout)


        LABEL_WIDTH = 120
        # ========== 2. Input Directory 行 ==========
        # input_dir_layout =  QHBoxLayout()
        self.input_dir_widget = QWidget()
        input_dir_layout = QHBoxLayout(self.input_dir_widget)
        input_dir_layout.setContentsMargins(0, 0, 0, 0)
        
        input_dir_label = QLabel("Input Directory:")
        input_dir_label.setFixedWidth(LABEL_WIDTH)  # 统一宽度
        input_dir_label.setAlignment(Qt.AlignVCenter)
        input_dir_layout.addWidget(input_dir_label)
        
        self.input_dir_field = QLineEdit()
        btn_browse_input = QPushButton("Browse")
        btn_browse_input.clicked.connect(self.select_input)
        input_dir_layout.addWidget(self.input_dir_field)
        input_dir_layout.addWidget(btn_browse_input)
        io_layout.addWidget(self.input_dir_widget)
        # io_layout.addLayout(input_dir_layout)

        # ========== 3. Output Directory 行 ==========
        output_dir_layout = QHBoxLayout()
        output_dir_layout.setContentsMargins(0, 0, 0, 0)
        output_dir_label = QLabel("Output Directory:")
        output_dir_label.setFixedWidth(LABEL_WIDTH)  # 统一宽度
        output_dir_layout.addWidget(output_dir_label)
        
        self.output_dir_field = QLineEdit()
        btn_browse_output = QPushButton("Browse")
        btn_browse_output.clicked.connect(self.select_output_dir)
        output_dir_layout.addWidget(self.output_dir_field)
        output_dir_layout.addWidget(btn_browse_output)
        io_layout.addLayout(output_dir_layout)

        # ========== 4. DIANN Software Path 行 ==========
        diann_path_layout = QHBoxLayout()
        diann_path_layout.setContentsMargins(0, 0, 0, 0)
        diann_label = QLabel("DIANN Path:")
        diann_label.setFixedWidth(LABEL_WIDTH)  # 统一宽度
        diann_path_layout.addWidget(diann_label)
        
        self.diann_path_field = QLineEdit()
        self.diann_path_field.setText("/data/xhuang/DIA-NN-2.0-Academia-Linux/diann-2.0/diann-linux")
        btn_browse_diann = QPushButton("Browse")
        btn_browse_diann.clicked.connect(self.select_diann_path)
        diann_path_layout.addWidget(self.diann_path_field)
        diann_path_layout.addWidget(btn_browse_diann)
        io_layout.addLayout(diann_path_layout)

        # ========== Input Files 行（默认隐藏） ==========
        self.input_files_widget = QWidget()
        input_files_layout = QVBoxLayout(self.input_files_widget)
        input_files_layout.setContentsMargins(0, 0, 0, 0)
        
        input_files_label = QLabel("Input Files:")
        input_files_label.setFixedWidth(LABEL_WIDTH)
        input_files_header = QHBoxLayout()
        input_files_header.addWidget(input_files_label)
        input_files_header.addStretch()# 左侧添加拉伸，让按钮整体居中
        input_files_layout.addLayout(input_files_header)
        
        self.input_files_list = QListWidget()
        self.input_files_list.setSelectionMode(QListWidget.ExtendedSelection)  # 支持多选
        # 键盘删除/退格键删除
        self.input_files_list.keyPressEvent = lambda event: self.handle_key_press(
            event, self.input_files_list, self.remove_input_files)
        # 右键菜单删除
        self.input_files_list.setContextMenuPolicy(Qt.CustomContextMenu)
        self.input_files_list.customContextMenuRequested.connect(
            lambda pos: self.show_context_menu(pos, self.input_files_list, self.remove_input_files))
        # 拖放重排（可选）
        self.input_files_list.setDragDropMode(QListWidget.InternalMove)
        input_files_layout.addWidget(self.input_files_list)

        # 按钮布局：Add + Remove 水平居中分布
        input_files_buttons_layout = QHBoxLayout()
        input_files_buttons_layout.addStretch()  # 左侧添加拉伸，让按钮整体居中
        btn_add_file = QPushButton("Add Files")
        btn_add_file.clicked.connect(self.add_input_files)
        btn_remove_file = QPushButton("Remove Files")  # 新增删除按钮
        btn_remove_file.clicked.connect(self.remove_input_files)  # 绑定删除逻辑
        input_files_buttons_layout.addWidget(btn_add_file)
        input_files_buttons_layout.addSpacing(10) 
        input_files_buttons_layout.addWidget(btn_remove_file)
        input_files_buttons_layout.addStretch()  # 右侧留白
        input_files_layout.addLayout(input_files_buttons_layout)  # 替换原单独的Add按钮
        
        io_layout.addWidget(self.input_files_widget)
        # 默认隐藏 Input Files
        self.input_files_widget.setVisible(False)
 
    def select_diann_path(self):
        """Select DIANN software path"""
        file_path, _ = QFileDialog.getOpenFileName(self, "Select DIANN Executable")
        if file_path:
            self.diann_path_field.setText(file_path)
        
    def create_pipeline_selector(self):
        """Create pipeline selector (FragPipe, SysteMHC or MaxQuant)"""
        pipeline_group = QGroupBox("Pipeline Used In Sample Library Generation")
        pipeline_layout = QHBoxLayout()
        pipeline_layout.setContentsMargins(0, 0, 0, 0)  # 每行也不留边距
        pipeline_layout.setSpacing(10)

        self.fragpipe_radio = QRadioButton("FragPipe")
        self.systemhc_radio = QRadioButton("SysteMHC-pipeline")
        self.maxquant_radio = QRadioButton("MaxQuant")  # 添加MaxQuant选项
        self.systemhc_radio.setChecked(True)  # Default to systemhc
        self.on_pipeline_changed()

        self.systemhc_radio.toggled.connect(self.on_pipeline_changed)
        self.fragpipe_radio.toggled.connect(self.on_pipeline_changed)
        self.maxquant_radio.toggled.connect(self.on_pipeline_changed)  # 连接新选项
        
        pipeline_layout.addWidget(self.fragpipe_radio)
        pipeline_layout.addWidget(self.systemhc_radio)
        pipeline_layout.addWidget(self.maxquant_radio)  # 添加到布局
        pipeline_group.setLayout(pipeline_layout)
        
        self.main_layout.addWidget(pipeline_group)
    
    def on_pipeline_changed(self):
        """Handle pipeline selection change"""
        if self.fragpipe_radio.isChecked():
            self.selected_pipeline = "FragPipe"
            # self.sample_lib_label.setText("Sample Library (tsv files):")
            self.library_filter = "Text Files (*.tsv)"

        elif self.systemhc_radio.isChecked():
            self.selected_pipeline = "SysteMHC"
            # self.sample_lib_label.setText("Sample Library (sptxt files):")
            self.library_filter = "SPTXT Files (*.sptxt)"
        else:  # MaxQuant
            self.selected_pipeline = "MaxQuant"
            # self.sample_lib_label.setText("Sample Library (msms.txt files):")
            self.library_filter = "MaxQuant Files (*.txt)"
    # def on_pipeline_changed(self, button):
    #     if button == self.fragpipe_radio:
    #         # self.expected_library_ext = "*.txt"      # 或 msms.txt
    #         self.library_filter = "Text Files (*.tsv)"
    #     elif button == self.systemhc_radio:
    #         # self.expected_library_ext = "*.sptxt"
    #         self.library_filter = "SPTXT Files (*.sptxt)"
    #     elif button == self.maxquant_radio:
    #         # self.expected_library_ext = "*.txt"
    #         self.library_filter = "MaxQuant Files (*.txt)"
    #     else:
    #         # self.expected_library_ext = "*.*"
    #         self.library_filter = "All Files (*)"

    def create_parameter_inputs(self):
        """Create library input sections with remove functionality"""
        param_group = QGroupBox("Library Input")
        grid = QGridLayout()
        grid.setContentsMargins(0, 4, 0, 4)  # ← 关键！移除 grid 自身的上下左右 margin
        grid.setSpacing(8)  # 可选：减小行/列间距（默认可能为 10~15）

        # Sample Library section
        self.sample_lib_label = QLabel("Sample Library:")
        grid.addWidget(self.sample_lib_label, 0, 0)
        
        self.parameters['sample_library_files'] = QListWidget()
        # 启用多选功能
        self.parameters['sample_library_files'].setSelectionMode(QListWidget.ExtendedSelection)
        # 连接键盘删除键
        self.parameters['sample_library_files'].keyPressEvent = lambda event: self.handle_key_press(
            event, self.parameters['sample_library_files'], self.remove_sample_library_files)
        # 右键菜单
        self.parameters['sample_library_files'].setContextMenuPolicy(Qt.CustomContextMenu)
        self.parameters['sample_library_files'].customContextMenuRequested.connect(
            lambda pos: self.show_context_menu(pos, self.parameters['sample_library_files'], self.remove_sample_library_files))
        # 启用拖放重排
        self.parameters['sample_library_files'].setDragDropMode(QListWidget.InternalMove)

        # grid.addWidget(self.parameters['sample_library_files'], 0, 1)
        self.parameters['sample_library_files'].setMinimumHeight(30)
        self.parameters['sample_library_files'].setMaximumHeight(50)
        self.parameters['sample_library_files'].setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)

        # 创建按钮布局
        sample_buttons_layout = QVBoxLayout()
        add_sample_btn = QPushButton("Add")
        remove_sample_btn = QPushButton("Remove")
        add_sample_btn.clicked.connect(self.add_sample_library_files)
        remove_sample_btn.clicked.connect(self.remove_sample_library_files)
        sample_buttons_layout.addWidget(add_sample_btn)
        sample_buttons_layout.addWidget(remove_sample_btn)
        
        grid.addWidget(self.parameters['sample_library_files'], 0, 1)
        grid.addLayout(sample_buttons_layout, 0, 2)
        
        # SysteMHC Library section
        grid.addWidget(QLabel("SysteMHC Library:"), 1, 0)
        self.parameters['systemhc_library_files'] = QListWidget()
        # 启用多选功能
        self.parameters['systemhc_library_files'].setSelectionMode(QListWidget.ExtendedSelection)
        # 连接键盘删除键
        self.parameters['systemhc_library_files'].keyPressEvent = lambda event: self.handle_key_press(
            event, self.parameters['systemhc_library_files'], self.remove_systemhc_library_files)
        # 右键菜单
        self.parameters['systemhc_library_files'].setContextMenuPolicy(Qt.CustomContextMenu)
        self.parameters['systemhc_library_files'].customContextMenuRequested.connect(
            lambda pos: self.show_context_menu(pos, self.parameters['systemhc_library_files'], self.remove_systemhc_library_files))
        # 启用拖放重排
        self.parameters['systemhc_library_files'].setDragDropMode(QListWidget.InternalMove)

        # grid.addWidget(self.parameters['systemhc_library_files'], 1, 1)
        self.parameters['systemhc_library_files'].setMinimumHeight(80)
        self.parameters['systemhc_library_files'].setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        # 创建按钮布局
        systemhc_buttons_layout = QVBoxLayout()
        remove_systemhc_btn = QPushButton("Remove")
        remove_systemhc_btn.clicked.connect(self.remove_systemhc_library_files)
        systemhc_buttons_layout.addWidget(remove_systemhc_btn)
        grid.addWidget(self.parameters['systemhc_library_files'], 1, 1)
        grid.addLayout(systemhc_buttons_layout, 1, 2)
        
        param_group.setLayout(grid)
        self.main_layout.addWidget(param_group)

    def add_sample_library_files(self):
        """Add sample library files based on selected pipeline"""
        # if self.selected_pipeline == "FragPipe":
        #     file_filter = "TSV Files (*.tsv)"
        #     valid_ext = ".tsv"
        # elif self.selected_pipeline == "SysteMHC":
        #     file_filter = "SPTXT Files (*.sptxt)"
        #     valid_ext = ".sptxt"
        # else:  # MaxQuant
        #     file_filter = "TXT Files (msms.txt)"
        #     valid_ext = ".txt"
        
        files, _ = QFileDialog.getOpenFileNames(self, "Add Sample Library Files", "", self.library_filter)
        
        # if not files:
        #     return
            
        # # Validate file extensions
        # invalid_files = []
        # for f in files:
        #     if self.selected_pipeline == "MaxQuant":
        #         # 对于MaxQuant，检查文件名是否为"msms.txt"
        #         if os.path.basename(f) != "msms.txt":
        #             invalid_files.append(f)
        #     elif not f.lower().endswith(valid_ext):
        #         invalid_files.append(f)
        
        # if invalid_files:
        #     QMessageBox.warning(
        #         self, 
        #         "Invalid Files", 
        #         f"The following files have invalid formats for {self.selected_pipeline} pipeline:\n" + 
        #         "\n".join(invalid_files)
        #     )
        #     files = [f for f in files if f not in invalid_files]
            
        if files:
            self.parameters['sample_library_files'].addItems(files)


    def handle_key_press(self, event, list_widget, remove_function):
        """Handle key press events for list widgets"""
        # 处理删除键和退格键
        if event.key() in (Qt.Key_Delete, Qt.Key_Backspace) and list_widget.selectedItems():
            remove_function()
        else:
            # 调用原始的keyPressEvent方法处理其他按键
            QListWidget.keyPressEvent(list_widget, event)

    def show_context_menu(self, pos, list_widget, remove_function):
        """Show context menu for list widgets"""
        # 只有在有选中项时才显示菜单
        if not list_widget.selectedItems():
            return
            
        context_menu = QMenu()
        remove_action = context_menu.addAction("Remove")
        
        # 显示菜单并获取选择的操作
        action = context_menu.exec_(list_widget.mapToGlobal(pos))
        
        if action == remove_action:
            remove_function()

    def remove_input_files(self):
        """Remove selected input files from the list"""
        selected_items = self.input_files_list.selectedItems()
        if not selected_items:
            QMessageBox.information(self, "Information", "Please select file(s) to remove")
            return
        
        # 从列表中移除选中的项
        for item in selected_items:
            self.input_files_list.takeItem(self.input_files_list.row(item))

    def remove_sample_library_files(self):
        """Remove selected sample library files from the list"""
        selected_items = self.parameters['sample_library_files'].selectedItems()
        if not selected_items:
            QMessageBox.information(self, "Information", "Please select file(s) to remove")
            return
        
        # 从列表中移除选中的项
        for item in selected_items:
            self.parameters['sample_library_files'].takeItem(
                self.parameters['sample_library_files'].row(item))
        
        self.output_area.append(f"Removed {len(selected_items)} sample library file(s)")

    def remove_systemhc_library_files(self):
        """Remove selected SysteMHC library files from the list"""
        selected_items = self.parameters['systemhc_library_files'].selectedItems()
        if not selected_items:
            QMessageBox.information(self, "Information", "Please select file(s) to remove")
            return
        
        # 从列表中移除选中的项
        for item in selected_items:
            self.parameters['systemhc_library_files'].takeItem(
                self.parameters['systemhc_library_files'].row(item))
        
        if self.main_logger:
            self.main_logger.write(f"Removed {len(selected_items)} SysteMHC library file(s)\n")
        else:
            self.output_area.append(f"Removed {len(selected_items)} SysteMHC library file(s)")

    def create_extra_params_section(self):
        """Create additional parameters section (three columns)"""
        self.extra_params_widget = QWidget()
        main_h_layout = QHBoxLayout()
        main_h_layout.setContentsMargins(0, 0, 0, 0)  # 每行也不留边距
        main_h_layout.setSpacing(0)
        # Create three vertical layout containers
        left_column = QVBoxLayout()
        middle_column = QVBoxLayout()
        right_column = QVBoxLayout()
        
        # Generate parameter input fields
        self.extra_param_inputs = {}
        params = list(self.default_params.items())
        
        # Calculate items per column (distribute evenly across 3 columns)
        items_per_column = (len(params) + 2) // 3

        LABEL_WIDTH = 90
        INPUT_WIDTH = 100

        for i, (param, default_value) in enumerate(params):
            param_layout = QHBoxLayout()
            param_layout.setContentsMargins(0, 0, 0, 0)  # 每行也不留边距
            param_layout.setSpacing(0)

            param_label = QLabel(f"{param}")
            param_label.setFixedWidth(LABEL_WIDTH)
            param_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)

            BOOLEAN_PARAMS = {
                'reanalyse',
                'xic-theoretical-fr',
            }

            if param in BOOLEAN_PARAMS:
                combo = QComboBox()
                combo.addItems(["false", "true"])  # 注意顺序：默认是 "false"
                combo.setFixedWidth(INPUT_WIDTH) 
                # combo.setCurrentText(default_value.lower())  # 确保匹配
                # 保存引用以便后续读取值
                self.parameters['reanalyse_combo'] = combo
                param_input = combo
                self.extra_param_inputs[param] = param_input

            else:
                param_input = QLineEdit()
                param_input.setText(default_value)
                param_input.setFixedWidth(INPUT_WIDTH)
                param_input.setAlignment(Qt.AlignLeft)  # ← 数值右对齐

                self.extra_param_inputs[param] = param_input
                
            param_layout.addWidget(param_label)
            param_layout.addWidget(param_input)
            param_layout.addStretch()
            # Add to appropriate column
            column_index = i // items_per_column
            if column_index == 0:
                left_column.addLayout(param_layout)
            elif column_index == 1:
                middle_column.addLayout(param_layout)
            else:
                right_column.addLayout(param_layout)
        
        # Add columns to main layout
        container = QWidget()
        container_layout = QHBoxLayout()
        
        container_layout.addLayout(left_column)
        container_layout.addLayout(middle_column)
        container_layout.addLayout(right_column)
        container.setLayout(container_layout)
        
        main_h_layout.addWidget(container)
        self.extra_params_widget.setLayout(main_h_layout)
        self.extra_params_widget.setVisible(True)  
        self.main_layout.addWidget(self.extra_params_widget)

    def toggle_extra_params_section(self):
        """Toggle visibility of additional parameters section"""
        self.extra_params_widget.setVisible(not self.extra_params_widget.isVisible())
    
    def create_allele_specific_layout(self):
        """Create layout for allele-specific library input with dual-control UX"""

        allele_specific_layout = QVBoxLayout()

        # label
        top_row_layout = QHBoxLayout()
        allele_label = QLabel("Input Alleles   ")
        top_row_layout.addWidget(allele_label)

        #class 
        self.allele_class_combo = QComboBox()  # ← 现在创建了！
        self.allele_class_combo.addItems(["ClassI", "ClassII"])
        top_row_layout.addWidget(QLabel("Class:"))
        top_row_layout.addWidget(self.allele_class_combo)

        # right： choose alleles
        self.allele_selector = QComboBox()
        top_row_layout.addWidget(QLabel("Allele:"))
        self.allele_selector.setEditable(True)
        self.allele_selector.setInsertPolicy(QComboBox.NoInsert)
        self.allele_selector.lineEdit().setPlaceholderText("Allele to search...")

        # completer
        self.allele_completer = QCompleter()
        self.allele_completer.setCompletionMode(QCompleter.UnfilteredPopupCompletion)
        self.allele_completer.setCaseSensitivity(Qt.CaseInsensitive)

        # ✅ 关键：启用子串匹配（SubstringMatch）
        self.allele_completer.setFilterMode(Qt.MatchContains)  # ← 这是核心！
        self.allele_selector.lineEdit().setCompleter(self.allele_completer)
        # self.allele_selector.setCompleter(self.allele_completer)
        self.update_allele_completer()
        top_row_layout.addWidget(self.allele_selector)


        # 可选：设置 stretch 让 selector 占更多空间
        top_row_layout.setStretchFactor(allele_label, 0)
        top_row_layout.setStretchFactor(self.allele_class_combo, 0)
        top_row_layout.setStretchFactor(self.allele_selector, 1)
        
        allele_specific_layout.addLayout(top_row_layout)


        bottom_row_layout = QHBoxLayout()
        # 左侧：显示已选等位基因（只读，美观）
        self.allele_display = QLineEdit()
        self.allele_display.setPlaceholderText("e.g., HLA-A01_01;HLA-B07_02")
        self.allele_display.setReadOnly(False)  # 可设为 False 允许手动编辑
        self.allele_display.setStyleSheet("QLineEdit { background-color: #f0f0f0; }")
        bottom_row_layout.addWidget(self.allele_display, stretch=3)

        # Download 按钮
        btn_download_allele = QPushButton("Download Allele Library")
        btn_download_allele.clicked.connect(self.download_and_add_allele_library)
        bottom_row_layout.addWidget(btn_download_allele, stretch=1)

        allele_specific_layout.addLayout(bottom_row_layout)
        
        self.main_layout.addLayout(allele_specific_layout)

        # 现在可以安全调用 update_allele_completer()
        self.update_allele_completer()

        # 连接信号（必须在创建后）
        self.allele_class_combo.currentTextChanged.connect(self.update_allele_completer)
        self.allele_class_combo.currentTextChanged.connect(self.on_allele_class_changed)

        # 连接选择信号
        self.allele_selector.lineEdit().returnPressed.connect(self.add_selected_allele)
        self.allele_selector.activated.connect(self.add_selected_allele)

    def on_allele_class_changed(self, text):
        # 可用于日志、重置等，目前可为空
        pass

    def add_selected_allele(self):
        """将右侧选择的等位基因添加到左侧显示框，并清空选择器"""
        text = self.allele_selector.currentText().strip()
        if not text:
            return

        # # 验证格式（可选）
        # import re
        # # if not re.match(r"^HLA-[A-Z][\d_]+$", text):
        # if not re.match(r"^[HLA|DR|DP|DQ]+$", text):
        #     QMessageBox.warning(self, "Invalid Format", f"Ignored invalid allele: {text}")
        #     self.allele_selector.clearEditText()
        #     return

        # 获取当前已选列表
        current = self.allele_display.text().strip()
        if current:
            new_text = current + "; " + text
        else:
            new_text = text

        self.allele_display.setText(new_text)
        self.allele_selector.clearEditText()  # 清空输入，准备下一次选择
        self.allele_selector.setFocus()       # 保持焦点在输入框

    def update_allele_completer(self):
        """Update autocomplete based on selected class"""
        selected_class = self.allele_class_combo.currentText()
        alleles = self.allele_list.get(selected_class, [])

        # self.allele_selector.addItems(alleles)

        model = QStringListModel(alleles)
        self.allele_completer.setModel(model)

        self.allele_selector.clear()

    def filter_completer(self, text):
        """Dynamically filter completer content based on input"""
        # 获取光标位置前的文本，拆分出当前正在输入的段
        line_edit = self.allele_selector.lineEdit()
        cursor_pos = line_edit.cursorPosition()
        selected_class = self.allele_class_combo.currentText()
        all_alleles = self.allele_list.get(selected_class, []) 
        
        # 分割多等位基因输入（以 ; 为分隔）
        parts = []
        start = 0
        for i, ch in enumerate(text):
            if ch == ';':
                parts.append((start, i))
                start = i + 1
        parts.append((start, len(text)))

        # 找到光标所在的 segment
        current_segment = ""
        for start, end in parts:
            if start <= cursor_pos <= end:
                current_segment = text[start:cursor_pos].strip()
                break

        # 过滤：只有非空输入才匹配
        if current_segment:
            filtered = [
                a for a in all_alleles
                if current_segment.lower() in a.lower()
            ]
        else:
            filtered = []  # 空输入 → 无下拉

        # 更新模型（空模型 = 不弹出）
        model = QStringListModel(filtered)
        self.allele_completer.setModel(model)

        # semicolon_positions = [-1]  # 起始位置视为 -1（虚拟分号）
        # for i, char in enumerate(text):
        #     if char == ';':
        #         semicolon_positions.append(i)

        # # 找到光标落在哪个 segment（即最后一个 <= cursor_pos 的分号索引）
        # last_semi_idx = -1
        # for pos in reversed(semicolon_positions):
        #     if pos < cursor_pos:
        #         last_semi_idx = pos
        #         break

        # # 当前段 = 从 last_semi_idx+1 到 cursor_pos 的子串
        # current_segment = text[last_semi_idx + 1:cursor_pos].strip()

        # # 过滤匹配项（不区分大小写）
        # filtered_alleles = [
        #     allele for allele in all_alleles
        #     if current_segment.lower() in allele.lower()
        # ]

        # # 更新模型
        # model = QStringListModel(filtered_alleles)
        # self.allele_completer.setModel(model)

        # if filtered_alleles and current_segment:
        #     self.allele_completer.setCompletionPrefix(current_segment)
        #     # 必须用 QTimer 才能可靠弹出补全框
        #     from PyQt5.QtCore import QTimer
        #     QTimer.singleShot(0, self.allele_completer.complete)
        # else:
        #     self.allele_completer.popup().hide()

    def on_completion_selected(self, selected_text):
        """补全选中后自动添加分号"""
        # 获取当前输入框文本和光标位置
        current_text = self.allele_specific_input.text()
        cursor_pos = self.allele_specific_input.cursorPosition()

        # 找到最后一个在光标前的分号位置
        last_semi = -1
        for i in range(cursor_pos - 1, -1, -1):
            if current_text[i] == ';':
                last_semi = i
                break

        # 构建新文本：保留光标前到 last_semi+1 的部分 + 选中项 + "; " + 光标后内容
        prefix = current_text[:last_semi + 1]
        if prefix and not prefix.endswith('; '):
            # 确保格式统一为 "HLA-XXX_XX; "
            if prefix.strip() and not prefix.rstrip().endswith(';'):
                prefix = prefix.rstrip() + "; "
            elif not prefix.strip():
                prefix = ""
        suffix = current_text[cursor_pos:].lstrip()

        new_text = prefix + selected_text + "; " + suffix
        self.allele_specific_input.setText(new_text)

        # 光标定位到新增项之后
        new_cursor_pos = len(prefix) + len(selected_text) + 2  # +2 for "; "
        self.allele_specific_input.setCursorPosition(new_cursor_pos)

    def eventFilter(self, obj, event):
        if obj == self.allele_specific_input and event.type() == QEvent.KeyPress:
            key = event.key()
            model = self.allele_completer.model()
            
            if key == Qt.Key_Tab:
                if model and model.rowCount() > 0:
                    first = model.index(0, 0).data()
                    self.on_completion_selected(first)
                    return True
                return False  # 允许 Tab 切焦点
            
            elif key in (Qt.Key_Enter, Qt.Key_Return):
                if model and model.rowCount() > 0:
                    first = model.index(0, 0).data()
                    self.on_completion_selected(first)
                    return True
                else:
                    # 阻止换行，但不清空输入
                    return True
                    
        return super().eventFilter(obj, event)

    def parse_allele_input(self):
        """从 allele_display 中解析等位基因列表"""
        input_text = self.allele_display.text().strip()
        if not input_text:
            return []
        
        alleles = [p.strip() for p in input_text.split(";") if p.strip()]
        # 去重（保持顺序）
        seen = set()
        unique = []
        for a in alleles:
            if a not in seen:
                unique.append(a)
                seen.add(a)
        
        # # 格式校验
        valid = []
        # # hla_pattern = re.compile(r"^[HLA]-[A-Z][\d_]+$")
        # hla_pattern = re.compile(r"^[HLA|DR|DP|DQ]+$")
        for allele in unique:
            valid.append(allele)
            # if hla_pattern.match(allele):
            #     valid.append(allele)
            # else:
            #     QMessageBox.warning(self, "Warning", f"Invalid allele ignored: {allele}")
        return valid
    
    def download_next_allele(self):

        attrs = self.batch_download_attrs
        
        if attrs['current_allele_idx'] >= len(attrs['alleles']):
            self.on_batch_download_finished()
            return
        
        # 下载当前等位基因
        allele = attrs['alleles'][attrs['current_allele_idx']]
        attrs['progress'].setLabelText(f"Downloading: {allele} ({attrs['current_allele_idx']+1}/{len(attrs['alleles'])})")
        
        # 启动单文件下载线程
        thread = SingleDownloadThread(allele, attrs['selected_class'], attrs['output_dir'])
        thread.download_done.connect(self.on_single_download_done)
        thread.start()
        attrs['download_threads'].append(thread)  # 保存线程引用
        attrs['current_allele_idx'] += 1

    def on_single_download_done(self, success, msg, allele):
        attrs = self.batch_download_attrs

        if success:
            attrs['success_count'] += 1
            dest_path = os.path.join(attrs['output_dir'], f"HCD_cons_{allele}_top12_bynam_ptm.tsv")
            self.parameters['systemhc_library_files'].addItem(dest_path)
            if self.main_logger:
                self.main_logger.write(f"{msg}\n")
            else:
                self.output_area.append(msg)
            
        else:
            attrs['failed_alleles'].append(allele)
            if self.main_logger:
                self.main_logger.write(f"Error: {msg}\n")
            else:
                self.output_area.append(f"Error: {msg}")
        
        # 更新进度
        attrs['progress'].setValue(attrs['current_allele_idx'])
        
        # # 继续下载下一个（直接调用类方法，无需self）
        # if not attrs.get('cancel_flag', False):
        #     QTimer.singleShot(100, self.download_next_allele)
        QTimer.singleShot(100, self.download_next_allele)

    def on_batch_download_finished(self):
        attrs = self.batch_download_attrs
        if not attrs:
            return
        attrs['progress'].close()  # 直接关闭
        QMessageBox.information(self, "Done", "Download completed.")

        # progress = attrs.get('progress')
        # if progress:  # and not progress.wasCanceled():
        #     progress.setCancelButtonText("OK")
        #     progress.setLabelText("Download completed!")

        # progress.setWindowFlags(progress.windowFlags() | Qt.WindowCloseButtonHint)
        # progress.show()

    def download_and_add_allele_library(self):
        """Download allele-specific library file directly with loading indicator"""

        # if hasattr(self, 'batch_download_attrs') and self.batch_download_attrs.get('progress'):
        #     self.batch_download_attrs['progress'].close()
        #     self.batch_download_attrs['cancel_flag'] = True
        #     for thread in self.batch_download_attrs['download_threads']:
        #         thread.terminate()
        #     self.batch_download_attrs['download_threads'] = []

        alleles = self.parse_allele_input()
        if not alleles:
            QMessageBox.warning(self, "Warning", "Please enter at least one allele name!")
            return

        output_dir = self.output_dir_field.text().strip()
        if not output_dir:
            QMessageBox.warning(self, "Warning", "Please select the output directory first!")
            return
            
        # Make sure output directory exists
        if not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to create output directory: {str(e)}")
                return
                
        # Get selected class
        selected_class = self.allele_class_combo.currentText()  # ClassI or ClassII

        progress = QProgressDialog(
            f"Downloading {len(alleles)} alleles...", 
            "", 
            0, 
            len(alleles), 
            self
        )
        progress.setWindowTitle("Downloading Libraries ...")
        progress.setWindowModality(Qt.WindowModal)
        progress.setMinimumDuration(0)
        progress.setCancelButton(None)
        # progress.setAutoClose(False)
        # progress.setAutoReset(False)

        progress.show()
        QApplication.processEvents()  # ✅ 强制刷新界面，显示进度条

        self.batch_download_attrs = {
            'alleles': alleles,
            'selected_class': selected_class,
            'output_dir': output_dir,
            'success_count': 0,
            'failed_alleles': [],
            'current_allele_idx': 0,
            'download_threads': [],
            'progress': progress  # ← 保存引用
        }
        
        self.download_next_allele()

    def merge_libraries(self):
        """Merge sample and SysteMHC libraries by directly calling the appropriate scripts"""
        output_dir = self.output_dir_field.text().strip()
        if not output_dir:
            raise Exception("Output directory not specified")
            
        # 创建日志文件路径
        log_file_path = os.path.join(output_dir, 'lib-base-result.log.txt')
        
        # 获取库文件
        sample_libs = [self.parameters['sample_library_files'].item(i).text() 
                    for i in range(self.parameters['sample_library_files'].count())]
        systemhc_libs = [self.parameters['systemhc_library_files'].item(i).text()
                        for i in range(self.parameters['systemhc_library_files'].count())]
                        
        if not sample_libs or not systemhc_libs:
            raise Exception("Both sample and SysteMHC libraries must be provided")
        
        # 选择第一个样本库
        sample_library_path = sample_libs[0]
        
        # # 创建日志重定向器
        # logger = TeeLogger(self.output_area, log_file_path)
        
        try:
            # 使用contextlib重定向标准输出和错误
            with redirect_stdout(self.main_logger), redirect_stderr(self.main_logger):
                try:
                    # 根据所选管道导入适当的模块
                    if self.selected_pipeline == "FragPipe":  
                        try:
                            self.main_logger.write(f"\n===== Merging libraries from SysteMHC and libraries from FragPipe =====\n")
                            merged_library = fragpipe_api.merge_libraries(
                                sample_library_path=sample_library_path,
                                systemhc_lib_paths=systemhc_libs,
                                output_dir=output_dir
                            )
                            self.main_logger.write(f"\n===== Merging libraries completed =====\n")
                            return merged_library
                        except ImportError:                           
                            self.main_logger.write(f"\n===== Merging libraries Failed =====\n")
                            return merged_library
                    elif self.selected_pipeline == "SysteMHC":
                        try:
                            self.main_logger.write(f"\n===== Merging libraries from SysteMHC and libraries from SysteMHC-pipeline =====\n")
                            merged_library = systemhc_api.merge_libraries(
                                sample_library_path=sample_library_path,
                                systemhc_lib_paths=systemhc_libs,
                                output_dir=output_dir
                            )
                            self.main_logger.write(f"\n===== Merging libraries completed =====\n")
                            return merged_library
                        except ImportError:
                            self.main_logger.write(f"\n===== Merging libraries Failed =====\n")
                            return merged_library
                    else:  # MaxQuant
                        try:
                            self.main_logger.write(f"\n===== Merging libraries from SysteMHC and libraries from MaxQuant =====\n")
                            merged_library = maxquant_api.merge_libraries(
                                sample_library_path=sample_library_path,
                                systemhc_lib_paths=systemhc_libs,
                                output_dir=output_dir
                            )
                            self.main_logger.write(f"\n===== Merging libraries completed =====\n")
                            return merged_library
                        except ImportError:
                            self.main_logger.write(f"\n===== Merging libraries Failed =====\n")
                            return merged_library
                except Exception as e:
                    # 捕获并记录在API调用中发生的任何错误
                    error_msg = f"Error during library merging: {str(e)}"
                    self.main_logger.write(error_msg + "\n")
                    # 重新抛出异常，让调用代码知道发生了错误
                    raise Exception(error_msg)
        finally:
            # 确保关闭日志器
            # self.main_logger.close()
            pass

    def run_rescore_step(self):
        """Execute the Rescore module after main pipeline success."""
        input_dir = self.input_dir_field.text().strip()
        output_dir = self.output_dir_field.text().strip()

        # 推断必要文件路径
        first_pass = os.path.join(output_dir, "lib-base-result-first-pass.parquet")
        main_result = os.path.join(output_dir, "lib-base-result.parquet")
        if os.path.exists(first_pass):
            report_file = first_pass          # 👉 如果存在，用 first-pass
        elif os.path.exists(main_result):
            report_file = main_result         # 👉 只有 first-pass 不存在时，才退而求其次用主结果
        # 4. 都没有就报错
        else:
            raise FileNotFoundError("No valid report file found for Rescore!")

        speclib_file = None
        for f in os.listdir(output_dir):
            if f.startswith("merged_Sample+SysteMHC_library") and f.endswith(".tsv"):
                speclib_file = os.path.join(output_dir, f)
                break

        if not os.path.exists(report_file):
            raise FileNotFoundError(f"Required report file not found: {report_file}")
        if not speclib_file or not os.path.exists(speclib_file):
            raise FileNotFoundError("Merged library file not found in output directory.")

        model = getattr(self, 'rescore_model', 'DNN')

        if not hasattr(self, 'main_logger') or self.main_logger is None:
            raise RuntimeError("Main logger not initialized!")

        try:
        # 重定向所有 print() 和 sys.stdout 到 GUI + 文件
            
            # 1. 保存原始的 sys.stdout/sys.stderr
            original_stdout = sys.stdout
            original_stderr = sys.stderr
            # 2. 重定向 sys.stdout/sys.stderr 到 main_logger
            sys.stdout = self.main_logger
            sys.stderr = self.main_logger

            # === 关键：临时重定向 logging 到 main_logger ===
            log_handler = QTextEditHandler(self.main_logger)
            log_formatter = logging.Formatter('%(asctime)s> %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
            log_handler.setFormatter(log_formatter)

            # 获取根 logger 并添加 handler
            root_logger = logging.getLogger()
            original_handlers = root_logger.handlers.copy()
            original_level = root_logger.level
            root_logger.handlers.clear()
            root_logger.setLevel(logging.INFO)
            root_logger.addHandler(log_handler)

            self.main_logger.write(f"▶ Performing Rescore with {model}...\n")
            QApplication.processEvents()

            # 动态导入 rescore 模块（避免启动时依赖）
            from dia_aspire.config import IOConfig, FineTuneConfig
            from dia_aspire.config.pipeline import Pipeline

            io_config = IOConfig(
                ms_file_dir=input_dir,
                report_file=report_file,
                ms_file_type="mzML",  # 可扩展支持 raw/d
                output_dir=output_dir,
            )

            finetune_config = FineTuneConfig(
                fdr_threshold=0.01,
                instrument='QE',
                nce=27,
                psm_num_to_train_ms2=8000,
                epoch_to_train_ms2=20,
                epoch_to_train_rt_ccs=25,
                train_verbose=True,
            )

            pipeline = Pipeline(
                io_config=io_config,
                finetune_config=finetune_config,
                feature_generators=["basic", "ms2", "rt", "xic"],
                model=model,
            )

            # 执行 rescore
            psm_df = pipeline.run()
            self.main_logger.write(f"\n✅ Rescore completed successfully. \n")
    
    # try:
    #     from dia_aspire_rescore.plot import plot_qvalues
    #     fig_path = os.path.join(output_dir, f"rescore_qvalue_{model.lower()}.png")
    #     plot_qvalues(psm_df, save_path=fig_path)
    #     self.output_area.append(f"📊 Q-value plot saved: {fig_path}")
    # except Exception as plot_e:
    #     self.output_area.append(f"⚠️ Plotting failed: {plot_e}")
            return None
        
        except Exception as e:
            error_msg = f"❌ Rescore failed: {str(e)}"
            self.main_logger.write(error_msg)
            import traceback
            self.main_logger.write(traceback.format_exc() + "\n")
            raise

        finally:
            # ========== 关键：恢复原始配置 ==========
            # 1. 恢复 sys.stdout/sys.stderr
            sys.stdout = original_stdout
            sys.stderr = original_stderr
            
            # 2. 恢复 logging 根 logger
            root_logger = logging.getLogger()
            root_logger.removeHandler(log_handler)
            root_logger.handlers = original_handlers
            root_logger.setLevel(original_level)
            
            # 3. 刷新 GUI
            QApplication.processEvents()


    def execute_command(self):
        if self.process.state() == QProcess.Running:
            QMessageBox.warning(self, 'Warning', 'A task is already running!')
            return
        
        self.rescore_enabled = self.rescore_checkbox.isChecked()
        if self.rescore_enabled:
            self.rescore_model = "DNN" if self.rescore_dnn_radio.isChecked() else "SVM"

        try:
            # Get input path or files
            input_type = self.input_type_combo.currentText()
            input_dir = self.input_dir_field.text().strip()
            input_files = [self.input_files_list.item(i).text() for i in range(self.input_files_list.count())]

            # Get output directory and DIANN path
            output_dir = self.output_dir_field.text().strip()
            diann_path = self.diann_path_field.text().strip()

            # Basic validation
            if input_type == "Folder" and not input_dir:
                QMessageBox.critical(self, 'Error', 'Please select an input directory!')
                return
            if input_type == "Files" and not input_files:
                QMessageBox.critical(self, 'Error', 'Please add input files!')
                return
            if not output_dir:
                QMessageBox.critical(self, 'Error', 'Please select an output directory!')
                return
            if not diann_path:
                QMessageBox.critical(self, 'Error', 'Please specify the DIANN software path!')
                return

            # Check if both library types have files
            sample_lib_count = self.parameters['sample_library_files'].count()
            systemhc_lib_count = self.parameters['systemhc_library_files'].count()
            
            if sample_lib_count == 0 or systemhc_lib_count == 0:
                QMessageBox.critical(self, 'Error', 'Both Sample Library and SysteMHC Library must contain files!')
                return

            # Check if paths exist
            if input_type == "Folder" and not os.path.exists(input_dir):
                QMessageBox.critical(self, 'Error', f'The input directory does not exist:\n{input_dir}')
                return
            if not os.path.exists(output_dir):
                try:
                    os.makedirs(output_dir)
                except Exception as e:
                    QMessageBox.critical(self, 'Error', f'Failed to create output directory:\n{str(e)}')
                    return

            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            self.gui_log_path = os.path.join(output_dir, f"DIA-Aspire_GUI_{timestamp}.log.txt")

            self.main_logger = TeeLogger(self.output_area, self.gui_log_path)

            self.main_logger.write(f"\n▶ Workflow Starting ...\n")
            self.output_area.verticalScrollBar().setValue(
                self.output_area.verticalScrollBar().maximum()
            )
            QApplication.processEvents()  # ⭐ 关键：强制刷新 GUI


            # 日志文件路径
            log_file_path = os.path.join(output_dir, 'lib-base-result.log.txt')
            
            # Copy irt_SYSTEMHC.csv to output directory if it exists
            # source_irt_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "constants", "irt_SYSTEMHC.csv")
            dest_irt_file = os.path.join(output_dir, "irt_SYSTEMHC.parquet")

            try:
                with importlib_resources.files('dia_aspire.constants').joinpath('irt_SYSTEMHC.parquet') as source_irt_file:
                    if source_irt_file.exists():
                        import shutil
                        shutil.copy2(source_irt_file, dest_irt_file)
                        # self.main_logger.write(f"Copied IRT file to: {dest_irt_file}\n")
            except Exception as e:
                self.main_logger.write(f"Warning: Failed to copy IRT file: {e}\n")

            # 记录使用的库文件（替代“Downloaded”日志）
            for i in range(systemhc_lib_count):
                lib_path = self.parameters['systemhc_library_files'].item(i).text()
                self.main_logger.write(f"Using SysteMHC library: {os.path.basename(lib_path)}\n")
            for i in range(sample_lib_count):
                lib_path = self.parameters['sample_library_files'].item(i).text()
                self.main_logger.write(f"Using Sample library: {os.path.basename(lib_path)}\n")


            # First merge libraries - 这部分已经在merge_libraries中处理了日志记录
            # self.main_logger.write("\nMerging Sample-specific and SysteMHC libraries ...\n")
            QApplication.processEvents()  # 可选：每步都刷新

            with redirect_stdout(self.main_logger), redirect_stderr(self.main_logger):
                merged_library = self.merge_libraries()
            self.main_logger.write(f"Libraries merged successfully: {merged_library} \n")


            # Build command with merged library
            with resources.files('dia_aspire.diann').joinpath('rundiann_file.sh') as script_path:
                script_path = str(script_path)
                os.chmod(script_path, 0o755)
                command = [script_path]

                # Add input path or files
                if input_type == "Folder":
                    command.extend(["--dir", input_dir])
                elif input_type == "Files":
                    for file in input_files:
                        command.extend(["--f", file])

                # Add merged library and output parameters
                command.extend(["--lib", merged_library])
                command.extend(["--out", "./lib-base-result.parquet"])
                command.extend(["--output-dir", output_dir])
                command.extend(["--log-file", log_file_path])  # 传递日志文件路径
                
                # Add DIANN path
                command.extend(["--diann-path", diann_path])

                # Add ALL additional parameters, including those with default values
                for param, input_widget in self.extra_param_inputs.items():
                    if isinstance(input_widget, QLineEdit):
                        value = input_widget.text().strip()
                    elif isinstance(input_widget, QComboBox):
                        value = input_widget.currentText().strip()
                    else:
                        # 其他控件类型（如 QCheckBox 等）可继续扩展
                        continue

                    # Always add the parameter, even if using default value
                    if value != 'false':
                        if value == 'true':
                            command.append(f"--{param}")
                        else:
                            command.append(f"--{param}")
                            command.append(value)
                # Display command
                display_cmd = ' '.join(command)
                # print(display_cmd)
                self.main_logger.write(f"[run command] {display_cmd}\n")
                QApplication.processEvents()

                # # 初始化 main_logger（追加模式）
                # self.main_logger = TeeLogger(self.output_area, self.gui_log_path)  # 'a' 模式
                self.process.setProcessChannelMode(QProcess.MergedChannels)
                # Start process
                self.process.setReadChannel(QProcess.StandardOutput)  # 确保优先读取标准输出
                self.process.start("bash", command)
                self.output_area.append("\n▶ Starting DIANN...\n")
                
        except Exception as e:
            error_msg = f"Failed to start workflow: {str(e)}"
            if hasattr(self, 'main_logger') and self.main_logger:
                self.main_logger.write(f"\n❌ {error_msg}\n")
                # self.main_logger.close()
            else:
                self.output_area.append(f'<span style="color: red;">{error_msg}</span>')
            QMessageBox.critical(self, 'Error', error_msg)

    def save_gui_log(self, output_dir=None):
        """保存GUI界面的日志到文件"""
        if not output_dir:
            output_dir = self.output_dir_field.text().strip()
            if not output_dir:
                QMessageBox.warning(self, "Warning", "Please specify an output directory first!")
                return False
        
        # 确保输出目录存在
        try:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to create output directory: {str(e)}")
            return False
        
        # 创建带有时间戳的日志文件名
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_filename = f"DIA-Aspire_GUI_log_{timestamp}.txt"
        log_filepath = os.path.join(output_dir, log_filename)
        
        # 保存日志
        try:
            with open(log_filepath, 'w', encoding='utf-8') as f:
                # 获取纯文本内容
                text_content = self.output_area.toPlainText()
                f.write(text_content)
            if self.main_logger:
                self.main_logger.write(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] GUI log saved to: {log_filepath}\n")
            else:
                self.output_area.append(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] GUI log saved to: {log_filepath}")
            return True
            
        except Exception as e:
            error_msg = f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Error saving GUI log: {str(e)}"
            if self.main_logger:
                self.main_logger.write(error_msg + "\n")
            else:
                self.output_area.append(error_msg)
            return False
        
    def stop_task(self):
        if self.process.state() == QProcess.Running:
            self.process.kill()
            self.output_area.append("\nThe task has been manually terminated")

    # def handle_output(self):
    #     data = self.process.readAllStandardOutput().data().decode()
    #     if data.strip():
            
    #         # 如果有日志文件，也写入日志
    #         if hasattr(self, 'main_logger'):
    #             self.main_logger.write(data)
    #         else:
    #             self.output_area.append(data.strip())

    # def handle_error(self):
    #     error = self.process.readAllStandardError().data().decode()
    #     if error.strip():
            
    #         # 如果有日志文件，也写入日志（纯文本格式）
    #         if hasattr(self, 'main_logger'):
    #             self.main_logger.write(f"ERROR: {error.strip()}\n")
    #         else:
    #             self.output_area.append(f'<span style="color: red;">{error.strip()}</span>')

    def handle_output(self):
        data = self.process.readAllStandardOutput().data()
        if data:
            try:
                text = data.decode('utf-8')
            except UnicodeDecodeError:
                text = data.decode('utf-8', errors='replace')

            if text.strip():
                if hasattr(self, 'main_logger'):
                    self.main_logger.write(text)
                    self.main_logger.flush()
                    QApplication.processEvents()
                else:
                    self.output_area.append(text.rstrip())

    def handle_error(self):
        data = self.process.readAllStandardError().data()
        if data:
            try:
                text = data.decode('utf-8')
            except UnicodeDecodeError:
                text = data.decode('utf-8', errors='replace')

            if text.strip():
                if hasattr(self, 'main_logger'):
                    self.main_logger.write(f"ERROR: {text}")
                    self.main_logger.flush()
                    QApplication.processEvents()
                else:
                    self.output_area.append(f'<span style="color: red;">{text.rstrip()}</span>')

    def task_finished(self, exit_code):

        # ========== 关键修复：强制读取所有剩余输出 ==========
        # 读取所有未处理的标准输出
        remaining_output = self.process.readAllStandardOutput().data()
        if remaining_output:
            try:
                text = remaining_output.decode('utf-8')
            except UnicodeDecodeError:
                text = remaining_output.decode('utf-8', errors='replace')
            if text.strip() and hasattr(self, 'main_logger'):
                self.main_logger.write(text)
        
        # 读取所有未处理的错误输出
        remaining_error = self.process.readAllStandardError().data()
        if remaining_error:
            try:
                text = remaining_error.decode('utf-8')
            except UnicodeDecodeError:
                text = remaining_error.decode('utf-8', errors='replace')
            if text.strip() and hasattr(self, 'main_logger'):
                self.main_logger.write(f"ERROR: {text}")
        
        # 强制刷新 GUI 事件循环
        QApplication.processEvents()


        if exit_code == 0:
            self.main_logger.write("\n✅ Library search completed successfully.\n")
            QApplication.processEvents()

            # 如果启用了 Rescore，则运行它
            if getattr(self, 'rescore_enabled', False):
                self.main_logger.write("\n▶ Starting Rescore module...\n")
                QApplication.processEvents()
                try:
                    self.run_rescore_step()
                    self.main_logger.write("\n✅ Entire workflow completed!\n")
                    QMessageBox.information(self, 'Success', 'Entire workflow completed!')
                except Exception as e:
                    import traceback
                    error_msg = f"Rescore failed: {str(e)}\n"
                    error_msg += traceback.format_exc() + "\n"
                    self.main_logger.write(error_msg)
                    self.output_area.append(error_msg)
                    QMessageBox.critical(self, 'Error', f"Rescore failed: {str(e)}")
            else:
                QMessageBox.information(self, 'Success', 'Library search completed!')

        else:
            error_msg = f'Task failed, exit code: {exit_code}\n'
            self.main_logger.write(error_msg)
            QMessageBox.critical(self, 'Error', error_msg)

        # 关闭日志文件
        if hasattr(self, 'main_logger'):
            self.main_logger.close()
            self.main_logger = None
        
        # 保存GUI日志
        # self.save_gui_log()

    def toggle_input_type(self):
        """Toggle between folder and file input modes"""
        input_type = self.input_type_combo.currentText()
        if input_type == "Folder":
            self.input_dir_widget.setVisible(True)
            self.input_files_widget.setVisible(False)
        elif input_type == "Files":
            self.input_dir_widget.setVisible(False)
            self.input_files_widget.setVisible(True)

    def select_input(self):
        """Select input directory"""
        dir_path = QFileDialog.getExistingDirectory(self, "Select Input Directory")
        if dir_path:
            self.input_dir_field.setText(dir_path)

    def add_input_files(self):
        """Add multiple input files"""
        files, _ = QFileDialog.getOpenFileNames(self, "Add Input Files", "", "Data Files (*.mzML *.raw)")
        if files:
            self.input_files_list.addItems(files)

    def select_output_dir(self):
        """Select output directory"""
        dir_path = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dir_path:
            self.output_dir_field.setText(dir_path)


def main():
    import sys
    from PyQt5.QtWidgets import QApplication

    app = QApplication(sys.argv)
    window = CommandLineGUI()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()