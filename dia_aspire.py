import sys
import os
import subprocess
import tarfile
import requests
import json
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QDesktopWidget,
                             QLabel, QLineEdit, QPushButton, QComboBox, QGroupBox, QGridLayout,
                             QListWidget, QMessageBox, QTextEdit, QFileDialog, QCheckBox, QRadioButton,
                             QMenu,QCompleter)
from PyQt5.QtCore import Qt, QStringListModel, QProcess


# 将src目录添加到导入路径
src_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "src")
sys.path.append(src_dir)

# 现在可以从src导入模块
try:
    from src import fragpipe_api
    from src import systemhc_api
    from src import sptxt2tsv
    from src import irt_alignment
except ImportError:
    # 如果添加路径后仍然无法导入，尝试直接导入
    import fragpipe_api
    import systemhc_api
    import sptxt2tsv
    import irt_alignment

class CommandLineGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.parameters = {
                    'sample_library_files': None,
                    'systemhc_library_files': None,  # New list for SysteMHC libraries
                    'extra_params': {}
                }
        self.default_params = {
            'threads': '32',
            'verbose':'5',
            'qvalue': '0.1',
            'matrix-qvalue': '0.01',
            'matrices':'true',
            'mass-acc': '15',
            'mass-acc-ms1': '15',
            'double-search': 'true',
            'reanalyse':'true',
            'no-prot-inf': 'true',
            'rt-profiling': 'true',
            'pg-level': '1',
            'report-lib-info': 'true'
        }
        self.allele_list = {}
        self.process = QProcess(self)
        self.main_layout = QVBoxLayout()
        self.extra_params_widget = None
        self.selected_pipeline = "FragPipe"  # Default pipeline
        self.initUI()
        self.load_allele_list()

        # Connect signals
        self.process.readyReadStandardOutput.connect(self.handle_output)
        self.process.readyReadStandardError.connect(self.handle_error)
        self.process.finished.connect(self.task_finished)

    def load_allele_list(self):
        """Load allele list from JSON file"""
        try:
            # 更新JSON文件路径
            json_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "src", "allele_list.json")
            with open(json_path, "r") as f:
                self.allele_list = json.load(f)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load allele list:\n{str(e)}")
            self.allele_list = {"ClassI": [], "ClassII": []}

    def initUI(self):
        self.setWindowTitle('DIA-Aspire')
        self.setGeometry(300, 300, 800, 600)
        self.center_window()
        
        # 更新logo路径
        logo_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "src", "DIA-Aspire_logo.jpg")
        if os.path.exists(logo_path):
            from PyQt5.QtGui import QIcon, QPixmap
            self.setWindowIcon(QIcon(logo_path))
            
            # 在应用程序顶部显示logo
            logo_layout = QHBoxLayout()
            logo_label = QLabel()
            pixmap = QPixmap(logo_path)
            # 将logo调整为合适的大小
            scaled_pixmap = pixmap.scaled(200, 100, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            logo_label.setPixmap(scaled_pixmap)
            logo_label.setAlignment(Qt.AlignCenter)
            
            logo_layout.addWidget(logo_label)
            self.main_layout.addLayout(logo_layout)

        # Create input/output paths
        self.create_io_layout()
        
        # Create pipeline selector
        self.create_pipeline_selector()
        
        # Create library input sections
        self.create_parameter_inputs()
        
        # Create Allele Specific section
        self.create_allele_specific_layout()

        # Create additional parameters section
        toggle_extra_params_btn = QPushButton("Parameters")
        toggle_extra_params_btn.clicked.connect(self.toggle_extra_params_section)
        self.main_layout.addWidget(toggle_extra_params_btn)
        self.create_extra_params_section()

        # Action buttons
        button_layout = QHBoxLayout()
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
            background-color: #1E1E1E;
            color: #D4D4D4;
            font-family: Consolas;
            font-size: 12pt;
            border: 2px solid #3E3E3E;
            padding: 5px;
        }
        """)

        self.main_layout.addLayout(button_layout)
        self.main_layout.addWidget(self.output_area)
        self.setLayout(self.main_layout)

    def center_window(self):
        screen_geometry = QDesktopWidget().screenGeometry()
        x = (screen_geometry.width() - self.width()) // 2
        y = (screen_geometry.height() - self.height()) // 2
        self.move(x, y)

    def create_io_layout(self):
        """Create input and output path layout"""
        # Input type selection
        input_type_layout = QHBoxLayout()
        input_type_layout.addWidget(QLabel("DIA Input (raw,mzML,diaPASEF.d):"))
        self.input_type_combo = QComboBox()
        self.input_type_combo.addItems(["Folder", "Files"])
        self.input_type_combo.currentTextChanged.connect(self.toggle_input_type)
        input_type_layout.addWidget(self.input_type_combo)
        self.main_layout.addLayout(input_type_layout)

        # Input directory selection (folder mode)
        self.input_dir_widget = QWidget()
        input_dir_layout = QHBoxLayout(self.input_dir_widget)
        input_dir_layout.addWidget(QLabel("Input Directory:"))
        self.input_dir_field = QLineEdit()
        btn_browse_input = QPushButton("Browse")
        btn_browse_input.clicked.connect(self.select_input)
        input_dir_layout.addWidget(self.input_dir_field)
        input_dir_layout.addWidget(btn_browse_input)

        # Input files selection (files mode)
        self.input_files_widget = QWidget()
        input_files_layout = QVBoxLayout(self.input_files_widget)
        input_files_layout.addWidget(QLabel("Input Files:"))
        self.input_files_list = QListWidget()
        input_files_layout.addWidget(self.input_files_list)
        btn_add_file = QPushButton("Add Files")
        btn_add_file.clicked.connect(self.add_input_files)
        input_files_layout.addWidget(btn_add_file)

        # Default to folder mode
        self.input_files_widget.setVisible(False)
        self.main_layout.addWidget(self.input_dir_widget)
        self.main_layout.addWidget(self.input_files_widget)

        # Output directory
        output_dir_layout = QHBoxLayout()
        output_dir_layout.addWidget(QLabel("Output Directory:"))
        self.output_dir_field = QLineEdit()
        btn_browse_output = QPushButton("Browse")
        btn_browse_output.clicked.connect(self.select_output_dir)
        output_dir_layout.addWidget(self.output_dir_field)
        output_dir_layout.addWidget(btn_browse_output)
        self.main_layout.addLayout(output_dir_layout)
        
        # DIANN Software Path
        diann_path_layout = QHBoxLayout()
        diann_path_layout.addWidget(QLabel("DIANN Software Path:"))
        self.diann_path_field = QLineEdit()
        self.diann_path_field.setText("/usr/diann/1.8.1/diann-1.8.1")  # Default path
        btn_browse_diann = QPushButton("Browse")
        btn_browse_diann.clicked.connect(self.select_diann_path)
        diann_path_layout.addWidget(self.diann_path_field)
        diann_path_layout.addWidget(btn_browse_diann)
        self.main_layout.addLayout(diann_path_layout)
        
    def select_diann_path(self):
        """Select DIANN software path"""
        file_path, _ = QFileDialog.getOpenFileName(self, "Select DIANN Executable")
        if file_path:
            self.diann_path_field.setText(file_path)
        
    def create_pipeline_selector(self):
        """Create pipeline selector (FragPipe or SysteMHC)"""
        pipeline_group = QGroupBox("Pipeline used in sample library generation")
        pipeline_layout = QHBoxLayout()
        
        self.fragpipe_radio = QRadioButton("FragPipe")
        self.systemhc_radio = QRadioButton("SysteMHC-pipeline")
        self.fragpipe_radio.setChecked(True)  # Default to FragPipe
        
        self.fragpipe_radio.toggled.connect(self.on_pipeline_changed)
        self.systemhc_radio.toggled.connect(self.on_pipeline_changed)
        
        pipeline_layout.addWidget(self.fragpipe_radio)
        pipeline_layout.addWidget(self.systemhc_radio)
        pipeline_group.setLayout(pipeline_layout)
        
        self.main_layout.addWidget(pipeline_group)
    
    def on_pipeline_changed(self):
        """Handle pipeline selection change"""
        if self.fragpipe_radio.isChecked():
            self.selected_pipeline = "FragPipe"
            self.sample_lib_label.setText("Sample Library (tsv files):")
        else:
            self.selected_pipeline = "SysteMHC"
            self.sample_lib_label.setText("Sample Library (sptxt files):")
            
    def create_parameter_inputs(self):
        """Create library input sections with remove functionality"""
        param_group = QGroupBox("Library Input")
        grid = QGridLayout()
        
        # Sample Library section
        self.sample_lib_label = QLabel("Sample Library (tsv files):")
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
        
        # 创建按钮布局
        systemhc_buttons_layout = QVBoxLayout()
        remove_systemhc_btn = QPushButton("Remove")
        remove_systemhc_btn.clicked.connect(self.remove_systemhc_library_files)
        systemhc_buttons_layout.addWidget(remove_systemhc_btn)
        
        grid.addWidget(self.parameters['systemhc_library_files'], 1, 1)
        grid.addLayout(systemhc_buttons_layout, 1, 2)
        
        param_group.setLayout(grid)
        self.main_layout.addWidget(param_group)

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
        
        self.output_area.append(f"Removed {len(selected_items)} SysteMHC library file(s)")

    def create_extra_params_section(self):
        """Create additional parameters section (two columns)"""
        self.extra_params_widget = QWidget()
        main_h_layout = QHBoxLayout()
        
        # Create two vertical layout containers
        left_column = QVBoxLayout()
        right_column = QVBoxLayout()
        
        # Generate parameter input fields
        self.extra_param_inputs = {}
        params = list(self.default_params.items())
        half_len = (len(params) + 1) // 2
        
        for i, (param, default_value) in enumerate(params):
            param_layout = QHBoxLayout()
            param_label = QLabel(f"{param}:")
            param_input = QLineEdit()
            param_input.setText(default_value)
            self.extra_param_inputs[param] = param_input
            param_layout.addWidget(param_label)
            param_layout.addWidget(param_input)
            
            # Add to appropriate column
            if i < half_len:
                left_column.addLayout(param_layout)
            else:
                right_column.addLayout(param_layout)
        
        # Add columns to main layout
        container = QWidget()
        container_layout = QHBoxLayout()
        container_layout.addLayout(left_column)
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
        """Create Allele Specific section"""
        allele_specific_layout = QHBoxLayout()

        # Section label
        allele_specific_label = QLabel("SysteMHC Allele-Specific Library:")
        allele_specific_layout.addWidget(allele_specific_label)

        # Class selection
        self.allele_class_combo = QComboBox()
        self.allele_class_combo.addItems(["ClassI", "ClassII"])
        self.allele_class_combo.currentTextChanged.connect(self.update_allele_completer)
        allele_specific_layout.addWidget(self.allele_class_combo)

        # Allele input with autocomplete
        self.allele_specific_input = QLineEdit()
        self.allele_specific_input.setPlaceholderText("Input allele name (HLA-A01_01)")
        allele_specific_layout.addWidget(self.allele_specific_input)

        # Autocomplete setup
        self.allele_completer = QCompleter()
        self.allele_completer.setCaseSensitivity(Qt.CaseInsensitive)
        self.allele_completer.setFilterMode(Qt.MatchContains)
        self.allele_specific_input.setCompleter(self.allele_completer)

        # Download button
        btn_download_allele = QPushButton("Download")
        btn_download_allele.clicked.connect(self.download_and_add_allele_library)
        allele_specific_layout.addWidget(btn_download_allele)

        self.main_layout.addLayout(allele_specific_layout)

        # Initialize completer
        self.update_allele_completer()

    def update_allele_completer(self):
        """Update autocomplete based on selected class"""
        selected_class = self.allele_class_combo.currentText()
        alleles = self.allele_list.get(selected_class, [])
        model = QStringListModel(alleles)
        self.allele_completer.setModel(model)

        # Connect text changed event for dynamic filtering
        self.allele_specific_input.textChanged.connect(self.filter_completer)

    def filter_completer(self, text):
        """Dynamically filter completer content based on input"""
        selected_class = self.allele_class_combo.currentText()
        all_alleles = self.allele_list.get(selected_class, [])
        filtered_alleles = [allele for allele in all_alleles if text.lower() in allele.lower()]
        model = QStringListModel(filtered_alleles)
        self.allele_completer.setModel(model)

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

    def add_sample_library_files(self):
        """Add sample library files based on selected pipeline"""
        file_filter = "TSV Files (*.tsv)" if self.selected_pipeline == "FragPipe" else "SPTXT Files (*.sptxt)"
        files, _ = QFileDialog.getOpenFileNames(self, "Add Sample Library Files", "", file_filter)
        
        if not files:
            return
            
        # Validate file extensions
        valid_ext = ".tsv" if self.selected_pipeline == "FragPipe" else ".sptxt"
        invalid_files = [f for f in files if not f.lower().endswith(valid_ext)]
        
        if invalid_files:
            QMessageBox.warning(
                self, 
                "Invalid Files", 
                f"The following files have invalid extensions for {self.selected_pipeline} pipeline:\n" + 
                "\n".join(invalid_files)
            )
            files = [f for f in files if f not in invalid_files]
            
        if files:
            self.parameters['sample_library_files'].addItems(files)

    def download_and_add_allele_library(self):
        """Download allele-specific library file directly"""
        allele_name = self.allele_specific_input.text().strip()
        if not allele_name:
            QMessageBox.warning(self, "Warning", "Please enter the allele name!")
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
        
        # Construct the URL based on allele name and class
        file_name = f"HCD_cons_{allele_name}_top12_bynam_ptm.tsv"
        url = f"https://systemhc.sjtu.edu.cn/data/Systemhc_v2_2023/Data/230623_build/SysteMHC_Library/{selected_class}/Allele-specific/{file_name}"
        
        # Destination path for the downloaded file
        dest_path = os.path.join(output_dir, file_name)
        
        try:
            # Download the file
            self.output_area.append(f"Downloading {file_name}...")
            response = requests.get(url, stream=True)
            
            if response.status_code != 200:
                raise Exception(f"Download Error, HTTP status code: {response.status_code}")
            
            # Save the file
            with open(dest_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            
            # Add to SysteMHC library list
            self.parameters['systemhc_library_files'].addItem(dest_path)
            
            self.output_area.append(f"Downloaded library file: {file_name}")
            QMessageBox.information(self, "Success", f"Added library file: {file_name}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Download error: {str(e)}")
            self.output_area.append(f"Error: {str(e)}")

    def merge_libraries(self):
        """Merge sample and SysteMHC libraries by directly calling the appropriate scripts"""
        output_dir = self.output_dir_field.text().strip()
        if not output_dir:
            raise Exception("Output directory not specified")
            
        # Get library files
        sample_libs = [self.parameters['sample_library_files'].item(i).text() 
                    for i in range(self.parameters['sample_library_files'].count())]
        systemhc_libs = [self.parameters['systemhc_library_files'].item(i).text()
                        for i in range(self.parameters['systemhc_library_files'].count())]
                        
        if not sample_libs or not systemhc_libs:
            raise Exception("Both sample and SysteMHC libraries must be provided")
        
        # Choose the first sample library
        sample_library_path = sample_libs[0]
        
        # Import the appropriate module based on pipeline selection
        if self.selected_pipeline == "FragPipe":
            # 从src导入
            try:
                from src import fragpipe_api
                self.output_area.append(f"Merging libraries using FragPipe pipeline...")
                merged_library = fragpipe_api.merge_libraries(
                    sample_library_path=sample_library_path,
                    systemhc_lib_paths=systemhc_libs,
                    output_dir=output_dir
                )
                return merged_library
            except ImportError:
                # 如果无法从src导入，尝试直接导入
                import fragpipe_api
                self.output_area.append(f"Merging libraries using FragPipe pipeline...")
                merged_library = fragpipe_api.merge_libraries(
                    sample_library_path=sample_library_path,
                    systemhc_lib_paths=systemhc_libs,
                    output_dir=output_dir
                )
                return merged_library
        else:
            # 从src导入
            try:
                from src import systemhc_api
                self.output_area.append(f"Merging libraries using SysteMHC pipeline...")
                merged_library = systemhc_api.merge_libraries(
                    sample_library_path=sample_library_path,
                    systemhc_lib_paths=systemhc_libs,
                    output_dir=output_dir
                )
                return merged_library
            except ImportError:
                # 如果无法从src导入，使用importlib
                import importlib.util
                spec = importlib.util.spec_from_file_location(
                    "systemhc_api", 
                    os.path.join(os.path.dirname(os.path.abspath(__file__)), "src", "systemhc_api.py")
                )
                systemhc_api = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(systemhc_api)
                
                self.output_area.append(f"Merging libraries using SysteMHC pipeline...")
                merged_library = systemhc_api.merge_libraries(
                    sample_library_path=sample_library_path,
                    systemhc_lib_paths=systemhc_libs,
                    output_dir=output_dir
                )
                return merged_library

    def execute_command(self):
        if self.process.state() == QProcess.Running:
            QMessageBox.warning(self, 'Warning', 'A task is already running!')
            return

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
                    
            # Copy irt_SYSTEMHC.csv to output directory if it exists
            source_irt_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "src", "irt_SYSTEMHC.csv")
            dest_irt_file = os.path.join(output_dir, "irt_SYSTEMHC.csv")
            
            try:
                if os.path.exists(source_irt_file):
                    import shutil
                    shutil.copy2(source_irt_file, dest_irt_file)
                    self.output_area.append(f"Copied irt_SYSTEMHC.csv to output directory")
                else:
                    self.output_area.append("Warning: irt_SYSTEMHC.csv not found in application directory")
            except Exception as e:
                self.output_area.append(f"Warning: Failed to copy irt_SYSTEMHC.csv: {str(e)}")
                        
            # First merge libraries
            self.output_area.append("Merging Sample and SysteMHC libraries...")
            merged_library = self.merge_libraries()
            self.output_area.append(f"Libraries merged successfully: {merged_library}")

            # Build command with merged library
            script_path = os.path.abspath(os.path.join("src", "rundiann_file.sh"))
            command = [script_path]

            # Add input path or files
            if input_type == "Folder":
                command.extend(["--dir", input_dir])
            elif input_type == "Files":
                for file in input_files:
                    command.extend(["--f", file])

            # Add merged library and output parameters
            command.extend(["--lib", merged_library])
            command.extend(["--out", "lib-base-result"])
            command.extend(["--output-dir", output_dir])
            
            # Add DIANN path
            command.extend(["--diann-path", diann_path])

            # Add ALL additional parameters, including those with default values
            for param, input_widget in self.extra_param_inputs.items():
                value = input_widget.text().strip()
                # Always add the parameter, even if using default value
                command.append(f"--{param}")
                command.append(value)

            # Display command
            display_cmd = ' '.join(command)
            self.output_area.append(f"[run command] {display_cmd}\n")

            # Start process
            self.process.start("bash", command)
            self.output_area.append("▶ Start processing data...\n")
            
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Execution failed: {str(e)}')

    def stop_task(self):
        if self.process.state() == QProcess.Running:
            self.process.kill()
            self.output_area.append("\nThe task has been manually terminated")

    def handle_output(self):
        data = self.process.readAllStandardOutput().data().decode()
        self.output_area.append(data.strip())

    def handle_error(self):
        error = self.process.readAllStandardError().data().decode()
        self.output_area.append(f'<span style="color: red;">{error.strip()}</span>')

    def task_finished(self, exit_code):
        if exit_code == 0:
            QMessageBox.information(self, 'Success', 'Mission accomplished!')
        else:
            QMessageBox.critical(self, 'Error', f'Task failed, exit code: {exit_code}')


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = CommandLineGUI()
    window.show()
    sys.exit(app.exec_())