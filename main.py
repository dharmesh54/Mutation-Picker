import sys
import os
import pandas as pd
import numpy
import numpy as np
from pandas.api.types import is_string_dtype
import warnings
import csv
import time
from PyQt5.QtGui import QPalette, QPixmap, QBrush
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QFileDialog, QLabel, QLineEdit, QCheckBox
from PyQt5.QtGui import QMovie
from PyQt5.QtCore import Qt
from PyQt5.QtCore import Qt, QThread, pyqtSignal
import warnings
warnings.filterwarnings("ignore")
from task import run_script

class Worker(QThread):
    progressChanged = pyqtSignal(int)
    
    def __init__(self, file_name, path, wild_id, total_gen):
        super().__init__()
        self.file_name = file_name
        self.path = path
        self.wild_id = wild_id
        self.total_gen = total_gen
    
    def stop(self):
        self.is_running = False

    def run(self):
        self.is_running = True
        self.run_script()
        self.is_running = False
         
    def run_script(self):
        print("Starting script execution")  
    
        file_name = self.file_name
        path = self.path
        Wild_id = self.wild_id
        total_gen = self.total_gen
        
        print("File Name:", file_name)
        print("Path:", path)
        print("Wild ID:", Wild_id)
        print("Total Gen:", total_gen)
       

        # Run the script in the worker thread
        print("Running script...")
        run_script(file_name, path, Wild_id, total_gen)

        # Re-enable the run button when the script finishes
        #print("Enabling run button")
        

        # Perform any necessary actions with the script output
        # ...

        # Emit a signal to indicate that the script has finished
        print("Script execution finished")
        self.progressChanged.emit(-1)
        
        #self.worker_finished = True
        #self.worker.stop()

class App(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.worker = self.worker = Worker("", "", "", 0)
        self.worker_finished = False  # Add the worker_finished attribute
        self.worker.progressChanged.connect(self.handleProgressChanged)
        self.worker_progress = -1
        


    def initUI(self):
	    
        palette = QPalette()
        background_image = QPixmap("codon_chart.jpg")
        palette.setBrush(QPalette.Window, QBrush(background_image))
        self.setPalette(palette)
        	
	
        self.file_label = QLabel('File Name:', self)
        self.file_label.move(50, 20)

        self.file_input = QLineEdit(self)
        self.file_input.move(120, 20)
        self.file_input.resize(150, 20)

        self.path_label = QLabel('Path:', self)
        self.path_label.move(50, 50)

        self.path_input = QLineEdit(self)
        self.path_input.move(120, 50)
        self.path_input.resize(150, 20)

        self.wild_id_label = QLabel('Wild ID:', self)
        self.wild_id_label.move(50, 80)

        self.wild_id_input = QLineEdit(self)
        self.wild_id_input.move(120, 80)
        self.wild_id_input.resize(150, 20)

        self.gen_label = QLabel('Gen:', self)
        self.gen_label.move(50, 110)

        self.gen_input = QLineEdit(self)
        self.gen_input.move(120, 110)
        self.gen_input.resize(150, 20)

        self.open_file_button = QPushButton('Open File', self)
        self.open_file_button.setToolTip('Select a file to process')
        self.open_file_button.move(50, 140)
        self.open_file_button.clicked.connect(self.open_file_dialog)

        self.check_box = QCheckBox('Enable Feature', self)
        self.check_box.move(50, 170)
        self.check_box.stateChanged.connect(self.on_check_box_changed)

        self.quit_button = QPushButton('Quit', self)
        self.quit_button.setToolTip('Exit application')
        self.quit_button.move(150, 210)
        self.quit_button.clicked.connect(self.close)

        self.run_button = QPushButton('Run', self)
        self.run_button.setToolTip('Run script')
        self.run_button.move(50, 210)
        self.run_button.clicked.connect(self.run_task)

        self.loading_label = QLabel(self)
        self.loading_label.setFixedSize(40, 40)
        self.loading_label.move(220, 210)
        self.loading_spinner = QMovie('pc7reEGKi.gif')
        self.loading_label.setMovie(self.loading_spinner)

        self.setGeometry(300, 300, 300, 250)
        self.setWindowTitle('Mutation Picker')
        self.show()
        
        # Prevent window from being maximized
        self.setFixedSize(self.size())


    def open_file_dialog(self):
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getOpenFileName(self,"Open File", "","All Files (*);;Python Files (*.py)", options=options)
        if file_path:
            self.file_input.setText(file_path)
            file_name = file_path.split("/")[-1]
            path = file_path[:len(file_path)-len(file_name)]
            self.path_input.setText(path)
        else:
            file_name = self.file_input.text()
            path = self.path_input.text()
            if file_name and path:
                file_path = path + '/' + file_name

    def on_check_box_changed(self, state):
        if state ==  Qt.Checked:
             #add script functionality when the checkbox is checked
             print('Feature Enabled')
        else:
             #add script functionality when the checkbox is unchecked
             print('Feature Disabled')
    
    
    def handleProgressChanged(self, value):
        if value == -1:
            # Script finished, update GUI accordingly
            self.loading_label.clear()
            self.run_button.setEnabled(True)  # Enable the run button
            self.worker.stop()  # Stop the worker thread
        else:
            # Update progress bar or timer
            self.loading_spinner.start()
            # Update other GUI elements as needed  
            
    def run_task(self):
        print("Processing")
          
        # Get input values from GUI
        file_name = self.file_input.text()
        if "/" in file_name:
            file_name = file_name.split("/")[-1]
        else:
            file_name =  File_name
        path = self.path_input.text()
        Wild_id = self.wild_id_input.text()
        total_gen = int(self.gen_input.text())
        
        #print("Creating worker instance...")
        # Create a Worker instance and pass the input values
        self.worker = Worker(file_name, path, Wild_id, total_gen)
        #self.worker.worker_progress.connect(self.handleProgressChanged)
        self.worker.progressChanged.connect(self.handleProgressChanged)
        
        # Disable the run button while the script is running
        #print("Disabling run button")
        self.run_button.setEnabled(False)
        
        self.loading_spinner.start()
        
        #print("Starting worker thread...")
        # Start the worker thread
        self.worker.start()
        #print("Worker thread finished")
        #self.worker.stop()
        
        #print("Enabling run button")
        #self.run_button.setEnabled(True)
       
        # Add a loop to wait for the worker thread to finish
        #while not self.worker_finished:
        #QApplication.processEvents()  # Process events to avoid freezing GUI
        
        #print("Enabling run button")
        #self.run_button.setEnabled(True)
        

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    ex.show()
    sys.exit(app.exec_())
