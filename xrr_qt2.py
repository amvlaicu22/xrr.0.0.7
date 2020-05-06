import sys
import time
import random
from threading import Thread
import numpy as np

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtCore import QThread
from PyQt5 import QtWidgets, uic

# ~ from queue import Queue
# ~ from PyQt5.QtCore import Qt, QObject, pyqtSlot
# ~ from PyQt5.QtGui import QIcon, QPalette
# ~ from PyQt5.QtWidgets import (QApplication, QMainWindow,  
# ~ QVBoxLayout, QHBoxLayout, QGridLayout, QFormLayout, QFormLayout, 
# ~ QGroupBox, QSizePolicy, QMessageBox, QDialogButtonBox, QLineEdit, 
# ~ QMenu, QMenuBar, QSpinBox, QTextEdit, QFrame, QWidget, QPushButton, 
# ~ QLabel, QComboBox, QDialog, QCheckBox, QRadioButton, QPushButton, 
# ~ QTableWidget, QSlider, QProgressBar)

# ~ from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
# ~ from matplotlib.backends.backend_qt4 import NavigationToolbar2QT as NavigationToolbar
# ~ from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

from XRR import XRR, _XLayer
# ~ from where import *
# ~ from functools import wraps

class Panel_Plot(QWidget):
    def __init__(self, parent, id=-1, dpi=None, figsize=(3,2), 
        tbar = 0, row=1, col=1, fplots=None):     #, **kwargs):
        # ~ def __init__(self, parent):
        # ~ super(QWidget, self).__init__(parent)
        QWidget.__init__(self, parent)
        # ~ wx.Panel.__init__(self, parent, id=id)  
        
        self.parent = parent
        self.row = row
        self.col = col
        self.plt = fplots
        # ~ 
        self.fig = Figure(dpi=dpi, figsize=figsize)
        # ~ self.cvs = mpl.figure.FigureCanvas(self, -1, self.fig)
        self.cvs = FigureCanvas(self.fig)
        self.cvs.setSizePolicy( QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.cvs.updateGeometry()
        # ~ for plt pyplot
        # ~ self.fig, self.axs = plt.subplots(row, col)
        # ~ for mpl.Figure
        self.axs = self.fig.subplots(row, col)   #, sharex=True)
        self.fig.subplots_adjust(top=1.0,
            bottom=0.1,
            left=0.12,
            right=0.8,
            hspace=0.3,
            wspace=0.2)
        # ~ # Adjust the scaling factor to fit your legend text completely outside the plot
        # ~ # (smaller value results in more space being made for the legend)

        self.lay = QVBoxLayout(self)
        self.setLayout(self.lay)
        self.lay.addWidget(self.cvs)
        if tbar:
            self.nav = NavigationToolbar(self.cvs, self)
            # ~ self.nav.Realize()
            self.lay.addWidget(self.nav)
            
        self.sizer_ = QVBoxLayout(parent)
        self.parent.setLayout(self.sizer_)
        self.sizer_.addWidget(self)
        
    def plot(self):
        axs = np.asarray(self.axs)
        # ~ print(axs.flat)
        # ~ input()
        for fplot, axs in zip(self.plt, axs.flat): fplot(axs)
        self.draw()
    def draw(self):
        # ~ self.fig.tight_layout()
        self.cvs.draw()
        # ~ self.fig.canvas.draw()
        # ~ self.fig.canvas.flush_events()

class SocketQ(QObject):
    '''
    use :sys.stdout = SocketQ(signal=self.OnPrint)
    '''
    signal = pyqtSignal(object)
    def write(self, obj):self.signal.emit(obj)
    def flush(self): pass


# GUI Frame class that spins off the worker thread
class MainFrame(QMainWindow, XRR):
    signal_run = pyqtSignal(object)
    signal_end = pyqtSignal(object)
    signal_upd = pyqtSignal(object)
    
    class _Socket(QObject):
        # wrapper to work like wx.PyEvent
        
        class _data():
            def __init__(self, data): 
                self.data = data
                
        signal = pyqtSignal(object)
        
        def __init__(self, dest=None, func=None):
            QObject.__init__(self)
            if dest and func:
                self.connect(dest, func)
        def connect(self, dest, func):
            self.dest = dest
            self.func = func
            self.signal.connect(func)
            return self
        def emit(self, data): 
            self.write(data)
        def write(self, data): 
            self.data = data
            self.signal.emit(self._data(data))
        def flush(self): pass
    

    def __init__(self):
        super(MainFrame, self).__init__()
        uic.loadUi('xrr_qt2gui.ui', self)
        self.show()
        # ~ self.B_Fit_Start = self.findChild(QtWidgets.QPushButton, 'B_Fit_Start')
        self.B_Fit_Start.   released.connect(self.OnButton_B_Fit_Start) 
        self.B_Fit_Stop.    released.connect(self.OnButton_B_Fit_Stop) 
        self.B_Data_Load.   released.connect(self.OnButton_B_Data_Load) 
        self.B_Layer_Load.  released.connect(self.OnButton_B_Layer_Load) 
        self.B_Layer_Save.  released.connect(self.OnButton_B_Layer_Save) 
        self.B_Layer_Insert.released.connect(self.OnButton_B_Layer_Insert) 
        self.B_Layer_Delete.released.connect(self.OnButton_B_Layer_Delete) 
        self.B_Layer_Update.released.connect(self.OnButton_B_Layer_Update) 
        
        self.Menu()   
        self.statusBar().showMessage('Message in statusbar.')
        
        self._plot()
        self.layer_to_grid()
        # ~ self.print_2socket()
        sys.stdout = SocketQ(signal=self.print_out)
        sys.stderr = SocketQ(signal=self.print_err)
        
        self.print_out("test out\n")
        self.print_err("test err\n")
        print('test')
        
            
    def closeEvent(self, event):
        # ~ Override the closeEvent method of QWidget in your main window.
        print( "closeEvent()")
        if self.DoYouWantTo(" SAVE results "):
            event.ignore()
            return
        if not self.DoYouWantTo( "continue FITTING"):
            print("closing")
            event.accept() # let the window close
        else:
            event.ignore()
    ######################################
    
    def print_color(self,r,g,b):
        self.T_Logger.moveCursor(QTextCursor.End)
        self.T_Logger.setTextColor( QColor(r,g,b) )
    
    def print_(self, text):
        font = QFont('Mono', 6, QFont.Light)
        self.T_Logger.setFont(font)
        self.T_Logger.insertPlainText( text )
        self.T_Logger.ensureCursorVisible()

    def print_out(self, text):
        self.print_(text)

    def print_err(self, text):
        self.print_color(255,0,0)
        self.print_(text)
        self.print_color(0,0,0)
        
    def _plot_init(self):
        # ~ self.plot_anchx = 1.02
        # ~ self.plot_anchy = 1.1
        # ~ self.plot_drag = True

        if not self.plot_plots:
            p1 = Panel_Plot(self.Panel1_Data, row=2, col=1, tbar=1,
                fplots=[self._plot_refl1, self._plot_refl2])
            p2 = Panel_Plot(self.Panel2_Profile, row=2, col=1, tbar=1,
                fplots=[self._plot_prof1, self._plot_prof2])
            p3 = Panel_Plot(self.Panel3_Error, row=1, col=1, 
                fplots=[self._plot_err])
            self.plot_plots = [p1, p2, p3]
        # ~ # ~~~~
        # ~ for p in self.plot_plots: p.plot()

       
    def fit_prestart(self, event):  # wxGlade: XRR_Frame.<event_handler>
        self.B_Fit_Start.setEnabled(False)
        self.B_Fit_N.setText('0')
        self.T_Logger.clear()

        self.layer_from_grid()
        self.layer_init()
        self.layer_to_grid()
        
        fitkeys = []
        if self.Chk_L_ab.isChecked():            fitkeys.append('ab')
        if self.Chk_L_dd.isChecked():            fitkeys.append('dd')
        if self.Chk_L_rh.isChecked():            fitkeys.append('rh')
        if self.Chk_L_sg.isChecked():            fitkeys.append('sg')
        
        self.fit_set_fitkeys(fitkeys)
        
        if self.RB_Fit_LM.isChecked() : self.fitmode = 0
        if self.RB_Fit_DE.isChecked() : self.fitmode = 1
        
        if self.RB_Fit_Layer.isChecked():    self.fitprof = 0
        if self.RB_Fit_Profile.isChecked():  self.fitprof = 1
        
        self.update = self.V_Fit_Update.value()
        self._pdz   = self.V_Profile_Step.value()
        
        fitkeysP = []
        if self.Chk_P_ab.isChecked():        fitkeysP.append('ab')
        if self.Chk_P_rh.isChecked():        fitkeysP.append('rh')
        if self.Chk_P_ac.isChecked():        fitkeysP.append('at')
        # ~ at = self.combo_box_atom.isChecked()
        # ~ atom = self.Chk_P_atom.GetValue()
        
               
    def OnButton_B_Fit_Start(self, event=None): 
        self.fit_start(event)

    def OnButton_B_Fit_Stop(self, event=None):  
        self.fit_stop(event)
        if self.fit_worker:
            if not self.fit_worker.is_alive():
                self.B_Fit_Start.setDisabled(False)
        else:   # if self.fit_worker == None:
            self.B_Fit_Start.setDisabled(False)
                

    
    def fit_run(self, event):
        XRR.fit_run(self, event)
        self.B_Fit_N.setText(str(self.fitn))
        
    def fit_upd(self, event):
        # ~ XRR.fit_upd(self, event)
        print( f'\nxrr>>: update {event.data}', flush=True)
        if not self.fitprof :   # layer
            self.layer_profile()
        self.layer_to_grid()
        self._plot()

        
    def fit_end(self, event):
        print( f'\nxrr>>: finish {event.data}', flush=True)
        self.fit_worker = None
        self.B_Fit_Start.setEnabled(True)
        self.layer_to_grid()
        self._plot()
        
 
    
    def layer_to_grid(self):
        data = []
        data.append(_XLayer._hd)
        for L in self.LL:
            row = [L._name, L._comp, "{:.1f}".format(L._dd), "{:.1f}".format(L._sg),
            "{:.4f}".format(L._rh), "{:.1f}".format(L._Mm), "{:.1f}".format(L._Vm) ]
            data.append(row)
        self.setData_rows(self.grid_1, data)
        
    def setData_rows(self, grid, 
        data_rows=[[   'name',     'comp.',        'd [A]',    's [A]',    'g/cm^3',   'g/mol',    'A^3'],
              [   'glass',    'Si;1;O;2',    '0',         '5.0',      '2.6',      '',         ''],
              [   'W',        'W;1',          '10.0',     '5.0',      '20.0',     '',         '']]): 
        # ~ self.data_row = data
        nr = len(data_rows) -1 # 1st row is for header
        nc = len(data_rows[0])
        grid.setRowCount(nr)
        grid.setColumnCount(nc)
        grid.setHorizontalHeaderLabels(data_rows[0])
        for r, row in enumerate(data_rows[1:], start=0):
            for c, item in enumerate(row):
                newitem = QTableWidgetItem(item)
                grid.setItem(r, c, newitem)
        grid.resizeColumnsToContents()
        grid.resizeRowsToContents()

    
    def layer_from_grid(self):
        for r, L in enumerate(self.LL):
            L._name  = self.grid_1.item( r, 0).text()
            L._comp = self.grid_1.item( r, 1).text() 
            L._dd   = float(self.grid_1.item( r, 2).text())
            L._sg   = float(self.grid_1.item( r, 3).text())
            L._rh   = float(self.grid_1.item( r, 4).text())
            # ~ print(r, L._name, L._comp, L._dd, L._sg, L._rh)
            
    def DoYouWantTo(self, action):
        # ~ def msgButtonClick(i):
            # ~ print("Button clicked is:",i.text())
        msgBox = QMessageBox()
        msgBox.setIcon(QMessageBox.Information)
        msgBox.setText("Do you want to "+action+" ... ?")
        msgBox.setWindowTitle("Please confirm")
        msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.No)  #Cancel
        # ~ msgBox.buttonClicked.connect(msgButtonClick)
        returnValue = msgBox.exec()
        if returnValue == QMessageBox.Ok:
            print( "Yes pressed")
            return True
        else:
            print( "No pressed")
            return False
            
    # ~ alert = QMessageBox()
    # ~ alert.setText('layer_load()')
    # ~ alert.exec_()
        
    # ~ self.grid_1.removeRow(row)
    # ~ self.grid_1.insertRow(row)
    # ~ self.grid_1.setCurrentCell(row,0)
    def Select_onerow_warn(self):
        QMessageBox.warning(None,       #parent=
            "Please confirm",           #title=
            "Select one (=1) row!!!",   #text=
            QMessageBox.Ok)             #buttons=
                
    def Select_onerow_get(self):
        row = self.grid_1.currentRow()
        if row <0 : 
            self.Select_onerow_warn()
            return -1
        sel = self.grid_1.selectedItems()
        row0 = self.grid_1.row(sel[0])
        rown = self.grid_1.row(sel[-1])
        if row0 != row or rown != row :
            Self.Select_onerow_warn()
            return -1
        return row
                
    def OnButton_B_Layer_Insert(self):  
        row = self.Select_onerow_get()
        if row < 0 : return
        # ~ else:
        self.LL.insert(row, _XLayer('NoName', 'Si;1', d=10.0, s=1.0, r=2.3))
        self.layer_to_grid()
        
        
    def OnButton_B_Layer_Delete(self):
        row = self.Select_onerow_get()
        if row < 0 : return
        # ~ else:
        del self.LL[row]
        self.layer_to_grid()
    


    def OnButton_B_Layer_Update(self):  
        self.layer_from_grid()
        self.layer_init()
        self.layer_to_grid()
        self.layer_profile()
        self.xrr_layer()
        self._plot()
        
        self.update = self.V_Fit_Update.value()
        self._pdz   = self.V_Profile_Step.value()
        
        print(self.atoms)
        self.Chk_P_atom.clear()
        # ~ for a in self.atoms:
            # ~ self.Chk_P_atom.addItem(a)
        self.Chk_P_atom.addItems(self.atoms)
        # ~ self.Chk_P_atom.SetSelection(0)


    def OnButton_B_Layer_Load(self):  
        if self.DoYouWantTo("Save Results") :
            return
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,
            "QFileDialog.getOpenFileName()", "",
            "XRR Files (*.xrr);;All Files (*)", 
            options=options)
        if not fileName:
            return
        #else:
        # Proceed loading the file chosen by the user
        self.pathname = fileName
        self.layer_load(self.pathname)
        self.layer_parse()
        self.layer_to_grid()
        self.layer_init()
        self.OnButton_B_Layer_Update()

            
    def OnButton_B_Layer_Save(self):  
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,
            "QFileDialog.getSaveFileName()","",
            "XRR Files (*.xrr);;All Files (*)", 
            options=options)
        if not fileName:
            return
         #else:
        # Proceed saving the file chosen by the user in self.directory
        self.layer_save(self.pathname)

            
    def OnButton_B_Data_Load(self):  
        if self.DoYouWantTo("Save Results") :
            return
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,
            "QFileDialog.getOpenFileName()", "",
            "Data Files (*.dat);;All Files (*)", 
            options=options)
        if not fileName:
            return
        #else:
        # Proceed loading the file chosen by the user 
        self.pathname = fileName
        self.T_Data_FName.setText(self.pathname)
        self.data_load(self.pathname)
        self.data_parse()
        self.data_init()
        self.plot_plots[0].plot()
        
    def Menu(self):
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('File')
        editMenu = mainMenu.addMenu('Edit')
        viewMenu = mainMenu.addMenu('View')
        searchMenu = mainMenu.addMenu('Search')
        toolsMenu = mainMenu.addMenu('Tools')
        helpMenu = mainMenu.addMenu('Help')
        
        exitButton = QAction(QIcon('exit24.png'), 'Exit', self)
        exitButton.setShortcut('Ctrl+Q')
        exitButton.setStatusTip('Exit application')
        # ~ exitButton.triggered.connect(self.close)
        exitButton.triggered.connect(self.close)
        fileMenu.addAction(exitButton)
        
        # ~ self.show()
        

    
        
        

class MyTableWidget(QTableWidget):
    data_row = []
    data_col = []
            
    def __init__(self, data, *args):
        QTableWidget.__init__(self, *args)
        self.resize(300,200)
        
        self.setShowGrid(True)
        self.setRowCount(0)
        self.setColumnCount(0)
        self.setData_rows()

    def setData_rows(self, 
        data=[[   'name',     'comp.',        'd [A]',    's [A]',    'g/cm^3',   'g/mol',    'A^3'],
              [   'glass',    'Si;1;O;2',    '0',         '5.0',      '2.6',      '',         ''],
              [   'W',        'W;1',          '10.0',     '5.0',      '20.0',     '',         '']]): 
        self.data_row = data
        nr = len(self.data_row) -1 # 1st row is for header
        nc = len(self.data_row[0])
        self.setRowCount(nr)
        self.setColumnCount(nc)
        self.setHorizontalHeaderLabels(self.data_row[0])
        for r, row in enumerate(self.data_row[1:], start=0):
            for c, item in enumerate(row):
                newitem = QTableWidgetItem(item)
                self.setItem(r, c, newitem)
        self.resizeColumnsToContents()
        self.resizeRowsToContents()
        
    @pyqtSlot()
    def on_click(self):
        print("\n")
        for item in self.tableWidget.selectedItems():
            print(item.row(), item.column(), item.text())
        
    def set_attrib(self):
        self.horizontalHeader().setDefaultSectionSize(80)
        self.horizontalHeader().resizeSection(0, 100)
        self.verticalHeader().setVisible(True)
        self.verticalHeader().setDefaultSectionSize(19)
        self.horizontalHeader().setSectionResizeMode(0, QHeaderView.Interactive)
        self.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.resizeColumnsToContents()
        self.resizeRowsToContents()
        self.setItem(1,2, QTableWidgetItem("Table Cell"))
    
    def setData_cols(self, 
        data={'col1':['1','2','3','4'],
        'col2':['1','2','1','3'],
        'col3':['1','1','2','1']}): 
        self.data_col = data
        horHeaders = []
        for n, key in enumerate(sorted(self.data_col.keys())):
            horHeaders.append(key)
            for m, item in enumerate(self.data[key]):
                newitem = QTableWidgetItem(item)
                self.setItem(m, n, newitem)
        self.setHorizontalHeaderLabels(horHeaders)



if __name__ == '__main__':

    app = QApplication(sys.argv)
    # ~ app = QApplication([])
    # ~ app.setStyle('Fusion')
    # ~ 'Fusion' , 'Windows' , 'WindowsVista' Windows only and 'Macintosh' Mac only.
    
    # ~ app.setStyleSheet("QPushButton { margin: 1ex; }")
    
    # ~ palette = QPalette()
    # ~ palette.setColor(QPalette.ButtonText, Qt.blue)
    # ~ app.setPalette(palette)

    frame = MainFrame()
    sys.exit(app.exec_())
    


"""
# ===============================================================
# The new Stream Object which replaces the default stream 
# associated with sys.stdout
# This object just puts data in a queue!

class WriteStream(object):
    def __init__(self,queue):   self.queue = queue
    def write(self, text):      self.queue.put(text)
    def flush(self):            pass;
        
# A QObject (to be run in a QThread) which sits waiting for data to come through a Queue.Queue().
# It blocks until data is available, and one it has got something from the queue, it sends
# it to the "MainThread" by emitting a Qt Signal 

class MyReceiver(QObject):
    mysignal = pyqtSignal(str)
    def __init__(self,queue,*args,**kwargs):
        QObject.__init__(self,*args,**kwargs)
        self.queue = queue
        
    @pyqtSlot()
    def run(self):
        while True:
            text = self.queue.get()
            self.mysignal.emit(text)
            
'''
    # Within your main window class...
    def __init__(self, parent=None, **kwargs):
        # ...
        # Install the custom output stream
        self.redirect()
        
    def __del__(self):
        # Restore sys.stdout
        sys.stdout = sys.__stdout__

    def redirect(self):
        # Create Queue and redirect sys.stdout to this queue
        self.queue = Queue()
        # Create thread that will listen on the other end of the queue, 
        # ~ and send the text to the textedit in our application
        self.thread = QThread()
        self.my_receiver = MyReceiver(self.queue)
        self.my_receiver.mysignal.connect(self.append_text)
        self.my_receiver.moveToThread(self.thread)
        self.thread.started.connect(self.my_receiver.run)
        self.thread.start()
        sys.stdout = WriteStream(self.queue)
        sys.stderr = WriteStream(self.queue)
'''
# ===============================================================

def _get_format_from_style(self, token, style):
        ''' Returns a QTextCharFormat for token by reading a Pygments style.
        '''
        result = QtGui.QTextCharFormat()
        for key, value in list(style.style_for_token(token).items()):
            if value:
                if key == 'color':
                    result.setForeground(self._get_brush(value))
                elif key == 'bgcolor':
                    result.setBackground(self._get_brush(value))
                elif key == 'bold':
                    result.setFontWeight(QtGui.QFont.Bold)
                elif key == 'italic':
                    result.setFontItalic(True)
                elif key == 'underline':
                    result.setUnderlineStyle(
                        QtGui.QTextCharFormat.SingleUnderline)
                elif key == 'sans':
                    result.setFontStyleHint(QtGui.QFont.SansSerif)
                elif key == 'roman':
                    result.setFontStyleHint(QtGui.QFont.Times)
                elif key == 'mono':
                    result.setFontStyleHint(QtGui.QFont.TypeWriter)
        return result 
        
# ~ class StdoutWrapper(object):
# ~ class StdoutWrapper(QTextEdit):
    # ~ def __init__(self, outwidget):
        # ~ super().__init__() # if using QTextEdit
        # ~ self.widget = outwidget
    # ~ def write(self, s):
        # ~ self.widget.moveCursor(QTextCursor.End)
        # ~ self.widget.insertText( s )
    # ~ def flush(self):
        # ~ pass;
# ~ use:
# ~ sys.stdout = StdoutWrapper(self.tx)
# ~ sys.stderr = StdoutWrapper(self.tx)
# similar for stderr, but you might want an error dialog or make 
# the text stand out using appendHtml
        
'''
import sys
orig_stdout = sys.stdout
f = file('out.txt', 'w')
sys.stdout = f
for i in range(2):
    print ('i = ', i)
sys.stdout = orig_stdout
f.close()
'''

class timer_class():
    def __init__():
        self.timer = QBasicTimer()
        self.step = 0
        
    def timerEvent(self, e):
        if self.step >= 100:
            # ~ self.timer.stop()
            # ~ self.btn.setText('Finished')
            self.step = 0
            return
            
        self.step = self.step + 1
        self.x1.setValue(self.step)
        
    def timeTest():
        if self.timer.isActive():
            self.timer.stop()
            self.step = 0
            self.b1.setText('Fit')
        else:
            self.timer.start(100, self)
            self.b1.setText('Stop')
        alert = QMessageBox()
        alert.setText('data_fit()')
        alert.exec_()
        return
    
        
"""

'''
# https://www.riverbankcomputing.com/static/Docs/PyQt5/signals_slots.html

class on_write_qt_emit(QObject):
    textWritten = pyqtSignal(str, int)
    def __init__(self, cbfunc, cval, *args,**kwargs):
        QObject.__init__(self,*args,**kwargs)
        self.textWritten.connect(cbfunc)
        self.cval = cval # 0:black, 1: red
    def write(self, text):
        self.textWritten.emit(str(text), self.cval)
    def flush(self):       pass;
'''

'''
# https://www.learnpyqt.com/examples/create-desktop-weather-app/
class WorkerSignals(QObject):

    #Defines the signals available from a running worker thread.
 
    finished = pyqtSignal()
    error = pyqtSignal(str)
    result = pyqtSignal(dict, dict)
    
class WeatherWorker(QRunnable):

    #Worker thread for weather updates.

    signals = WorkerSignals()
    is_interrupted = False

    def __init__(self, location):
        super(WeatherWorker, self).__init__()
        self.location = location

    @pyqtSlot()
    def run(self):
        try:
            params = dict(
                q=self.location,
                appid=OPENWEATHERMAP_API_KEY
            )

            url = 'http://api.openweathermap.org/data/2.5/weather?%s&units=metric' % urlencode(params)
            r = requests.get(url)
            weather = json.loads(r.text)

            # Check if we had a failure (the forecast will fail in the same way).
            if weather['cod'] != 200:
                raise Exception(weather['message'])

            url = 'http://api.openweathermap.org/data/2.5/forecast?%s&units=metric' % urlencode(params)
            r = requests.get(url)
            forecast = json.loads(r.text)

            self.signals.result.emit(weather, forecast)

        except Exception as e:
            self.signals.error.emit(str(e))

        self.signals.finished.emit()
        
class MainWindow(QMainWindow, Ui_MainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setupUi(self)
        self.pushButton.pressed.connect(self.update_weather)
        self.threadpool = QThreadPool()
        self.show()
        
    def update_weather(self):
        worker = WeatherWorker(self.lineEdit.text())
        worker.signals.result.connect(self.weather_result)
        worker.signals.error.connect(self.alert)
        self.threadpool.start(worker)

    def alert(self, message):
        alert = QMessageBox.warning(self, "Warning", message)
        
    def weather_result(self, weather, forecasts):
        self.latitudeLabel.setText("%.2f °" % weather['coord']['lat'])
        self.longitudeLabel.setText("%.2f °" % weather['coord']['lon'])

        self.windLabel.setText("%.2f m/s" % weather['wind']['speed'])

        self.temperatureLabel.setText("%.1f °C" % weather['main']['temp'])
        self.pressureLabel.setText("%d" % weather['main']['pressure'])
        self.humidityLabel.setText("%d" % weather['main']['humidity'])

        self.weatherLabel.setText("%s (%s)" % (
            weather['weather'][0]['main'],
            weather['weather'][0]['description']
        )

'''

'''
# https://stackoverflow.com/questions/43964766/pyqt-emit-signal-with-dict
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QApplication


class Emiterer(QtCore.QThread):
    f = QtCore.pyqtSignal(dict)

    def __init__(self):
        super(Emiterer, self).__init__()

    def run(self):
        self.f.emit({"2": {}})
        # self.f.emit({2: {}})  < == this don't work!


class Main(QtWidgets.QMainWindow):

    def __init__(self):
        super(Main, self).__init__()
        self.e = Emiterer()
        self.e.f.connect(self.finised)
        self.e.start()

    def finised(self, r_dict):
        print(r_dict)


if __name__ == '__main__':
    import sys
    app = QApplication(sys.argv)
    m = Main()
    m.show()
    sys.exit(app.exec_())
    


Use object instead of dict in the pyqtSignal definition. E.g.

class Emiterer(QtCore.QThread):
    f = QtCore.pyqtSignal(object)

The reason for this is that signals defined as pyqtSignal(dict) 
are actually interpreted the same as pyqtSignal('QVariantMap') 
by PyQt5 and QVariantMap can only have strings as keys.

You can check this (for your specific class) with

m = Emiterer.staticMetaObject
method = m.method(m.methodCount() - 1)  # the last defined signal or slot
print(method.methodSignature())

This would print PyQt5.QtCore.QByteArray(b'f(QVariantMap)')

'''

'''

class on_write_qt_emit(QObject):
    textWritten = pyqtSignal(str, int)
    def __init__(self, cbfunc, cval, *args,**kwargs):
        QObject.__init__(self,*args,**kwargs)
        self.textWritten.connect(cbfunc)
        self.cval = cval # 0:black, 1: red
    def write(self, text):
        # self.textWritten.emit(str(text), self.cval)
    def flush(self):       pass;
    

    # Within your main window class...
    def __init__(self, parent=None, **kwargs):
        # ...
        # Install the custom output stream
        sys.stdout = on_write_qt_emit(textWritten=self.append2QTextEdit)

    def __del__(self):
        # Restore sys.stdout
        sys.stdout = sys.__stdout__

    def append2QTextEdit(self, text):
        """Append text to the QTextEdit."""
        # Maybe QTextEdit.append() works as well, but this is how I do it:
        cursor = self.textEdit.textCursor()
        cursor.movePosition(QtGui.QTextCursor.End)
        cursor.insertText(text)
        self.textEdit.setTextCursor(cursor)
        self.textEdit.ensureCursorVisible()
'''

'''
print('___________________________________________')
import sys
import inspect
# https://www.oreilly.com/library/view/python-cookbook/0596001673/ch14s08.html
# ~ sys._getframe. 
# ~ This function returns a frame object whose attribute f_code is a code object 
# ~ and the co_name attribute of that object is the function name

this_function_name = sys._getframe(  ).f_code.co_name
this_line_number = sys._getframe(  ).f_lineno
this_filename = sys._getframe(  ).f_code.co_filename

def whoami(  ):
    import sys
    return sys._getframe(1).f_code.co_name
def callersname(  ):
    import sys
    return sys._getframe(2).f_code.co_name


# ~ https://medium.com/@vadimpushtaev/name-of-python-function-e6d650806c4

def fname(func):
    @wraps(func)
    def tmp(*args, **kwargs):
        print('>>', func.__name__ , '<<')
        return func(*args, **kwargs)
    return tmp

fname1 = lambda n=0: sys._getframe(n + 1).f_code.co_name
# for current func name, specify 0 or no argument.
# for name of caller of current func, specify 1.
# for name of caller of caller of current func, specify 2. etc.

fname2 = lambda: inspect.stack()[1][3]

def fname3(): 
    print(inspect.stack()[0][3])
    print(inspect.stack()[1][3]) #will give the caller of foos name, if something called foo
    
    print(inspect.stack()[0].function, 'inspect.stack[0].function')
    print(inspect.stack()[1].function, 'inspect.stack[1].function') 
    
    print( fname1(),  'currentFuncName() with sys._getframe(n)')
    print( fname1(1), 'currentFuncName(1)' )   
    
    print(inspect.currentframe().f_code.co_name, 'inspect.currentframe().f_code.co_name')
    print(inspect.currentframe().f_back.f_code.co_name, 'inspect.currentframe().f_back.f_code.co_name')
    # ~ print( "my name is '{}'".format( variable )

@fname
def my_funky_name():
    print( "STUB" )
    fname1()
    fname2()
    fname3()

my_funky_name()

# ~ $ python -m timeit -s 'import inspect, sys' 'inspect.stack()[0][0].f_code.co_name'
# ~ 1000 loops, best of 3: 499 usec per loop
# ~ $ python -m timeit -s 'import inspect, sys' 'inspect.stack()[0][3]'
# ~ 1000 loops, best of 3: 497 usec per loop
# ~ $ python -m timeit -s 'import inspect, sys' 'inspect.currentframe().f_code.co_name'
# ~ 10000000 loops, best of 3: 0.1 usec per loop
# ~ $ python -m timeit -s 'import inspect, sys' 'sys._getframe().f_code.co_name'
# ~ 10000000 loops, best of 3: 0.135 usec per loop

def safe_super(_class, _inst):
    """safe super call"""
    try:
        return getattr(super(_class, _inst), _inst.__fname__)
    except:
        return (lambda *x,**kx: None)


def with_name(function):
    def wrap(self, *args, **kwargs):
        self.__fname__ = function.__name__
        return function(self, *args, **kwargs)
    return wrap

class A(object):

    def __init__():
        super(A, self).__init__()

    @with_name
    def test(self):
        print 'called from A\n'
        safe_super(A, self)()

class B(object):

    def __init__():
        super(B, self).__init__()

    @with_name
    def test(self):
        print 'called from B\n'
        safe_super(B, self)()

class C(A, B):

    def __init__():
        super(C, self).__init__()

    @with_name
    def test(self):
        print 'called from C\n'
        safe_super(C, self)()

# ~ testing it :

a = C()
a.test()

# ~ Inside each @with_name decorated method you have access to 
# ~ self.__fname__ as the current function name.

def foo():
    """foo docstring"""
    print(eval(sys._getframe().f_code.co_name).__doc__)
    
foo()


# ~ print('___________________________________________')
'''

        # ~ radiobutton = QRadioButton("Australia")
        # ~ radiobutton.setChecked(True)
        # ~ radiobutton.country = "Australia"
        # ~ radiobutton.toggled.connect(self.onClicked)
        # ~ self.show()

# ~ class MyTabWidget(QWidget):
    
    ##def __init__(self, parent):
        ##super(QWidget, self).__init__(parent)
    # ~ def __init__(self):
        # ~ super(QWidget, self).__init__()
        
        # Initialize tab screen
        # ~ self.tabs = QTabWidget()
        # ~ self.tabs.resize(300,200)
        # ~ self.tabs.setTabPosition(2) # Left-side
        
        # Add tabs to widget
        # ~ self.layout = QVBoxLayout(self)
        # ~ self.layout.addWidget(self.tabs)
        # ~ self.setLayout(self.layout)
        
        # Add tabs
        # ~ self.tab1 = QWidget()
        # ~ self.tab2 = QWidget()
        # ~ self.tab3 = QWidget()
        # ~ self.tabs.addTab(self.tab1,"I vs 2theta")
        # ~ self.tabs.addTab(self.tab2,"I*Q^4 vs Q")
        # ~ self.tabs.addTab(self.tab3,"Layer profile")
        
        # Create first tab
        # ~ self.lay1 = QVBoxLayout(self)
        # ~ self.tab1.setLayout(self.lay1)
        # ~ self.cnv1 = PlotCanvas(self, title='c1')
        # ~ self.nav1 = NavigationToolbar(self.cnv1, self)
        # ~ self.lay1.addWidget(self.cnv1)
        # ~ self.lay1.addWidget(self.nav1)
        
        
        # Create second tab
        # ~ self.lay2 = QVBoxLayout(self)
        # ~ self.tab2.setLayout(self.lay2)
        # ~ self.cnv2 = PlotCanvas(self, title='c2')
        # ~ self.nav2 = NavigationToolbar(self.cnv2, self) 
        # ~ self.lay2.addWidget(self.cnv2)
        # ~ self.lay2.addWidget(self.nav2)
        

        # Create third tab
        # ~ self.lay3 = QVBoxLayout(self)
        # ~ self.tab3.setLayout(self.lay3)
        # ~ self.cnv3 = PlotCanvas(self, title='c2')
        # ~ self.nav3 = NavigationToolbar(self.cnv3, self) 
        # ~ self.lay3.addWidget(self.cnv3)
        # ~ self.lay3.addWidget(self.nav3)
