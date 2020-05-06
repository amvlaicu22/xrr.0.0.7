from PyQt5 import QtWidgets, uic
import sys

class Ui(QtWidgets.QMainWindow):
    def __init__(self):
        super(Ui, self).__init__()
        uic.loadUi('xrr_qt2gui.ui', self)
        self.show()
        # ~ self.B_Fit_Start = self.findChild(QtWidgets.QPushButton, 'B_Fit_Start')
        self.B_Fit_Start.released.connect(self.B_Fit_Start_) 
        # Remember to pass the definition/method, not the return value!
    
    def B_Fit_Start_(self):
        print("B_Fit_Start")

app = QtWidgets.QApplication(sys.argv)
window = Ui()
app.exec_()
