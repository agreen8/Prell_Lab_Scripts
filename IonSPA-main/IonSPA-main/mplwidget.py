# Imports
from PyQt5 import QtWidgets
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as Canvas, NavigationToolbar2QT as NavigationToolbar

import matplotlib

# Ensure using PyQt5 backend
matplotlib.use('QT5Agg')

# Matplotlib canvas class to create figure
class MplCanvas(Canvas):
    def __init__(self):
        self.fig = Figure()
        self.fig.set_size_inches(10,7)
        self.axes = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)
        super(MplCanvas, self).setSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
        super(MplCanvas, self).updateGeometry()

# Matplotlib widget
class MplWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super(MplWidget, self).__init__(parent)   # Inherit from QWidget
        self.canvas = MplCanvas()                 # Create canvas object

        # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
        toolbar = NavigationToolbar(self.canvas, self)

        self.vblayout = QtWidgets.QVBoxLayout()              # Set box for plotting
        self.vblayout.addWidget(toolbar)
        self.vblayout.addWidget(self.canvas)
        self.setLayout(self.vblayout)
        self.initgeom = self.geometry()

    # self.plot = self.canvas.axes.plot()



