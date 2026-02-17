import sys
from os import path
import json
import pandas as pd

from PyQt5 import QtWidgets, QtCore
        # , QtGui #pyqt stuff

QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True) #enable highdpi scaling
QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True) #use highdpi icons

from spa_avggui import Ui_spa_avg_qt
from ionspa import loadjfile    #, fracloss, lossrate, fracremains
import spa_avg as spaavg
from ionspa import makecell
import argparse

celld_def = {'type': 'cell', 'gas':'N2', 'L': 0.18, 'E': 38.89, 'pressure': 0.00002416, 'T': 300.0, 'Voffset': 0.0}
iond_def = {'name': 'peptide', 'mass': 17.594, 'charge': 9, 'num_atoms': 2411, 'CCS': 20.2, 'T0': 300.0, 'hcprofname': 'peptide'}
rund_def = {}

def getplotinfo(tlogger, xaxis, yaxis):
    if xaxis == 't':
        xvals = tlogger.t
        pltxlabel = 't [s]'
    elif xaxis == 'coll':
        xvals = tlogger.coll
        pltxlabel = 'collisions'
    elif xaxis == 'zp':
        xvals = tlogger.zp
        pltxlabel = 'z position [m]'

    # add other possible y values options here
    if yaxis == 'Ti':
        yvals = tlogger.Ti
        pltylabel = 'T [K]'
    elif yaxis == 'coll':
        yvals = tlogger.coll
        pltylabel = 'collisions'
    elif yaxis == 'KE':
        yvals = tlogger.KE
        pltylabel = 'KE [eV]'
    elif yaxis == 'ui':
        yvals = tlogger.ui
        pltylabel = 'ui [eV]'
    elif yaxis == 'chi':
        yvals = tlogger.chi
        pltylabel = 'chi'
    elif yaxis == 'zp':
        yvals = tlogger.zp
        pltylabel = 'z position [m]'

    return xvals, yvals, pltxlabel, pltylabel
'''
paramset = {'paramd': {'cell': {'type': 'cell', 'Mgamu': 28.0, 'L': 0.18, 'E': 38.89, 'pressure': 0.02416, 'T': 300.0, 'Voffset': 0.0}, 'ion': {'name': 'peptide', 'mass': 17.594, 'charge': 9.0, 'num_atoms': 2411, 'CCS': 20.2, 'T0': 300.0, 'hcprofname': 'peptide'}, 'injectV': 10.5, 'option': 'IICT', 'plotfile': 'none'}, 'plotax': {'xaxis': 't', 'yaxis': 'Ti'}}
plot_data = [
    {
        'params': json.dumps(paramset),
        'coll': [1, 2, 3],
        't': [0, 1, 2],
        'zp': [0.1, 0.2, 0.3],
        'KE': [10, 20, 30],
        'ui': [5, 15, 25],
        'Ti': [100, 200, 300],
        'chi': [0.01, 0.02, 0.03]
    },
    {
        'params': json.dumps(paramset),
        'coll': [1, 2, 3, 4],
        't': [0, 1, 2, 3],
        'zp': [0.1, 0.2, 0.3, 0.4],
        'KE': [10, 20, 30, 40],
        'ui': [5, 15, 25, 35],
        'Ti': [100, 200, 300, 400],
        'chi': [0.01, 0.02, 0.03, 0.04]
    }
]
'''
class MainWindow(QtWidgets.QMainWindow, Ui_spa_avg_qt):
    def __init__(self, *args, obj=None, **kwargs):
        # with pre-saved sp_avg_qtwin.py made with command line pyuic5 tool
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setupUi(self)

        # or do this direct from the ui file
#        from PyQt5.QtWidgets import QUiLoader
#        self.ui = QUiLoader().load('spa_avg_qtsrc.ui')  # Load your UI file
#        self.setCentralWidget(self.ui)
        # intended to give same results either way

        self.outlines = ["output text:"]
        self.plotdata = []
        self.butn_plot.clicked.connect(self.onplot)
        self.butn_addplot.clicked.connect(self.onaddplot)
        self.plaxes = self.plot_widget.canvas.axes
        self.cvsizedelta = None
        self.setWindowTitle(__file__)
        self.setctrlvalues(gparamd)
        self.actionOpen_cell_ion_params.triggered.connect(self.open_cellionparams)
        self.actionExport_plot_data.triggered.connect(self.export_plotdata)
        self.actionSave_run_params.triggered.connect(self.save_params)

#        print(help(self.actionOpen_cell_ion_params.checked))
#        print(dir(QtWidgets.QAction))
#        qw = QtWidgets.QAction
        

#        print(dir(self.actionExport_plot_data))

    def save_params(self):
        fname, filt = QtWidgets.QFileDialog.getSaveFileName(caption='Save params', filter='txt file (*.txt)')
        if fname == '':
            return
        print(f'save params requested to {repr(fname)}')

        paramd, plotax = self.getctrlvalues()
        with open(fname, 'w') as fp:
            json.dump(paramd, fp, indent=4)       # probably needs formatting to match input file structure

    def open_cellionparams(self):
        fnames, filt = QtWidgets.QFileDialog.getOpenFileNames(caption='Open cell/ion file(s)')
        # more than one can be chosen.  returns a list of full absolute path names.
        if len(fnames) < 1:
            return
        paramd = loadjfile(fnames)
        self.setctrlvalues(paramd)
        print(f'open requested with {fnames}\n loaded: {paramd}')


    def export_plotdata(self):
        fname, filt = QtWidgets.QFileDialog.getSaveFileName(caption='Save plot data', filter='csv file (*.csv)')
#        base,ext = path.splitext(fname)
        # the filter will ensure .csv extension unless user specifies differently.
        if fname == '':
            return

        print(f'export plot data requested.   file = {repr(fname)}')
        
        # rearrange the plotdata list (1 item for each plot) into a list of dicts.
        # the dicts have keys: plotname plus the elements from each plotdata item.
        # a header line onl includes plotname and plotdata[i]['params']
        # then N rows with plotname, empty params field and corresponding plotdata[i][field][j] values.
        # this is then assembled to a dataframe and exported using df.to_csv(fname, index=False)

        # Prepare data for export
        export_data = []
        for i, plot in enumerate(self.plotdata):
            plot_name = f'plot{i+1}'
            # Add a header row for the plot parameters
            export_data.append({
                'plot': plot_name,
                'params': plot['params'],
                'coll': '',
                't': '',
                'zp': '',
                'KE': '',
                'ui': '',
                'Ti': '',
                'chi': ''
            })
            # Add the data rows for the plot
            for j in range(len(plot['coll'])):
                row = {
                    'plot': plot_name,
                    'params': '',
                    'coll': plot['coll'][j],
                    't': plot['t'][j],
                    'zp': plot['zp'][j],
                    'KE': plot['KE'][j],
                    'ui': plot['ui'][j],
                    'Ti': plot['Ti'][j],
                    'chi': plot['chi'][j]
                }
                export_data.append(row)
        # Create DataFrame
        df = pd.DataFrame(export_data)
        # Write to CSV file
        df.to_csv(fname, index=False)




    def resizeEvent(self, event):
        self.resizeFunction()  # Your custom logic here
        super().resizeEvent(event)  # Call the base class implementation
        
    def printsizes(self, headtxt):
        sw = self.centralwidget.width(); sh = self.centralwidget.height()
        pcw = self.plot_widget.canvas.width(); pch = self.plot_widget.canvas.height()
        pww = self.plot_widget.width(); pwh = self.plot_widget.height()
        print(f'{headtxt}:  s:[{sw},{sh}]   pc:[{pcw},{pch}]   pw:[{pww},{pwh}]')
    
    def resizeFunction(self):
        # Resize the Matplotlib canvas (MplWidget) to match the window size
        #self.ui.MplWidget.setGeometry(0, 0, self.width(), self.height() - 20)
        if self.cvsizedelta:
            cw = self.width() - self.cvsizedelta[0]
            ch = self.height() - self.cvsizedelta[1]
            self.plot_widget.canvas.resize(cw,ch)
            parent = self.plot_widget.canvas.parentWidget()
            parent.resize(cw+18, ch+58)
#        self.printsizes('resize')

#    def resize(self, xsize, ysize):
#        #print(f'{self.centralwidget}')
#        print(f'caught a resize: [{xsize}, {ysize}]')

    def setctrlvalues(self, paramd):
#        print(gparamd)
        jsonfile = paramd.get('jsonfile', None)
        if jsonfile:
            print('jsonfile:', jsonfile)
        for key, value in paramd.items():
            print(key, value)
        celld = paramd.get('cell', celld_def)
        print('cell: ', celld)
        iond = paramd.get('ion', iond_def)
        print('ion:', iond)

        if celld:
            # let cell class decode pressure vs density input
            cc = makecell(celld)
            pressure = cc.P * 1000.0    # our pressure is in mbar, not bar
            
            gas = celld.get('gas', None)
            glist = 'Ar He N2'.split()
            if gas not in glist:
                cMgamu = celld.get('Mgamu', 28)
                if cMgamu == 40: gas = 'Ar'
                elif cMgamu == 28: gas = 'N2'
                elif cMgamu == 4: gas = 'He'
                else: gas = 'Ar'

            self.label_ctypevalue.setText(celld['type'])
            self.c_gas.setCurrentText(gas)
            self.c_L.setText(str(celld['L']))
            self.c_E.setText(str(celld['E']))
            self.c_P.setText(f"{pressure:0.4g}")
            self.c_Tg.setText(str(celld['T']))
            self.c_Voff.setText('0')                # we might add this to command line param if desired.

        if iond:
            self.i_name.setText(iond['name'])
            self.i_Mi.setText(str(iond['mass']))
            self.i_charge.setText(str(iond['charge']))
            self.i_atoms.setText(str(iond['num_atoms']))
            self.i_CCS.setText(str(iond['CCS']))
            self.i_T0.setText(str(iond['T0']))
            self.i_hctype.setCurrentText(iond['hcprofname'])          #setCurrentIndex(5)

    def getctrlvalues(self):
        celld = {
            "type": self.label_ctypevalue.text(),
            "gas": self.c_gas.currentText(),
            "L": float(self.c_L.text()), 
            "E": float(self.c_E.text()),
            "pressure": float(self.c_P.text())/1000.0,
            "T": float(self.c_Tg.text()),
            "Voffset": float(self.c_Voff.text())
        }
        iond = {
            "name": self.i_name.text(),
            "mass": float(self.i_Mi.text()),
            "charge": float(self.i_charge.text()),
            "num_atoms": int(self.i_atoms.text()),
            "CCS": float(self.i_CCS.text()),
            "T0": float(self.i_T0.text()),
            "hcprofname": self.i_hctype.currentText()
        }
        paramd = {
            "cell": celld,
            "ion": iond,
            "injectV": float(self.run_injectV.text()),
            "option": self.run_calcoption.currentText(),
            "plotfile": "none"
        }
        plotax = {
            "xaxis": self.run_xaxis.currentText(),
            "yaxis": self.run_yaxis.currentText()
        }

        return paramd, plotax
    
    def onplot(self):
        self.plaxes.cla()
        self.outlines = []
        self.plotdata = []
        self.doplot()
    
    def onaddplot(self):
        self.doplot()

    def doplot(self):
        #print(dir(self))
        paramd, plotax = self.getctrlvalues()
        pds = dict(paramd=paramd, plotax=plotax)
        self.outlines.append(f'{pds}')
        self.statustext.setPlainText('\n'.join(self.outlines))
        xaxis = plotax['xaxis']
        yaxis = plotax['yaxis']
        ionname = paramd['ion']['name']
        retvals, postdlist = spaavg.main(**paramd)
        self.plotdata.append( dict(params=json.dumps(pds), 
                coll=retvals.coll, 
                t=retvals.t,
                zp=retvals.zp,
                Ti=retvals.Ti,
                KE=retvals.KE,
                ui=retvals.ui,
                chi=retvals.chi)
            )
        xvals, yvals, pltxlabel, pltylabel = getplotinfo(retvals, xaxis, yaxis)
        # to save plot date, we need to save the above 4 items, xvals, yvals, xlabel, ylabel
        # in addition, each set should be linked to the pds dict for the parameters used
        # and we can order by plot number.  a possible structure:
        # list(dict(plotnum=i, plotparms=pds, xvals, yvals, xlabel, ylabel))
        # and the list is cleared when called from onplot and appended in both onplot and onaddplot.
        # exporting to a cvs or other format is TBD
        # some generalized ways of doing this are limited: plots might have 
        # been added with different xlabel, ylabel combinations.  Not generally useful, 
        # but there might be a reason to keep this.
        
        self.plaxes.plot(xvals, yvals)
        self.plaxes.set_title(ionname)
        self.plaxes.set_xlabel(pltxlabel)
        self.plaxes.set_ylabel(pltylabel)
        self.plot_widget.canvas.draw()
#        self.plot_panel.plot(xvals, yvals, title=ionname, xlabel=pltxlabel, ylabel=pltylabel)

gparamd = None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='python spa_avg_qt.py',
        description='''Launch the GUI for spa_avg.  An optional json input file may be specified to initialize the cell and ion parameters.''')
    parser.add_argument('jsonfile', nargs='*')
#    parser.add_argument('--vin', default=10) #   , nargs=1)
#    parser.add_argument('--xname', choices=['t', 'z', 'collnum'], default='t')
    args = parser.parse_args()  # uses sys.argv internally

    gparamd = loadjfile(args.jsonfile)  # can be used to get initial cell, ion values
#    gparamd.update(colname=args.vname, Vin=args.vin)


    app = QtWidgets.QApplication(sys.argv)

    window = MainWindow()
    window.show()
    # window.printsizes('mains')  # initial size adjustments were done with show()
    dw = window.width() - window.plot_widget.canvas.width()
    dh = window.height() - window.plot_widget.canvas.height()
    window.cvsizedelta = (dw,dh)
    # print(f'cvdelta: {window.cvsizedelta}')
    app.exec()