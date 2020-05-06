#!/usr/bin/python3
# -*- coding: utf-8 -*-
# -*- coding: iso-8859-15 -*-
# -*- coding: latin-1 -*-
# -*- coding: ascii -*-
# ~ from __future__ import print_function
import numpy as np
import time, sys
from threading import Thread
import wx
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar

from mpl_toolkits.mplot3d import Axes3D

from XRR import XRR     #, _XLayer
from xrr_wx2gui import XRR_Frame
# ~ from xrr_wx2xrc import xrc_MainFrame

class Panel_Plot(wx.Panel):
    def __init__(self, parent=None, id=-1, dpi=None, figsize=(3,2), 
        tbar = 0, row=1, col=1, fplots=None):     #, **kwargs):
        wx.Panel.__init__(self, parent=parent, id=id)  
        self.parent = parent
        self.row = row
        self.col = col
        self.plt = fplots
        
        self.fig = mpl.figure.Figure(dpi=dpi, figsize=figsize)
        self.cvs = FigureCanvas(self, -1, self.fig)
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

        self._sizer = wx.BoxSizer(wx.VERTICAL)  # HORIZONTAL (bad!)
        self._sizer.Add(self.cvs, 1, wx.EXPAND)
        if tbar: 
            self.toolbar = NavigationToolbar(self.cvs)
            self.toolbar.Realize()
            self._sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        self.SetSizer(self._sizer)
        
        self.sizer_ = wx.BoxSizer(wx.VERTICAL)
        self.sizer_.Add(self, -1, wx.ALL|wx.EXPAND, 5)
        self.parent.SetSizer(self.sizer_)
        
    def plot(self):
        axs = np.asarray(self.axs)
        for fplot, axs in zip(self.plt, axs.flat): fplot(axs)
        self.draw()
    def draw(self):
        # ~ self.fig.tight_layout()
        # ~ self.fig.canvas.draw()
        # ~ self.fig.canvas.flush_events()
        self.cvs.draw()




    
    
# GUI Frame class that spins off the worker thread
# ~ class MainFrame(wx.Frame):
class MainFrame( XRR_Frame, XRR):
# ~ class MainFrame(xrc_MainFrame, XRR):

    class _Socket(wx.PyEvent):
        def __init__(self, dest=None, func=None):
            wx.PyEvent.__init__(self)  
            self.id = wx.NewId()
            self.SetEventType(self.id)
            if dest and func:
                self.connect(dest, func)
        def connect(self, dest, func):
            self.dest = dest
            self.func = func
            self.dest.Connect(-1, -1, self.id, self.func )
            return self
        def emit(self, data): 
            self.write(data)
        def write(self, data): 
            self.data = data
            wx.PostEvent(self.dest, self)    
        def flush(self): pass
    
    def __init__(self, *args, **kwds):
        # ~ super().__init__(*args, **kwds)
        XRR_Frame.__init__(self, *args, **kwds)
        # ~ xrc_MainFrame.__init__(self, None)
        self.Show()
        XRR.__init__(self)
        
        self.print_2sys()
        self.contentNotSaved = True
        self.SetTitle("X-ray reflectivity - Parrat algorithm")    
        
        self.layer_to_grid()
        self._plot()
        
        print(self.atoms)
        self.Chk_P_atom.Clear()
        for a in self.atoms:
            self.Chk_P_atom.Append(a)
        self.Chk_P_atom.SetSelection(0)
        
        self.V_Fit_Update.SetValue(self.update) 
        self.update = self.V_Fit_Update.GetValue()
        
        self._pdz   = self.V_Profile_Step.GetValue()
        self.fitprof = self.RB_Fit_Profile.GetValue()
        

        import  wx.html2
        # ~ self.help_html = wx.html.HtmlWindow(self.P_Help)
        self.help_html = wx.html2.WebView.New(self.notebook_2_Help)
        sizer_16 = wx.BoxSizer(wx.VERTICAL)
        sizer_16.Add(self.help_html, 1, wx.EXPAND, 0)
        self.notebook_2_Help.SetSizer(sizer_16)
        self.help_html.SetPage(self.help_str, "")
        # ~ html.LoadPage(help_str)
        
        
        
        self.print_2socket()
        print('test')
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        



        
    def print_out(self, event):
        font = wx.Font( wx.FontInfo(8).Family(wx.TELETYPE) )
        self.T_Logger.SetForegroundColour(wx.BLACK)
        self.T_Logger.SetFont(font)
        self.T_Logger.write(event.data)
        
    def print_err(self, event):
        font = wx.Font( wx.FontInfo(8).Bold().Family(wx.TELETYPE) )
        self.T_Logger.SetForegroundColour(wx.RED)
        self.T_Logger.SetFont(font)
        self.T_Logger.write(event.data)
        
        # ~ wx.CallAfter(self.T_Logger.write,     event)    #thread safe
        # ~ wx.CallAfter(self.T_Logger.WriteText, event)    #thread safe
        # ~ wx.Yield()

        # ~ font = wx.Font(pointSize=10, 
            # ~ family=wx.MODERN, # DEFAULT, DECORATIVE, ROMAN, SCRIPT, SWISS, MODERN, TELETYPE, MAX
            # ~ style=wx.NORMAL, # FONTSTYLE_ITALIC, FONTSTYLE_SLANT, FONTSTYLE_MAX
            # ~ weight=wx.NORMAL, # THIN, EXTRALIGHT, LIGHT, NORMAL, MEDIUM, SEMIBOLD, BOLD, EXTRABOLD, HEAVY, EXTRAHEAVY, MAX
            # ~ underline=False, 
            # ~ faceName=u'Console',
            # ~ encoding=FONTENCODING_DEFAULT)
        
        # Create a font using wx.FontInfo
        # ~ FontInfo:
        # ~ AntiAliased     Set anti-aliasing flag.
        # ~ Bold            Use a bold version of the font.
        # ~ Encoding        Set the font encoding to use.
        # ~ FaceName        Set the font face name to use.
        # ~ Family          Set the font family.
        # ~ GetWeightClosestToNumericValue        Get the symbolic weight closest to the given raw weight value.
        # ~ Italic        Use an italic version of the font.
        # ~ Light        Use a lighter version of the font.
        # ~ Slant        Use a slanted version of the font.
        # ~ Strikethrough        Use a strike-through version of the font.
        # ~ Style        Specify the style of the font using one of FontStyle constants.
        # ~ Underlined        Use an underlined version of the font.
        # ~ Weight          Specify the weight of the font.

        # ~ wx.FONTFAMILY_DEFAULT        Chooses a default font.
        # ~ wx.FONTFAMILY_DECORATIVE    A decorative font.
        # ~ wx.FONTFAMILY_ROMAN         A formal, serif font.
        # ~ wx.FONTFAMILY_SCRIPT        A handwriting font.
        # ~ wx.FONTFAMILY_SWISS         A sans-serif font.
        # ~ wx.FONTFAMILY_MODERN        A fixed pitch font.
        # ~ wx.FONTFAMILY_TELETYPE      A teletype (i.e. monospaced) font.
        # ~ wx.FONTFAMILY_MAX

        # ~ font.setWeighht(wx.FONTWEIGHT_BOLD if cval else wx.FONTWEIGHT_NORMAL)
        # ~ wx.FONTWEIGHT_THIN        Thin font (weight = 100).
        # ~ wx.FONTWEIGHT_EXTRALIGHT    Extra Light (Ultra Light) font (weight = 200).
        # ~ wx.FONTWEIGHT_LIGHT         Light font (weight = 300).
        # ~ wx.FONTWEIGHT_NORMAL    Normal font (weight = 400).
        # ~ wx.FONTWEIGHT_MEDIUM    Medium font (weight = 500).
        # ~ wx.FONTWEIGHT_SEMIBOLD    Semi Bold (Demi Bold) font (weight = 600).
        # ~ wx.FONTWEIGHT_BOLD    Bold font (weight = 700).
        # ~ wx.FONTWEIGHT_EXTRABOLD    Extra Bold (Ultra Bold) font (weight = 800).
        # ~ wx.FONTWEIGHT_HEAVY    Heavy (Black) font (weight = 900).
        # ~ wx.FONTWEIGHT_EXTRAHEAVY    Extra Heavy font (weight = 1000).
        # ~ wx.FONTWEIGHT_MAX
        


    def fit_upd(self, event):  # wxGlade: XRR_Frame.<event_handler>
        # ~ XRR.fit_upd(self, event)
        print( f'\nxrr>>: update {event.data}', flush=True)
        if self.fitprof == 0:
            self.layer_profile()
        self.layer_to_grid()
        self._plot()

        
    def fit_run(self, event):  # wxGlade: XRR_Frame.<event_handler>
        XRR.fit_run(self, event)
        self.B_Fit_N.SetLabel(str(self.fitn))

        
    def fit_end(self, event):  # wxGlade: XRR_Frame.<event_handler>
        # ~ XRR.fit_end(self, event)
        print( f'\nxrr>>: finish {event.data}', flush=True)
        self.fit_worker = None
        self.B_Fit_Start.Enable()
        self.fit_upd(event)


        
    def fit_prestart(self, event):  # wxGlade: XRR_Frame.<event_handler>
        self.B_Fit_N.SetLabel('0')
        self.T_Logger.Clear()

        self.layer_from_grid()
        self.layer_init()
        self.layer_to_grid()
        
        fitkeys = []
        if self.Chk_L_ab.GetValue():            fitkeys.append('ab')
        if self.Chk_L_dd.GetValue():            fitkeys.append('dd')
        if self.Chk_L_rh.GetValue():            fitkeys.append('rh')
        if self.Chk_L_sg.GetValue():            fitkeys.append('sg')
        
        self.fit_set_fitkeys(fitkeys)
        
        if self.RB_Fit_LM.GetValue():       self.fitmode = 0
        if self.RB_Fit_DE.GetValue():       self.fitmode = 1
        if self.RB_Fit_Layer.GetValue():    self.fitprof = 0
        if self.RB_Fit_Profile.GetValue():  self.fitprof = 1
        
        self.update = self.V_Fit_Update.GetValue()
        self._pdz   = self.V_Profile_Step.GetValue()
        # ~ if p == 1 : profile
        fitkeysP = []
        if self.Chk_P_ab.GetValue():        fitkeysP.append('ab')
        if self.Chk_P_rh.GetValue():        fitkeysP.append('rh')
        if self.Chk_P_ac.GetValue():        fitkeysP.append('at')
        atom = self.Chk_P_atom.GetValue()
        # ~ self.button_fit_start.SetLabel(' ... ')
        self.B_Fit_Start.Disable()
        # ~ myobject = event.GetEventObject()
        # ~ myobject.Disable()
        

                
        
    
        
    ################################################################
    
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
        # ~ else:
        # ~ for p in self.plot_plots: p.plot()
        
    def layer_to_grid(self):
        # ~ print_(where())
        self.grid_1.ClearGrid()
        nr = self.grid_1.GetNumberRows()
        self.grid_1.DeleteRows( pos=0, numRows=nr)  #, updateLabels=True )
        self.grid_1.SetDefaultCellAlignment(wx.ALIGN_RIGHT, wx.ALIGN_CENTRE)
        self.grid_1.SetRowLabelSize(20)
        
        row = 0
        for L in self.LL:
            self.grid_1.AppendRows(numRows=1)   #, updateLabels=True)
            self.grid_1.SetRowLabelValue(row, str(row))
            self.grid_1.SetCellValue( row, 0, L._name)
            self.grid_1.SetCellValue( row, 1, L._comp) 
            self.grid_1.SetCellValue( row, 2, "{:.1f}".format(L._dd)) 
            self.grid_1.SetCellValue( row, 3, "{:.1f}".format(L._sg)) 
            self.grid_1.SetCellValue( row, 4, "{:.4f}".format(L._rh)) 
            self.grid_1.SetCellValue( row, 5, "{:.1f}".format(L._Mm)) 
            self.grid_1.SetCellValue( row, 6, "{:.1f}".format(L._Vm))
            row = row +1
        # ~ self.grid_1.AutoSize()
    
    def layer_from_grid(self):
        r = 0
        for L in self.LL:
            L._name  = self.grid_1.GetCellValue( r, 0)
            L._comp = self.grid_1.GetCellValue( r, 1) 
            L._dd   = float(self.grid_1.GetCellValue( r, 2))
            L._sg   = float(self.grid_1.GetCellValue( r, 3))
            L._rh   = float(self.grid_1.GetCellValue( r, 4))
            r = r + 1
            # ~ print(L._name, L._comp, L._dd, L._sg, L._rh)


            
    def wxFileDialog(self, ftype = ".abc", save=0): # 0 = open, 1 = save
        d = {0:"Load", 1:"Save"}
        msg = d[save] + " *"+ftype+" file"
        descript = "Files (*"+ftype+")"
        wildc = descript+"|*"+ftype
        style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
        if save is 1:
            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
            
        with wx.FileDialog(self, message = msg, 
            wildcard=wildc, style=style) as fileDialog:
            answ = fileDialog.ShowModal()
            self.pathname    = fileDialog.GetPath()
            self.directory   = fileDialog.GetDirectory()
            self.filename    = fileDialog.GetFilename()
        if answ == wx.ID_CANCEL:   return 0
        else:
            if save is 0: 
                return  1
            else:
                if ftype not in self.filename:
                    self.filename += ftype
                    self.pathname += ftype                    
                return 1
    
    def OnValue_V_Fit_Update(self, event):  # wxGlade: XRR_Frame.<event_handler>
        # ~ print("Event handler 'OnValue_V_Fit_Update' not implemented!")
        self.update = self.V_Fit_Update.GetValue()
        print("New update=", self.update)
        event.Skip()
        
    def OnButton_B_Data_Load(self, event):  # wxGlade: XRR_Frame.<event_handler>
        # ~ print("Event handler 'ev_data_load' not implemented!")
        if self.DoYouWantTo("Save Results") :
            return
        if not self.wxFileDialog(ftype = ".dat", save=0) :
            event.Skip()
            return
        #else:
        # Proceed loading the file chosen by the user in self.directory
        self.T_Data_FName.write(self.pathname)
        self.data_load(self.pathname)
        self.data_parse()
        self.data_init()
        self.plot_plots[0].plot()
        # ~ self.data_plot1()
        # ~ self.data_plot2()
        event.Skip()

    def OnButton_B_Layer_Load(self, event):  # wxGlade: XRR_Frame.<event_handler>
        # ~ print("Event handler 'ev_layers_load' not implemented!")
        if not self.wxFileDialog(ftype = ".xrr", save=0) :
            event.Skip()
            return
        #else:
        # Proceed loading the file chosen by the user in self.directory
        self.layer_load(self.pathname)
        self.layer_parse()
        self.layer_to_grid()
        self.layer_init()
        self.layer_profile()
        self.layer_plot()
        event.Skip()

    def OnButton_B_Layer_Save(self, event):  # wxGlade: XRR_Frame.<event_handler>
        # ~ print("Event handler 'ev_layers_save' not implemented!")
        if not self.wxFileDialog(ftype = ".xrr", save=1) :
            event.Skip()
            return
        #else:
        # Proceed saving the file chosen by the user in self.directory
        self.layer_save(self.pathname)
        event.Skip()
        
    def DoYouWantTo(self, action):
        quest = "Do you want to "+action+" ... ?"
        title = "Please confirm"
        parent = None
        result =  wx.MessageBox(quest, title, wx.YES_NO, parent)
        # ~ result =  wx.MessageDialog(
            # ~ parent, quest, title, wx.YES_NO).ShowModal()
        if result == wx.YES:
            print( "Yes pressed")
            return True
        else:
            print( "No pressed")
            return False

    def Select_onerow_warn(self):
        quest = "Select one (=1) row!!!"
        title = "Please confirm"
        parent = None
        result =  wx.MessageBox(quest, title, wx.YES, parent)
                
    def Select_onerow_get(self):
        rows = self.grid_1.GetSelectedRows()
        if len(rows) != 1 : 
            self.select_onerow_warn()
            return -1
        return rows[0]

    def OnButton_B_Layer_Insert(self, event):  
        row = self.Select_onerow_get()
        if row < 0 : return
        # ~ else:
        self.LL.insert(row, _XLayer('NoName', 'Si;1', d=10.0, s=1.0, r=2.3))
        self.layer_to_grid()
        event.Skip()

    def OnButton_B_Layer_Delete(self, event):  
        row = self.Select_onerow_get()
        if row < 0 : return
        # ~ else:
        del self.LL[row]
        self.layer_to_grid()
        event.Skip()

    def OnButton_B_Layer_Update(self, event):  # wxGlade: XRR_Frame.<event_handler>
        # ~ print("Event handler 'ev_profile_update' not implemented!")
        self.update = self.V_Fit_Update.GetValue()
        self._pdz   = self.V_Profile_Step.GetValue()
        print(self.atoms)
        self.Chk_P_atom.Clear()
        for a in self.atoms:
            self.Chk_P_atom.Append(a)
        self.Chk_P_atom.SetSelection(0)
        
        self.layer_from_grid()
        self.layer_init()
        self.layer_to_grid()
        self.layer_profile()
        self.xrr_layer()  
        self._plot()
        event.Skip()
        



    def OnButton_B_Fit_Start(self, event):  # wxGlade: XRR_Frame.<event_handler>
        # ~ print("Event handler 'OnButton_B_Fit_Start' not implemented!")
        self.fit_start(event)
        event.Skip()

    def OnButton_B_Fit_Stop(self, event):  # wxGlade: XRR_Frame.<event_handler>
        # ~ print("Event handler 'OnButton_B_Fit_Stop' not implemented!")
        self.fit_stop(event)
        event.Skip()
        
        
        
        
    def OnMenu_Close(self, event):  # wxGlade: XRR_Frame.<event_handler>
        print("Event handler 'OnMenu_Close' not implemented!")
        event.Skip()
        self.Close()
        
    def OnClose(self, evt):
        # ~ self.Bind(wx.EVT_CLOSE, self.OnClose)
        print( "OnClose()")
        if self.DoYouWantTo(" SAVE results "):
            return
        if not self.DoYouWantTo( "continue FITTING"):
            evt.Skip()
            # ~ evt.StopPropagation()
            # ~ self.Close()
            # ~ self.Destroy()

class MainApp(wx.App):
    """Class Main App."""
    def OnInit(self):
        """Init Main App."""
        self.frame = MainFrame(None, wx.ID_ANY, "")
        self.SetTopWindow(self.frame)
        self.frame.Show(True)
        return True
        
        
        
    

if __name__ == '__main__':
    print("XRR WX START")
    app = MainApp(0)
    app.MainLoop()

        

# ~ insert test 
# ~ a = [0,1,2,3,4]
# ~ a.insert(1, 1.5)
# ~ b = a[:2] + [3, 3.1, 3.2] + a[2:]
# ~ c = a[:]
# ~ c[4:4] = [4.1, 4.2, 4.3]
# ~ print(a)
# ~ print(b)
# ~ print(c)
# ~ index = first_list.index(item1)
# ~ del first_list[index]
# ~ del other_list[index]

#---------------------------------------------------------------------------
# ~ class PyEvDat(wx.PyEvent):
        # ~ """Simple event to carry arbitrary result data."""
        # ~ def __init__(self, EVT_ID, data):
            # ~ """Init Result Event."""
            # ~ wx.PyEvent.__init__(self)
            # ~ self.SetEventType(EVT_ID)
            # ~ self.data = data
            
            
#===============================================
# ~ https://stackoverflow.com/questions/26312061/what-is-the-difference-between-wx-lib-newevent-newevent-and-wx-neweventtype
# ~ A wx.lib.newevent.NewEvent() is just an easier wxpython way thats been added to make a wx.NewEventType().
# ~ if you have a look at the code in the module newevent you will see what it does.
# ~ """Easy generation of new events classes and binder objects"""
# ~ __author__ = "Miki Tebeka <miki.tebeka@gmail.com>"

# ~ def NewEvent():
    # ~ """Generate new (Event, Binder) tuple
        # ~ e.g. MooEvent, EVT_MOO = NewEvent()
    # ~ """
    # ~ evttype = wx.NewEventType()

    # ~ class _Event(wx.PyEvent):
        # ~ def __init__(self, **kw):
            # ~ wx.PyEvent.__init__(self)
            # ~ self.SetEventType(evttype)
            # ~ self.__dict__.update(kw)

    # ~ return _Event, wx.PyEventBinder(evttype)
    
# ~ https://discuss.wxpython.org/t/wx-pyevent-losing-my-mind/34462/2
# ~ >>> import wx.lib.newevent as ne
# ~ >>> MooEvent, EVT_MOO = ne.NewEvent()
# ~ >>> 
# ~ >>> evt = MooEvent(data=dict(a=1, b=2, c=3))
# ~ >>> 
# ~ >>> evt.data
# ~ {'a': 1, 'b': 2, 'c': 3}
# ~ >>> 
# ~ key = 0
# ~ d = {0:RUN, -1:ABORT, None:END}
# ~ if key in d:
    # ~ func = d[key]
# ~ else:
    # ~ func = RUN
# ~ func()
# ~ d.get(key, RUN)()   # default = RUN
