#!/usr/bin/env python
"""
A Pyface widget to embed a matplotlib plot widget in a pyface form.

"""
# Author: Prabhu Ramachandran <prabhu@aero.iitb.ac.in>
# Copyright (c) 2006, Prabhu Ramachandran
# License: BSD Style.

import wx

import matplotlib
# Uncomment/comment appropriately to use different backends.
#matplotlib.use('WX')
#from matplotlib.backends.backend_wx import FigureCanvasWx as FigureCanvas
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas

from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

# Enthought imports.
from enthought.traits import Any, Instance, Float, Str
from enthought.pyface import Widget

import pdb
##########################################################################
# The Matplotlib PyFace widget.
##########################################################################
class MPLWidget(Widget):
    """ A MatPlotLib PyFace Widget """
    # Public traits
    figure = Instance(Figure)
    axes = Instance('matplotlib.axes.Axes')
    mpl_control = Instance(FigureCanvas)

    # Private traits.
    _panel = Any
    _sizer = Any
    _toolbar = Any
    
    def __init__(self, parent, **traits):
        """ Creates a new matplotlib widget. """
        # Calls the init function of the parent class.
        super(MPLWidget, self).__init__(**traits)

        self.control = self._create_control(parent)


    def _create_control(self, parent):
        """ Create the toolkit-specific control that represents the widget. """
        #pdb.set_trace()
        # The panel lets us add additional controls.
        if isinstance(parent,wx.Panel):
            self._panel = parent
        else:
            self._panel = wx.Panel(parent, -1, style=wx.CLIP_CHILDREN)
        self._sizer = wx.BoxSizer(wx.VERTICAL)
        self._panel.SetSizer(self._sizer)
        # matplotlib commands to create a figure, and add an axes object
        self.figure = Figure()
        self.figure.set_facecolor((0.9296875, 0.91015625, 0.91015625))
        # Yes we want to allow small figures
        self.axes = self.figure.add_axes([0.05, 0.04, 0.9, 0.92])
        self.mpl_control = FigureCanvas(self._panel, -1, self.figure)
        self._sizer.Add(self.mpl_control, 1, wx.LEFT | wx.TOP | wx.GROW)
        self._toolbar = self._get_toolbar()
        self._sizer.Add(self._toolbar, 0, wx.EXPAND)
        self._sizer.Layout()
        self.figure.canvas.SetMinSize((100,100))
        return self._panel
    

    def _get_toolbar(self):
        return NavigationToolbar2Wx(self.mpl_control)


##########################################################################
# A sample application.
##########################################################################
try:
    from enthought.pyface.api import SplitApplicationWindow, PythonShell, GUI
except ImportError:
    from enthought.pyface import SplitApplicationWindow, PythonShell, GUI

class MainWindow(SplitApplicationWindow):

    mpl = Instance(MPLWidget)
    shell = Instance(PythonShell)

    ratio = Float(0.75)
    direction = Str('horizontal')

    def _create_lhs(self, parent):
        self.mpl = MPLWidget(parent)
        return self.mpl.control

    def _create_rhs(self, parent):
        shell = self.shell = PythonShell(parent)
        shell.bind('mpl', self.mpl)
        return self.shell.control

if __name__ == '__main__':
    gui = GUI()
    window = MainWindow()
    window.open()
    window.size = (800, 700)
    # Simple test case.
    from pylab import arange, sin
    x = arange(1,10,0.1)
    axes = window.mpl.axes.plot(x,sin(x))
    gui.start_event_loop()

