#!/usr/bin/env python
"""
This file defines a pyface widget class to display a matrix in a pseudo color
plot.
"""
import pdb
from numpy import zeros 

# My imports
from mplwidget import MPLWidget

import matplotlib

from enthought.traits import *
from enthought.traits.ui import View, Item, Group, EnumEditor,\
        ButtonEditor, spring
try:
    from enthought.pyface.api import GUI
except ImportError:
    from enthought.pyface import GUI

##########################################################################
# A cursor object 
##########################################################################
class StaticCursor(object):
    """ A cross object that is displayed in axes, severely copied from 
        matplotlib.widgets.Cursor
    """
    # A boolean indicating whether cursor is visible or not.
    # Added because self.visible never seems to be False, and I couldn't
    # figure out which attribute of 'axes' would tell me this.
    is_visible = false

    def __init__(self, ax, useblit=False, **lineprops):
        """
        Add a cursor to ax.  If useblit=True, use the backend
        dependent blitting features for faster updates (GTKAgg only
        now).  lineprops is a dictionary of line properties.  See
        examples/widgets/cursor.py.
        """
        self.ax = ax
        self.canvas = ax.figure.canvas

        self.visible = True
        self.horizOn = True
        self.vertOn = True
        self.useblit = useblit

        self.lineh = ax.axhline(0, visible=False, **lineprops)
        self.linev = ax.axvline(0, visible=False, **lineprops)

        self.background = None
        self.needclear = False

    def clear(self, event):
        'clear the cursor'
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.linev.set_visible(False)
        self.lineh.set_visible(False)
        self.is_visible = False # a cheap hack
        self._update()

    def draw(self, x, y):
        'draw the cursor if visible'
        if not self.visible: return
        self.linev.set_xdata((x, x))

        self.lineh.set_ydata((y, y))
        self.linev.set_visible(self.visible and self.vertOn)
        self.lineh.set_visible(self.visible and self.horizOn)
        self.is_visible = True # a cheap hack
        self._update()

    def _update(self):
        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            self.ax.draw_artist(self.linev)
            self.ax.draw_artist(self.lineh)
            self.canvas.blit(self.ax.bbox)
        else:
            self.canvas.draw_idle()

        return False

##########################################################################
# Traited imshow view
##########################################################################
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from matplotlib.cm import jet, gray, gray_r, copper, gist_earth, pink

class MatrixView( HasTraits ):
    """ A MatPlotLib Widget displaying a pseudo color view of a 
        2D numpy array.
    """
    # A dictionary of 2D arrays to be plotted.
    data = Dict

    # Attributs starting by "_" are internal to the object and should not be
    # accessed. 
    _image_list = List(['image1',],)

    selected_image = Str('image1', desc="the image shown on the screen",
                                   label='Display')

    # Bolean that decides whether the colormap should be autoscaled or not.
    # 'true' is a trait!
    autoscale_colormap = true(label = "Auto",  
                              desc="if the colormap is autoscaled")

    # Event that is represented as a button.
    rescale = Event(label = "Rescale")
    
    # The max for the colormap
    z_max = CFloat

    # The min for the colormap
    z_min = CFloat
  
    def get_z_extrema(self):
        """Return max and min values of the currently selected image 
           for use in colour mapping. Use current data limits!
           Return min, max.
        """
        return 0.001*int(self.data[self.selected_image]\
            [self.y_min:self.y_max, self.x_min:self.x_max].min()*1000), \
            0.001*int(self.data[self.selected_image]\
            [self.y_min:self.y_max, self.x_min:self.x_max].max()*1000)
    # Create GUI "doubles" of z_max and z_min so that scaling only 
    # occurs when "Rescale" buton pressed, rather than as soon as new 
    # values are typed in. 
    z_max_gui = CFloat(35000.0, desc="max. value of colour map",
                             label = "max")
    z_min_gui = CFloat(600.0, desc="min. value of colour map",
                                label = "min")

    colormaps =  {   
                    'gist_earth': gist_earth,
                    'jet': jet,
                    'grey':  gray,
                    'grey_r': gray_r,
                    'copper': copper,
					'pink': pink
                 }

    selected_colormap = Str('gist_earth',
                        label='Cmap',
                        desc='colormap used to display the image')

    _colormap_list = List()

    def __colormap_list_default(self):
        """ Initializes the colormap enum when the object is created """
        return self.colormaps.keys()
    
    # Enum of data-subset-making tools: cursor, or region-of-interest (ROI) sum
    data_subset = Enum("ROI sum", "cursor",
            desc="whether to fit to ROI summed data, or slices at the\
                cursor position",
            label="data subset")
    # NB. Cursor positions run from 1 to the number of rows/columns in the
    # data array. They do not start counting from zero.
    # Horizontal (x cursor position
    cursor_x_pos = CFloat(0.0, desc="cursor x position on image canvas",
                            label="x")
    
    # Vertical (y) position
    cursor_y_pos = CFloat(0.0, desc="cursor y position on image canvas",
                            label="y")
    
    # Data value at cursor
    cursor_data = CFloat(0.0, desc="value of image array at cursor",
                             label="value")

    # Data limtis for image display (MatrixView) and analysis
    x_min = CFloat(desc="lower x limit on image data", label="xmin")
    x_max = CFloat(desc="upper x limit on image data", label="xmax")
    y_min = CFloat(desc="lower y limit on image data", label="ymin")
    y_max = CFloat(desc="upper y limit on image data", label="ymax")

    # Maximum image data size based on kinetics mode window size.
    # These attributes are touched in the acquisition job in pycamera.py
    # (Minimum index values are zero.)
    x_MAX = Int(1024)
    y_MAX = Int(512)

    # Create GUI objects using 'View'
    view = View(Group(
            Group( Item('selected_image', 
                            editor=EnumEditor(name='_image_list'), 
                            show_label=False),
                    'autoscale_colormap',
                    Item('selected_colormap', 
                            editor=EnumEditor(name='_colormap_list')),
                    orientation='horizontal' ),
            Group(  Item('z_min_gui', show_label=True, width=-50),
                            #style='readonly'), 
                    Item('z_max_gui', show_label=True, width=-50),
                            #style='readonly'), 
                    Item('rescale', editor=ButtonEditor(), show_label=False, 
                            width=-45),
                    Item('data_subset', show_label=False),
                     orientation='horizontal' ),
            Group(  Item('cursor_x_pos', show_label=True, width=-40, 
                            style='readonly'), 
                    spring,
                    Item('cursor_y_pos', show_label=True, width=-40, 
                            style='readonly'), 
                    spring, 
                    Item('cursor_data', show_label=True, width=-60, 
                            style='readonly'), 
                    spring,
                       orientation='horizontal', show_border=True, 
                       label="Cursor"),
                   ))
    
    def _selected_image_changed(self):
        """ Executed when user selects new image to display.
        """
        self.update()

    def _selected_colormap_changed(self):
        """ Executed when user selects new colormap.
        """
        self.update()

    def _rescale_fired(self):
        """Executed when 'Rescale' button clicked by user. """
        #self.autoscale_colormap = False
        self.update()

    def do_rescale_colormap(self):
        """ Recompute limits for colour scaling.
        """
        if self.autoscale_colormap:
            # Set scale based on max and min in image.
            self.z_min, self.z_max = self.get_z_extrema()
            self.z_min_gui = self.z_min
            self.z_max_gui = self.z_max
        else:
            # Set scale based on gui values.  
            # Check that max >= min
            if self.z_max_gui >= self.z_min_gui:
                self.z_max = self.z_max_gui 
                self.z_min = self.z_min_gui
            else:
                # Revert to autoscale if max < min.
                # Same as case when self.autoscale_colormap == True
                self.autoscale_colormap=True
                self.z_min, self.z_max = self.get_z_extrema()
                self.z_min_gui = self.z_min
                self.z_max_gui = self.z_max
  
    def __init__(self, parent):
        """Executed upon instantiation of MatrixView class """
        matplotlib.rcParams.update({'xtick.major.pad':-20,
                                    'ytick.major.pad':-35,
                                    'xtick.color': 'w',
                                    'ytick.color': 'w',
                                     })
        self.widget = MPLWidget(parent)
        self.widget.figure.clear()
        self.axes = self.widget.figure.add_axes([0., 0., 1., 1.])
        self.canvas = self.widget.figure.canvas
        # Initialize cursor object (private)
        self._cursor = StaticCursor(self.axes, color='red', linewidth=2,
                                            linestyle = "--")
        self.canvas.mpl_connect('button_press_event',
                            self.mouse_button_callback)
        # Initialize the data and axes limits
        self.data_axes_init()
        # Update canvas 
        self.update()
        matplotlib.rcdefaults()
   
    def data_axes_init(self):
        """ Initialize data array and data/axes limits.  Called by __init__
            and in pycamera.py ."""
        # Initialize the data - y_MAX and x_MAX set by pycamera.py to include
        # the proper binning factors.
        self.data={ 'image1':zeros((self.y_MAX, self.x_MAX), dtype='i') }
        # Initialize the data limits and axes limits
        self.y_max, self.x_max = (self.y_MAX, self.x_MAX) 
        self.y_min, self.x_min = (0.0, 0.0)
        self.axes.set_xlim((self.x_min, self.x_max))
        self.axes.set_ylim((self.y_min, self.y_max))
   
    def mouse_button_callback(self, event):
        """ Fired by mpl when there is a mouse click on the canvas.
        """
        if not event.inaxes:
            return
        if event.button == 3 : 
            # Right-click - draw cursor lines 
            self._cursor_position = (event.xdata, event.ydata)
            self._cursor.draw(event.xdata, event.ydata)
            # Get cursor position and value of image array at that position.
            self.cursor_x_pos = int(self._cursor_position[0])
            self.cursor_y_pos = int(self._cursor_position[1])
            self.cursor_data = 0.001*int(self.data[self.selected_image]\
                [self.cursor_y_pos, self.cursor_x_pos]*1000)
        elif event.button == 1: 
            # Left-click - clear the cursor
            self._cursor.clear(event)
        
    def update(self):
        """ Updates the display in a clean way (event-loop friendly)
        """
        self.axes.images=[]
        if self.selected_image not in self.data:
            self.selected_image = self.data.keys()[0]
        # Get data limits from current axes limits 
        (self.x_min, self.x_max) = self.axes.get_xlim()
        (self.y_min, self.y_max) = self.axes.get_ylim()
        # Re-calculate color scaling limits before displaying image.
        self.do_rescale_colormap()
        # Display image (nicely formatted, a matplotlib method).
        self.axes.imshow(self.data[self.selected_image],
                            cmap=self.colormaps[self.selected_colormap],
                            vmin=self.z_min, vmax=self.z_max,
                            aspect='auto',
                            interpolation='nearest',
                            origin='lower') # Key to vert. coordinates!!!!
        # Set data limits for refreshed image 
        self.axes.set_xlim((self.x_min, self.x_max))
        self.axes.set_ylim((self.y_min, self.y_max))
        # Update value under cursor, if cursor is visible.
        if self._cursor.is_visible: 
            try:
                self.cursor_data = 0.001*int(self.data[self.selected_image]\
                    [self.cursor_y_pos, self.cursor_x_pos]*1000)
            except Exception, inst:
                print "******** Cursor out of bounds.  Try zooming out *****"
                print inst
        # Draw and update canvas.
        GUI.invoke_later(self.canvas.draw) 
        GUI.invoke_later(self.update_list)

    def update_list(self):
        """ Updates image list. Must be called from the event-loop !
        """
        if not self._image_list == self.data.keys():
            self._image_list = self.data.keys()
