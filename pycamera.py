#!/usr/bin/env python
"""
file `pycamera.py`

pycamera: absorption imaging, pixis camera
=============================================

This a the main file of the application. It build a series of objects that it
is going to represent in a window.

The windows is itself an object, it has two panesl, each objects. 

The left hand side panel is divided in two and has to sons: the upper panel and
the lower panel. Each are objects derived from MPLWidget: the MatPlotLib
Widget, wich is a wrapper to a MatPlotLib canvas and allows plotting of
scientific graphs.

The lower left panel is a simple MPLWidget object, and exposes the `axes` of a
MatPlotLib graph. Similar to MatLab axes, these allow to plot in them.

The upper left panel is a `MatrixView` object. It has a "data" attribute,
which should be a dictionary of 2D arrays, and it creates a pseudo-color
representation of these matrices with its "update()" method.

The right panel contains a representation of the different objects used to 
control the camera and to describe the experiment. The representation is
automatically created from the object's attributs using the traitsUI package.
see http://code.enthought.com/traits/ .
"""
##########################################################################
# Imports 
##########################################################################
# To have Traits understand numpy
import os
os.environ['NUMERIX']='numpy'

# Need 'io' from 'scipy' for the 'savemat' function
from scipy.io import savemat

#### Enthought imports. 
# Traits: allows manifest typing of object's attributs and events dispatch 
# while setting attributes.
from enthought.traits import *
# TraitsUI: allows graphical representation of traited objects
from enthought.traits.ui import View, Group, Item, ui, ButtonEditor, \
        InstanceEditor, spring, VGroup, TextEditor
# pyface: the abstraction layer to the GUI windows and panels.
# We have to try/except this import as it has changed in the latest version
# of the enthought library.
try:
    from enthought.pyface.api import SplitApplicationWindow, GUI, SplitPanel
except ImportError:
    from enthought.pyface import SplitApplicationWindow, GUI, SplitPanel

# Display the traits errors on the console:
import enthought.traits.ui.wx.view_application
enthought.traits.ui.wx.view_application.redirect_filename = None

# Use the FBI debugger to catch exceptions:
# In the latest version of the enthought librairy this allows exception to 
# raise the "Frame Based Inspector" debugger.
try:
    # We need to try/except this, as it is only in recent svn.
    from enthought.traits.api import push_exception_handler
    from enthought.debug.fbi import fbi

    push_exception_handler ( handler = lambda o,t,ov,nv: fbi(),
                            main = True,
                            locked = True )
except:
    print "No FBI exception handler."

# Numerical and mathematical tools
from numpy import *

# For labelling saved files
from time import localtime

##########################################################################
# My objects 
##########################################################################
# A thread wrapper that allows a representation of the thread to
# Start/Abort it.
from thread_runner import ThreadRunner
# A software emulation of the camera.
# The camera object
from experiment2 import Camera as PixisCamera
# The plot window with the false color image, and its UI controls.
from matrix_view import MatrixView
# The object containing all the experimental data, and the object performing
# and representing the analysis.
from experiment2 import Experiment, Analysis, Splitting
from data_analysis2 import make_OD

# The object to discribe the camera.
camera = PixisCamera()
# The object to discribe the experiment.
experiment = Experiment()
# The object that deals with details of RF splitting analysis
splitting = Splitting()
# The object that performs the analysis of the data on the fly.
analysis = Analysis(experiment=experiment, camera=camera, splitting=splitting)

#import pdb

##########################################################################
# Open the camera
camera.open()
camera.open_shutter()
# Build the acquisition thread
from time import sleep, localtime
def acquisition_job(thread):
    """ Function that runs the acquisition loop. Will be called in a non
        blocking way (in a seperate thread), and aborts when 
        thread.wants_abort is true.
    """
    # A NOTE ON THE CAMERA SHUTTER:  Shutter opening/closing not controlled
    # manually, but rather by the shutter mode defined in lib_pixis.cpp
    # Shutter opens and remains open during camera trigger.  Shutter closes 
    # some milliseconds after acquisition is complete, and remains closed
    # during readout time (which could be up to 10.5 seconds!)
    
    # If number of kinetics mode sub-images or the binning has changed from 
    # previous acquisition loop, reset the data arrays and axes limits.
    if (camera.num_kin_shots is not camera.old_num_kin_shots) or\
       (camera.binning_X is not camera.old_binning_X) or\
       (camera.binning_Y is not camera.old_binning_Y):
        # Reset max image sizes - mind the binning!
        analysis.matrix_view.y_MAX = camera._kin_win_size / camera.binning_Y
        # NB. vertical binning taken care of in 'camera.num_kin_shots_changed()'
        analysis.matrix_view.x_MAX = camera._CCD_size / camera.binning_X
        # Re-initialize data and axes limits 
        analysis.matrix_view.data_axes_init()
        # NB. The 'old' variables are updated by traits magic in PixisCamera().
    
    # The acquisition loop itself:
    ####################################################################
    while not thread.wants_abort:
        # If the parameters have changed, adjust the buffersize.
#        camera.get_buffer_size()
        data = camera.acquire()
        if thread.wants_abort:
            camera.close_shutter()
            return
        ############################################################
        # Build a dictionary with the different images and construct
        # the absorption image.  Determine number of kinetics mode 
        # sub-images and their limits. Create the divided image 
        ############################################################
        #
        # Integer vertical size of each kinetics sub image (number of rows).
        Ysize = int(data.shape[0]/camera.num_kin_shots) # integer-ize
        #
        # 2 sub-images (OD ~ image2/image2) 
        if camera.num_kin_shots == 2:
            data = {'image1':data[:Ysize,:], 'image2':data[Ysize:,:]}
            data['O.D.'] = make_OD(data['image1'],580,data['image2'],580)
            # define 'image_data_to_save' for disk-space-saving.
            # Save 'signal' and 'reference' images.
            image_data_to_save = data
        #
        # 4 sub-images (OD ~ image2/image3)
        elif camera.num_kin_shots == 4:
            data = {'image1':data[:Ysize,:], 'image2':data[Ysize:2*Ysize,:],\
                 'image3':data[2*Ysize:3*Ysize,:],\
                 'image4':data[3*Ysize:,:]}
            data['O.D.'] = make_OD(data['image2'],data['image4'],data['image3'],data['image4'])
            # define 'image_data_to_save' for disk-space-saving.
            # Save 'signal', 'reference', and 'dark' images.
            image_data_to_save = {'image2':data['image2'],\
               'image3':data['image3'], 'image4':data['image4']}
        #
        ## Do something very bad here: ASSUME Ysize=1024 and hard-code
        ## the sub-image sizes (num. rows) for 5 and 6 subimages only.
        ## Do this because 1024 doesn't divide evenly into 5 or 6.
        #
        # 5 sub-images - use for dual-species Rb-K imaging
        # (OD1 ~ image2/image4, OD2 ~ image3/image5)
        elif camera.num_kin_shots == 5:
            # 205*4 + 204 = 1024
            data = {'image1':data[1:204,:], 'image2':data[204:409,:],\
                 'image3':data[409:614,:],'image4':data[614:819,:],\
                 'image5':data[819:1024,:]}
            data['O.D.1'] = make_OD(data['image2'],580,data['image4'],580)
            data['O.D.2'] = make_OD(data['image3'],580,data['image5'],580)
            # define 'image_data_to_save' for disk-space-saving.
            # Save 'signal', 'reference', and 'dark' images, if applicable.
            image_data_to_save = {'image2':data['image2'],\
               'image3':data['image3'], 'image4':data['image4'],\
               'image5':data['image5']}
        # 6 sub-images - use for dual-species Rb-K imaging
        # (OD1 ~ image2/image4, OD2 ~ image3/image5)
        elif camera.num_kin_shots == 6:
            # 170*5 + 174 = 1024
            data = {'image1':data[1:174,:], 'image2':data[174:344,:],\
                 'image3':data[344:514,:],'image4':data[514:684,:],\
                 'image5':data[684:854,:],'image6':data[854:1024,:]}
            data['O.D.1'] = make_OD(data['image2'],580,data['image4'],580)
            data['O.D.2'] = make_OD(data['image3'],580,data['image5'],580)
            # define 'image_data_to_save' for disk-space-saving.
            # Save 'signal', 'reference', and 'dark' images, if applicable.
            image_data_to_save = {'image2':data['image2'],\
               'image3':data['image3'], 'image4':data['image4'],\
               'image5':data['image5'], 'image5':data['image6']}
        #
        # 8 sub-images 
        # 8 July 2008.  Revert to old 8-shot mode, using just one O.D. image
        # from images 6 and 7.  Resave as 'pycamera-8Jul08.py'
        elif camera.num_kin_shots == 8:
            data = {'image1':data[:Ysize,:], 'image2':data[Ysize:2*Ysize,:],\
             'image3':data[2*Ysize:3*Ysize,:],'image4':data[3*Ysize:4*Ysize,:],\
             'image5':data[4*Ysize:5*Ysize,:],'image6':data[5*Ysize:6*Ysize,:],\
             'image7':data[6*Ysize:7*Ysize,:],'image8':data[7*Ysize:,:]}
            data['O.D.'] = make_OD(data['image6'],data['image8'],data['image7'],data['image8'])
            # define 'image_data_to_save' for disk-space-saving.
            # Save 'signal', 'reference', and 'dark' images.
            image_data_to_save = {'image6':data['image6'], 'image7':data['image7'],\
              'image8':data['image8']}
        #
        # 10 sub-images 
        elif camera.num_kin_shots == 10:
            data = {'image1':data[:Ysize,:], 'image2':data[Ysize:2*Ysize,:],\
              'image3':data[2*Ysize:3*Ysize,:],'image4':data[3*Ysize:4*Ysize,:],\
              'image5':data[4*Ysize:5*Ysize,:],'image6':data[5*Ysize:6*Ysize,:],\
              'image7':data[6*Ysize:7*Ysize,:],'image8':data[7*Ysize:8*Ysize,:],\
              'image9':data[8*Ysize:9*Ysize,:],'image10':data[9*Ysize:,:]}
            data['O.D.'] = make_OD(data['image2'],580,data['image3'],580)
            # define 'image_data_to_save' for disk-space-saving.
            # Save 'signal', 'reference', and 'dark' images.
            image_data_to_save = {'image8':data['image8'], 'image9':data['image9'],\
                'image10':data['image10']}

        else:
            print "*** camera.num_kin_size should be 2,4,8 or 10 ******"
        ######################################################
        # Save data to file, if analysis.save_files == True
        ######################################################
        if analysis.save_files:
            if not analysis.file_name:
                print "**************** No directory specified **************"
                # File dialog opened by analysis object if no file name. 
                # (See experiment.py, "_save_files_changed" function)
            try:
                # Create filename. Use shot numbers with leading zeros.
                # 'shot_number+1'  since it won't have been incremented yet.
                file_name = analysis.file_name +\
                        str("_%04d" % (analysis.shot_number+1) )
                # FIXME: if history buffer runs out, clearing history buffer 
                # and continuing to save will over-write file names (shot number
                # back to zero.  Need to manually update file_name.
                # Not that likely to happen, but still...
               
                # Save some or all kinetics-mode images.
                if analysis.save_all_subimages:
                    # Save ALL kinetics mode subimages.
                    # Overwrite 'image_data_to_save' definitions to "all".
                    image_data_to_save = {} 
                    # Build dic. with image data - don's save absn. data.
                    for key, value in data.iteritems():
                        # save only those images not called "O.D"
                        if (key != 'O.D.') & (key!='O.D.1') & (key!='O.D.2') :
                            image_data_to_save[key] = value
                
                # If analysis.save_al_subimages is False, save only selected
                # kinetics-mode sub-images. The relevant images have already 
                # been selected: after O.D. construction in the 'if' statement 
                # above.  Do nothing!

                # Save data in .mat format (dictionary --> structure)
                # "Orders of magnitude faster than saving in ascii." -G.V.
                # NB. 'savemat' imported from scipy.io
                savemat(file_name, image_data_to_save)
                
                # Record relevant image aquisition and experimental parameters.
                # Array-ize the values to allow clean export with io.savemat
                # Save them in a separate file; include date/time in filename.
                settings_data_to_save = {
                    'ADC_speed_1_2MHz_0_100kHz':array([camera.ADC_speed[\
                        camera.selected_ADC_speed] ], dtype=float32),
                    'row_shift_time_us': 
                        array([camera.row_shift_time],dtype=float32),
                    'gain': 
                        array([camera.gain], dtype=float32),
                    'exposure_time_us': 
                        array([camera.exposure], dtype=float32),
                    'num_kin_shots': 
                        array([camera.num_kin_shots], dtype=float32),
                    'binning_X': 
                        array([camera.binning_X], dtype=float32),
                    'binning_Y': 
                        array([camera.binning_Y], dtype=float32),
                    'CCD_temp_centiC': 
                        array([camera.CCD_temp], dtype=float32),
                    'probe_detuning_MHz':
                        array([experiment.detuning], dtype=float32),
                    'TOF_ms': 
                        array([experiment.tof], dtype=float32),
                    'trap_freq_axial_Hz':
                        array([experiment.trap_freq_ax], dtype=float32),
                    'trap_freq_radial_Hz':
                        array([experiment.trap_freq_rad], dtype=float32),
                 # FIXME: Need to convert these strings to floats
                 #   'view_1_axial_0_radial':
                 #       array([experiment.view_point], dtype=float32),
                 #   'imaging_1_fluor_0_absn':
                 #       array([experiment.imaging], dtype=float32),
                 #   'species_1_K_0_Rb':
                 #       array([experiment.atoms], dtype=float32),
                        }
                # Include time of day in filename. ('join' is a str method.)
                file_name_2 =\
                    file_name+'_settings_' + '_'.join(map(str,localtime()[:-4]))
                savemat(file_name_2, settings_data_to_save)

                print "Saved the data."
            except Exception, inst:
                print "**************** Couldn't save file ******************"
                print inst
        # This starts the analysis ( method analysis._current_data_changed() )
        analysis.current_data = data
        # Update canvas
        try:
            update_canvas(data)
        except ValueError:
            print " can't update canvas"
            pass
        # Do nothing until anal. thread is finished.
        # Important to avoid accessing 'stored_data' too soon.
        # FIXME: THIS SLOWS THE AQUISITION LOOP CONSIDERABLY.
        while analysis.analysis_thread.isAlive(): 1
        # Update stored data (history). Check that history is not 
        # full before updating it
        if analysis.shot_number < analysis._history_len:
            for key, value in analysis.phys_params.iteritems():
                analysis.stored_data[key][analysis.shot_number] = value

# The acquisition thread object
acquisition_thread = ThreadRunner(job = acquisition_job, )


##########################################################################
# The panel class
class MyPanel( HasTraits ):
    """ A wrapper class, to display several objects on the same panel. It 
        allows the graphical representation of several different objects
        in the same panel by using these objects as its attributes.
    """
    
    acquisition_thread = Instance(ThreadRunner)

    experiment = Instance(Experiment)

    camera = Instance(PixisCamera)

    analysis = Instance(Analysis)

    matrix_view = Instance(MatrixView)

    splitting = Instance(Splitting)

    view = View( Group( 
              Group(
                    Item('acquisition_thread',style='custom',
                        show_label=False, ),
                    '_',
                    Item('matrix_view', style='custom', show_label=False),
                    Item('analysis', style='custom', show_label=False, ),
                     label='Acquisition', ),
              Group(Item('experiment', show_label=False, style="custom" ),
                     label='Experiment', ),
              Group(Item('camera', style='custom', show_label=False, ),
                    label='Camera'),
              Group(Item('splitting', style='custom', show_label=False, ),
                    label='Splitting'),
              layout='tabbed',
                ),
            )

# Create the right-hand-side panel object.
my_panel = MyPanel(experiment = experiment, 
            acquisition_thread = acquisition_thread,
            camera = camera,
            analysis = analysis, 
            splitting = splitting,
            )

##########################################################################
# The left hand side split panel
from mplwidget import MPLWidget

class LeftPanel(SplitPanel):
    """ This is a subclass of pyface's SplitPanel, and is therefore a panel
        containing two other widgets, created by self._create_lhs and 
        self._create_rhs.
    """
    
    matrix_view = Instance(MatrixView)
    plot_window = Instance(MPLWidget)

    ratio = Float(0.75)
    direction = Str('horizontal')
    
    def _create_lhs(self, parent):
        """ Upper left-hand display (CCD image). """
        self.matrix_view = MatrixView(parent)
        return self.matrix_view.widget.control

    def _create_rhs(self, parent):
        """ Lower left-hand display (analysis plots). """
        # On the right hand side we use a raw matplotlib widget, this allows
        # more flexibility at the cost of exposing a bit its internals.
        self.plot_window = MPLWidget(parent)
        return self.plot_window.control

##########################################################################
# The application window
##########################################################################
class MainWindow(SplitApplicationWindow):
    """ This is the main window of the program. It is a subclass of pyface's
        SplitApplicationWindows and contains two panels, one that we use to
        represent all the objects of our program, and the other that we use
        to plot graphs.
    """
    left_panel = Instance(LeftPanel)
    traits_panel = Instance(MyPanel)

    ratio = Float(0.75)
    direction = Str('vertical')
    title = Str("PyCamera: Absorption imaging - Pixis Camera")
    size = Tuple(1024, 740)

    def _create_lhs(self, parent):
        self.left_panel = LeftPanel(parent)
        # Now that the plot window exists, add it to the analysis object, 
        # so that it can  plot to it.
        analysis.plot_window = self.left_panel.plot_window
        return self.left_panel.control

    def _create_rhs(self, parent):
        # Link matrix_view object to my_panel (GUI) for cursor interaction. 
        my_panel.matrix_view = self.left_panel.matrix_view
        # Add it to analysis object too, for setting data limits
        analysis.matrix_view = self.left_panel.matrix_view
        self.traits_panel = my_panel
        # Return the control to the traitsUI panel created by the method
        # "edit_traits". That way the panel is embedded in the application
        # window.
        return self.traits_panel.edit_traits(parent = parent,
                                             kind='subpanel',
                                            ).control

    def _on_close(self, event):
        """ Called when the frame is being closed. """
        # Turn the acquisition thread off
        acquisition_thread.abort()
        # FIXME: This is ugly, it seems that the _running boolean of the thread
        # is not properly updated sometimes, so I am limiting the time
        # available for the acquisition thread to die to 3s.
        for i in xrange(1,30):
            if acquisition_thread._running:
                break
            sleep(0.1)
        camera.close()
        self.close()

        return


##########################################################################
# And now, start all this !
if __name__ == '__main__':
    gui = GUI()
    window = MainWindow()
    window.open()
    # Now that the window has been created, we can define the update_canvas
    # function.
    def update_canvas(data):
        window.left_panel.matrix_view.data = data
        window.left_panel.matrix_view.update()

    gui.start_event_loop()
