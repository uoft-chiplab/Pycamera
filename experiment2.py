"""
This file creates the experiment, analysis, and splitting classes.
Why is this called experiment2.py? -- KX 2026
I think this thing uses Python 2.5 or something crazy, common keywords don't work
"""
# To have Traits understand numpy
import os
os.environ['NUMERIX']='numpy'

# Import scipy, the grand unified scientific module !
from scipy import *

# Import figure, plotting commands from plyab for 'analysys window'
#from pylab import figure, show, subplot, imshow

# Enthought imports.
from enthought.traits import *
from enthought.traits.ui import View, Group, Item, ui, ButtonEditor, \
        InstanceEditor, spring, VGroup, TextEditor, EnumEditor
from enthought.traits.ui.menu import \
        OKButton, CancelButton, RevertButton, UndoButton, HelpButton

from enthought.pyface import FileDialog, GUI

file_dialog = FileDialog()

# Import the general purpose threading library
from threading import Thread
from pixis_camera import PixisCamera as Camera
#from mock_camera import PixisCamera as Camera
from mplwidget import MPLWidget
from matrix_view import MatrixView

from time import localtime

# Import analysis functions from data_analysis.py
from data_analysis2 import gaussian1D, dic_to_string, fit_prep,\
        crunch_params, gfit1D, seed_gaussian1D, pretty_print,\
        crunch_split_params, fit_prep2, compute_fit_quality

import pdb

base_dir = '\\\\UNOBTAINIUM\\Carmen_Sandiego\\Data'

##########################################################################
# The experiment object (define BEFORE analysis object!)
##########################################################################
class Experiment( HasTraits ):
    """ A class with all the different experimental parameters. This class
        serves as a representation to the experiment.
    """
    detuning = CFloat(0, 
        desc="detuning of the probe laser to the atomic transition, in Mhz",
        label="Detuning")
    
    tof = CFloat(2.25, 
        desc="time of flight, in ms",
        label = "TOF")

    atoms = Enum("Rb", "K", 
        desc="the species probed",
        label = "Species")

    trap_freq_ax = CFloat(200, 
        desc="the trap frequency in the axial (longitudinal) direction, in Hz",
        label="Axial")

    trap_freq_rad = CFloat(860, 
        desc="the trap frequency in the radial (transverse) direction, in Hz",
        label="Radial")
    
    imaging = Enum("absorption", "fluoescence", 
        desc="the method of imaging the atoms - fluorescence or absorption",
        label="Imaging")

    view_point = Enum("radial", "axial",
        desc="the imaging orientation with respect to atom cloud",
        label="View")
    
    # The view attribute, used by traitsUI to create a representation.
    view = View(Group(
             Group( Item('imaging', width=-120),
                    Item('atoms', width=-60),
                    Item('detuning', width=-40),
                    Item('tof', width=-40),
                   ),
             spring,
             Group( Item('view_point', show_label=True, width=-80),
                    Group(
                        Item('trap_freq_ax', width=-50),
                        spring,
                        Item('trap_freq_rad', width=-50), 
                        orientation='horizontal'),
                    label='Trap frequencies', 
                    show_border = True,),
             ),
                ) 

##########################################################################
# The splitting object (define BEFORE analysis object!)
##########################################################################
class Splitting( HasTraits ):
    """A class with some analysis parameters related to RF splitting
       and atom number counting.
    """
   
    # Define data limits for the left- and right-hand analysis regions
    # for splitting analysis.  These parameters are used in the 'fit_prep'
    # method of data_analysis.py
    left_horiz_left = Int(460, 
        desc="Bottom analysis box, Left-hand Horizontal boundary in pixels",
        label="BHL")
    
    left_horiz_right = Int(575, 
        desc="Bottom analysis box, Right-hand Horizontal boundary in pixels",
        label="BHR")
    
    left_vert_bot = Int(140, 
        desc="Bottom analysis box, Bottom Vertical boundary in pixels",
        label="BVB")
    
    left_vert_top = Int(190, 
        desc="Bottom analysis box, Top Vertical boundary in pixels",
        label="BVT")
    
    right_horiz_left = Int(460, 
        desc="Top analysis box, Left-hand Horizontal boundary in pixels",
        label="THL")
    
    right_horiz_right = Int(575, 
        desc="Top analysis box, Right-hand Horizontal boundary in pixels",
        label="THR")
    
    right_vert_bot = Int(190, 
        desc="Top analysis box, Bottom Vertical boundary in pixels",
        label="TVB")
    
    right_vert_top = Int(240, 
        desc="Top analysis box, Top Vertical boundary in pixels",
        label="TVT")

    # Pixel sum atom number from left and right regions; split fraction
    N_ps_left = CFloat(0.0,
        desc="Bottom analysis box pixel sum atom number",
        label="N_ps_B")
    
    N_ps_right = CFloat(0.0,
        desc="Top analysis box pixel sum atom number",
        label="N_ps_T")

    p_L = CFloat(0.0,
        desc="Fraction of atoms on Bottom", label="p_B")

    L_y = CFloat(0.0,
        desc="Bottom box cloud centre y-coordinate", label="B_y")
    
    L_x = CFloat(0.0,
        desc="Bottom box cloud centre x-coordinate", label="B_x")

    R_y = CFloat(0.0,
        desc="Top box cloud centre y-coordinate", label="T_y")
    
    R_x = CFloat(0.0,
        desc="Top box cloud centre x-coordinate", label="T_x")

    delta_y = CFloat(0.0,
        desc="Difference between the top and bottom center y-coordinates", label="delta_y")

    delta_x = CFloat(0.0,
        desc="Difference between the top and bottom center x-coordinates", label="delta_x")
    
    # The view attribute, used by traitsUI to create a representation.
    view = View(
            Group(
             Group(
                Group(
                    Item('left_horiz_left', width=-40),
                    Item('left_horiz_right', width=-40),
                    Item('left_vert_bot', width=-40),
                    Item('left_vert_top', width=-40),
                     orientation='vertical',
                     ),
                Group(
                    Item('right_horiz_left', width=-40),
                    Item('right_horiz_right', width=-40),
                    Item('right_vert_bot', width=-40),
                    Item('right_vert_top', width=-40),
                     orientation='vertical',
                     ),
                 label='T and B region definitions',
                 orientation='horizontal',
                 show_border = True,
             ),
             Group(
                 Item('N_ps_left', width=-60),
                 Item('N_ps_right', width=-60),
                 Item('p_L', width=-60),
##                 Item('L_y', width=-60),
##                 Item('L_x', width=-60),
##                 Item('R_y', width=-60),
##                 Item('R_x', width=-60),
##                 Item('delta_y', width=-60),
##                 Item('delta_x', width=-60),
                  orientation='horizontal',
                  label='T and B pix sum atom numbers',
                  show_border = True,
             ),
                orientation='vertical',
             ),

             )
##########################################################################
# The analysis object
##########################################################################
class Analysis( HasTraits ):
    """ A class that processes the data in a separate thread and displays 
        the results. This class automatically spawns a new analysis thread 
        (calling its method "do_analysis") when its attribute current_data 
        is changed.
    """

    # This is where you put the images to be processed. For the time being this
    # should be a dictionnary of 2D arrays
    current_data = Any
   
    # Length of stored data array.
    _history_len = Int(2000)
    
    # Physical parameters returned from analysis function.
    # Store them in when analysis is completed. 
    phys_params = Dict
    
    # and units for these physical parameters
    phys_units = {  'N_fit':'',
                    'TxoTf':'',
                    'N_pix_sum':'',
                    'Tx':'kelvin',
                    'Ty':'kelvin',
                    'sigmaX':'pixels',
                    'sigmaY':'pixels',
                    'x_centre':'pixels',
                    'y_centre':'pixels',
                    'psd':'',
                    'N_L':'',
                    'N_R':'',
                    'p_L':'',
                    'L_y':'pixels',
                    'L_x':'pixels',
                    'R_y':'pixels',
                    'R_x':'pixels',
                    'delta_y':'pixels',
                    'delta_x':'pixels',
                    'z':'(p_R - p_L)',
                    'RMSE':'',
                    'goodFit':''}

   

    # Stored_data: use this as a dictionary of arrays to store results from
    # analysis for each run. Also known as the "history".
    # N.B. The stored_data arrays are actually updated by the 'aquisition_job'
    # method in pycamera.py
    stored_data = Any

    def _stored_data_default(self):
        """ Builds the data structure to store the computed data. This is 
            nothing more than a dictionary of arrays.
        """
        return {'N_fit':zeros(self._history_len),\
                'TxoTf':zeros(self._history_len),\
                'N_pix_sum':zeros(self._history_len),\
                'Tx':zeros(self._history_len),\
                'Ty':zeros(self._history_len),\
                'sigmaX':zeros(self._history_len),\
                'sigmaY':zeros(self._history_len),\
                'x_centre':zeros(self._history_len),\
                'y_centre':zeros(self._history_len),\
                'psd':zeros(self._history_len),\
                'N_L':zeros(self._history_len),\
                'N_R':zeros(self._history_len),\
                'p_L':zeros(self._history_len),\
                'L_y':zeros(self._history_len),\
                'L_x':zeros(self._history_len),\
                'R_y':zeros(self._history_len),\
                'R_x':zeros(self._history_len),\
                'delta_y':zeros(self._history_len),\
                'delta_x':zeros(self._history_len),\
                'z':zeros(self._history_len),\
                'RMSE':zeros(self._history_len),\
                'goodFit':zeros(self._history_len),\
                'fit':0}
    # EDIT HERE
    
    selected_stored_data = Str('fit', label="<-- Plot",
        desc="time series data to be plotted in lower left-hand canvas")

    _stored_data_list = List()

    def __stored_data_list_default(self):
        """ Initializes the stored_data enum when the object is created.
        """
        return self.stored_data.keys()
    
    # Manage y-axis scaling of stored (history) data display.
    history_ymin = CFloat(0.0, desc='min. value of plot y-axis',
            label='Plot y-axis min')
    
    history_ymax = CFloat(1.0, desc='max. value of plot y-axis',
            label='max')

    history_autoscale = true(label='Auto',
            desc='if the plot y-axis is autoscaled')

    # Max RMSE of the 1D gaussian fits for a fit to count as "good".
    # Default is permissive; tune against real data.
    rmse_tol = CFloat(1e9, desc='max RMSE for a fit to count as good',
            label='RMSE tol')

    # Event (button) to clear the stored_data buffer when it is full.
    clear_history = Event
        
    # Event (button) to clear the stored_data buffer when it is full.
    clear_history = Event

    # Event (button) to initiate re-analysis of current data.
    refit = Event
    
    # The index of the run in the sequence. It is incremented at each run, and
    # is put to zero when the data is saved.
    shot_number = -1 
    
    # A pointer to the object where the experimental parameters are stored.
    experiment = Instance(Experiment)

    # A pointer to the object where the camera parameters are stored.
    camera = Instance(Camera)

    # A pointer to the object containing splitting analysis info is stored.
    splitting = Instance(Splitting)

    # Thread object that is going to be run to do the analysis.
    analysis_thread = Instance(Thread)
    
    # A pointer to MPLWidget to do some plotting (lower-left axes).
    plot_window = Instance(MPLWidget)
    
    # A pointer to MatrixView for image display (upper-left axes).
    matrix_view = Instance(MatrixView)

    # The text displayed on screen.
    displayed_text = String(desc="the results of the analysis")

    # The file to which save the data.
    file_name = File(desc="the file that the data should be saved to.")
    # 
    # For keeping track of "old" filename whenever file_name is changed.
    #old_file_name = File()
    
    # Boolean type on GUI to indicate whether files should be saved or not.
    save_files = false(desc="whether the images are going to be saved to files",
                        label="Save images")

    # Boolean type on GUI to indicate whether to save ALL kinetics mode
    # images, or only those relevant to data.
    save_all_subimages = true(desc="whetherALL kinetics images should be saved",
                        label="ALL subimgs")
    # Boolean type on GUI to indicate whether save path should be generated
    # automatically
    save_autopath = true(desc="whether to auto-generate the savepath", label="Autopath")
    
    # Event that triggers the save history. This is represented as a button.
    save_history = Event
    
    # user-specified param for cost writing
    cost_string = Str('', desc="fit parameter to be used as cost function for ML optimization.", label="")
    
    # The view attribute, used by traitsUI to create a representation.
    view = View(
            Item('displayed_text',show_label=False, springy=False, 
                style='custom', editor=TextEditor(),),
            Group(
                'save_files',
                'save_all_subimages',
                'save_autopath',
                Item('save_history', editor=ButtonEditor(), show_label=False,
                    width=-70),
                Item('clear_history', editor=ButtonEditor(), show_label=False,
                    visible_when='shot_number >= _history_len - 1',),
                Item('cost_string', show_label=False, springy=False, width=-30),
                orientation = 'horizontal',
            ),
            Item('file_name', show_label=False, springy=False, ),
            '_',
            Group(
                Item('history_ymin', show_label=True, width=-40), spring,
                Item('history_ymax', show_label=True, width=-40), spring,
                Item('history_autoscale', show_label=True),
                orientation = 'horizontal',
            ),
            Group(
                Item('rmse_tol', show_label=True, width=-60),
                orientation = 'horizontal',
            ),
            Group(
                Item('selected_stored_data', 
                    editor=EnumEditor(name='_stored_data_list')),
                spring,
                Item('refit', editor=ButtonEditor(), show_label=False ),
                orientation = 'horizontal',
            )
           )  

    def _save_history_fired(self):
        """ Saves the data.
        """
        data_to_save = {}
        # Build a dictionary with the stored data truncated to the last shot
        for key, value in self.stored_data.iteritems():
            if key != 'fit':
                # go to "shot_number+1" to accomodate python array indexing
                data_to_save[key] = value[:self.shot_number+1].astype(float32)
        # Check that there is a valid filename
        while not self.file_name:

            file_dialog.open()
            self.file_name = file_dialog.path
        # Set the file name for saving 
        file_name = \
            self.file_name + '_history_' + '_'.join(map(str,localtime()[:-4])) 
        # Save the data in .mat format (dictionary --> structure)
        io.savemat(file_name, data_to_save)
        self.display_line('\n History data saved: %d shots.' %(self.shot_number+1))

    def _clear_history_fired(self):
        """ Clear the stored_data buffer
        """
        # Reset stored data arrays to zeros
        self._stored_data_default()
        # Re-start shot counter
        self.display_line("History cleared.")
        self.shot_number = -1

    def _refit_fired(self):
        """ Executed by traits magic when redo_analysis button pressed.
            Force fresh analysis.
        """
        # FIXME: much clicking of this button seems to crash pycamera!
        self.do_analysis()

    def _analysis_thread_default(self):
        """ Creates an empty Thread object
        """
        return Thread()

    def _save_files_changed(self, new_save_files):
        """ Check that the file_name string is not empty
        """
        # Keep track of changes in file name.  Will be used to automatically
        # save "history" if filename is changed.
        old_file_name = self.file_name

        # NB. if new_save_files == False, do nothing.
        # aquisition_job in pycamera.py will take care to stop saving.
        if new_save_files:
            # Check for valid (non-empty) file name.
            
            if self.save_autopath:
                #i.e \\UNOBTAINIUM\Carmen_Sandiego\Data\2026
                path = os.path.join(base_dir, str(localtime().tm_year))
                #because we make folders by date and letter, we can recursively look for the last subdir
                # WARNING THIS WILL BREAK IF AN UNEXPECTED FOLDER APPEARS
                for i in range(3):
                    dirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
                    path = os.path.join(path, dirs[-1])
                    if i==2:
                        fn = dirs[-1]
                path = os.path.join(path, 'imgs')
                self.file_name = os.path.join(path, fn)
                
            if not self.file_name:
                # Set save_files flag to False temporarily while the file name
                # and path get sorted out.  Don't want to save any files without
                # a proper file name.
                self.save_files = False
                file_dialog.open()
                self.file_name = file_dialog.path
                # Check that valid path has been entered.
                if not self.file_name:
                    # Still no valid file name.
                    print "You haven't selected a correct path"
                    # At this point, just resume event loop execution.
                    # Don't save until a valid path is obtained.
                else:
                    self.save_files = True
                    self._clear_history_fired()
                    # Clearing the display to give a visual hint that
                    # only the new data will be saved.
                    GUI.invoke_later(setattr, self, 'displayed_text', '')
                    print "Starting to save data now. Thanks for filename."
            else:
                # Valid path - start saving.
                self._clear_history_fired()
                # Clearing the display to give a visual hint that the 
                # data will not be saved.
                GUI.invoke_later(setattr, self, 'displayed_text', '')
                print "Starting to save data now."

    def _current_data_changed(self, new_data):
        """ Called automaticaly by the Traits magic when the 
            current_data attribute is changed. Starts the analysis in a new 
            thread.
        """
        # Increment the shot number.
        self.shot_number += 1
        #print "after incr., shot number = %d" % self.shot_number
        if not self.analysis_thread.isAlive():
            # We have to create a new thread object, and start it:
            self.analysis_thread = Thread()
            self.analysis_thread.run = self.do_analysis
            self.analysis_thread.start()
        else:
            print "+++ The analysis takes too long, it skipped a beat"
    
    def do_analysis(self):
        """ This is the function that is called to do the analysis. It should 
            perform the analysis using the parameters stored in the camera and
            experiment objects, and the experimental results, stored in the
            current_data attribute of this object.
        """
        #pdb.set_trace()
        # Fit data, analyze, spit out results and plot.
        # Retrive experimental parameters experiment, matrix_view, and
        # camera pointers.
        ######################################################################
        # Manually update data limits before performing analysis.
        # See update() method in MatrixView class, for instance.
        self.matrix_view.x_min, self.matrix_view.x_max =\
                self.matrix_view.axes.get_xlim()
        self.matrix_view.y_min, self.matrix_view.y_max =\
                self.matrix_view.axes.get_ylim()
        # If user has 'dragged' axes out of the image bounds
        if (self.matrix_view.x_min < 0) | (self.matrix_view.y_min < 0)\
            | (self.matrix_view.y_max > self.matrix_view.y_MAX)\
            | (self.matrix_view.x_max > self.matrix_view.x_MAX):
            print "You dragged image out of bounds..."
            return
        #pdb.set_trace()
        # Set up data for 1D fitting and pixel sum counting.
        # This is the original fit prep 
        h_xdata, h_gdata, v_xdata, v_gdata, pix_sum, bg,\
            pix_sum_left, pix_sum_right = \
                fit_prep(self.current_data, self.matrix_view,\
                    self.splitting)

        #prep data for 1D fits for both left and right boxes

        L_h_xdata, R_h_xdata, L_v_xdata, R_v_xdata, L_h_gdata,\
           L_v_gdata, R_h_gdata, R_v_gdata = fit_prep2(self.current_data, self.matrix_view,\
                    self.splitting)

        # Seed the 1D fits.
        h_gseed, v_gseed = seed_gaussian1D(h_xdata, h_gdata, v_xdata, v_gdata,
                                            self.matrix_view)
        # Seed 1D fits for left and right boxes
        #May require matrix argument change to sp.

        L_h_gseed, L_v_gseed = seed_gaussian1D(L_h_xdata, L_h_gdata, L_v_xdata, L_v_gdata,\
                                            self.matrix_view)
        R_h_gseed, R_v_gseed = seed_gaussian1D(R_h_xdata, R_h_gdata, R_v_xdata, R_v_gdata,\
                                            self.matrix_view)
        
        # Perform the fits.
        h_fit = gfit1D(h_xdata, h_gdata, h_gseed)
        v_fit = gfit1D(v_xdata, v_gdata, v_gseed)

        L_h_fit = gfit1D(L_h_xdata, L_h_gdata, L_h_gseed)
        L_v_fit = gfit1D(L_v_xdata, L_v_gdata, L_v_gseed)

        R_h_fit = gfit1D(R_h_xdata, R_h_gdata, R_h_gseed)
        R_v_fit = gfit1D(R_v_xdata, R_v_gdata, R_v_gseed)

        # Compute the relevant physical parameters.
        self.phys_params = crunch_params(self.experiment, self.camera,\
            self.matrix_view, h_fit, v_fit, L_h_fit, L_v_fit, R_h_fit, R_v_fit,\
                                         pix_sum, pix_sum_left, pix_sum_right, self._stored_data_list)
        # Compute atom number in left, right regions; fraction too
        self.splitting.N_ps_left, self.splitting.N_ps_right, \
            self.splitting.p_L = \
            crunch_split_params(self.experiment, self.camera,\
                pix_sum_left, pix_sum_right)
        # Goodness-of-fit of the 1D gaussian fits against the profile data.
        rmse, goodFit = compute_fit_quality(h_xdata, h_gdata, h_fit,\
            v_xdata, v_gdata, v_fit, self.rmse_tol)
        self.phys_params['RMSE'] = rmse
        self.phys_params['goodFit'] = goodFit
        # Display the params.
        self.display_line(self.make_output_string())
        # Write cost from specified param if exists
        if self.cost_string in self.phys_params and self.save_files and self.file_name:
            self.write_cost(self.cost_string)
        # Clear lower left and canvas for plotting
        self.plot_window.figure.clf()
        # Plot one of the stored_data dic. items OR fits to the cloud.
        if self.selected_stored_data == 'fit':
            # Left-hand plot
            plot_axes1 = self.plot_window.figure.add_subplot(1,2,1) 
            # clear these axes
            plot_axes1.cla()
            # plot horiz. secition and fit of gaussian data 
            plot_axes1.plot(h_xdata, h_gdata, "b-")
            plot_axes1.plot(h_xdata, gaussian1D(h_xdata, h_fit), "r")
            plot_axes1.set_title('horizontal [X]')
            plot_axes1.grid()
            # Right hand plot
            plot_axes2 = self.plot_window.figure.add_subplot(1,2,2) 
            plot_axes2.cla()
            plot_axes2.plot(v_xdata, v_gdata, "g-")
            plot_axes2.plot(v_xdata, gaussian1D(v_xdata, v_fit), "r")
            plot_axes2.set_title('vertical [Y]')
            plot_axes2.grid()
        else:
            # Define some axes and set them to be 'current axes'
            self.plot_window.figure.add_axes()
            plot_axes = self.plot_window.figure.gca()
            # Get history data 
            data_to_plot = self.stored_data[self.selected_stored_data]\
                [:self.shot_number+1]
            # Plot the last 50 shots in red
            image_nums = arange(0, self.shot_number+1)
            plot_axes.plot(image_nums[-50:], data_to_plot[-50:], 'r',\
                linewidth=3)
            # User-defined y-axis scaling
            if self.history_autoscale == False:
                plot_axes.set_ylim((self.history_ymin, self.history_ymax))
            else:
                # force axis autoscaling
                plot_axes.set_ylim((\
                    data_to_plot[-50:].min(),data_to_plot[-50:].max() ))
            # Label plot, add grid
            plot_title = self.selected_stored_data + " vs. image number"
            plot_axes.set_title(plot_title)
            plot_axes.grid()
        # And update the display
        GUI.invoke_later(self.plot_window.figure.canvas.draw)
        print "Analysis done"
        
    def display_line(self, string):
        """ Thread safe display event, this adds a line to the displayed text
            in a thread safe way. Always call this function: do not modify
            the displayed sting yourself.
        """
        new_text = "%s\n%s\n\n" % (string, self.displayed_text)
        # Limit the length of the report, to prevent crashes
        new_text = new_text[:3000]
        GUI.invoke_later(setattr, self, 'displayed_text', new_text)
        
    def make_output_string(self):
        """ Manually construct a formatted string of analysis for output."""
        output_string = \
                "\n*****   Image no. %d   ***** \n\
                Atom number (fit):  %s\n\
                Atom number (pixel sum):  %s\n\
                Tx:  %s kelvin\n\
                Ty:  %s kelvin\n\
                Tx/TF:  %s\n\
                sigmaX:  %s pixels\n\
                sigmaY:  %s pixels\n\
                Xc:  %s pixels\n\
                Yc:  %s pixels\n\
                PSD:  %s\n\
                N_B:  %s\n\
                N_T:  %s\n\
                p_B:  %s\n\
                B_y:  %s\n\
                B_x:  %s\n\
                T_y:  %s\n\
                T_x:  %s\n\
                delta_y:  %s\n\
                delta_x:  %s\n\
                z:  %s\n\
                RMSE:  %s\n\
                goodFit:  %s" % (\
                self.shot_number,\
                pretty_print(self.phys_params['N_fit']),\
                pretty_print(self.phys_params['N_pix_sum']),\
                pretty_print(self.phys_params['Tx']),\
                pretty_print(self.phys_params['Ty']),\
                pretty_print(self.phys_params['TxoTf']),\
                pretty_print(self.phys_params['sigmaX']),\
                pretty_print(self.phys_params['sigmaY']),\
                pretty_print(self.phys_params['x_centre']),\
                pretty_print(self.phys_params['y_centre']),\
                pretty_print(self.phys_params['psd']),\
                pretty_print(self.phys_params['N_L']),\
                pretty_print(self.phys_params['N_R']),\
                pretty_print(self.phys_params['p_L']),\
                pretty_print(self.phys_params['L_y']),\
                pretty_print(self.phys_params['L_x']),\
                pretty_print(self.phys_params['R_y']),\
                pretty_print(self.phys_params['R_x']),\
                pretty_print(self.phys_params['delta_y']),\
                pretty_print(self.phys_params['delta_x']),\
                pretty_print(self.phys_params['z']),\
                pretty_print(self.phys_params['RMSE']),\
                pretty_print(self.phys_params['goodFit']) )
        return output_string
        
    def write_cost(self, param_str):
        """
            For ML optimization. After a fit, write selected parameter to cost_ddd.txt
            Only write the real value when the fit is good; otherwise write 0.
        """
        cost = self.phys_params[param_str]
        if not self.phys_params.get('goodFit'):
            cost = 0
        index = self.shot_number
        filename = "cost_%03d.txt" % index
        filepath = os.path.dirname(self.file_name)
        writepath = os.path.join(filepath, filename)
        wp = open(writepath, 'w')
        try:
            wp.write(str(cost))
        finally:
            wp.close()
