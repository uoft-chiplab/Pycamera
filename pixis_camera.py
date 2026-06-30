"""
This file creates a camera class.

This implementation is for the Pixis camera.
"""

# To have Traits understand numpy
import os
os.environ['NUMERIX']='numpy'

from numpy import random, zeros, array

# Enthought imports.
from enthought.traits import *
from enthought.traits.ui import View, Group, Item, ui, ButtonEditor, \
        InstanceEditor, spring, VGroup, EnumEditor

# Import the low level interface to C.
import pixis_interface

from time import sleep

##########################################################################
# The camera object
##########################################################################

class PixisCamera( HasTraits ):
    """ Class that implements a Pixis camera. The attributs of the class 
        describe the camera and are used to create a graphical representation.

        High level methods are provided to control the camera.
    """
    ADC_speed = { '100kHz':0, '2MHz':1 }

    selected_ADC_speed = Str('100kHz', label='ADC',
                desc='the speed of ADC digitiziation / readout')
    
    _ADC_list = List()

    def __ADC_list_default(self):
        """ Initialized the ADC speed enum when the object is created.
        """
        return self.ADC_speed.keys()
    
    row_shift_time = Enum(3.2, 5.0, 15.0, 30.2, label="Shift",
            desc="the time in microseconds to shift one row of pixels\
            in kinetics mode")

    gain = Enum(1,2,3,
        desc="the gain index of the camera",
        label="gain", )

    exposure = CInt(100,
        desc="the exposure time, in microseconds",
        label="t_exp", )
    
    number_images = Int(1,
        desc="the number of images acquired in one sequence",
        label="Number of images", ) 
        # This is not shown as I do not want people to add
        # more Images and screw the acquisition process.
    
    _CCD_size = Int(1024)   # The Pixis 1024BR has a 1024 x 1024 pixel array
   
    # Kinetics mode acquisition settings.
    _kin_win_size = Int(512,
        desc="the kinetics windows size, in the parallel direction",
        label="Kinetics window", )
        # This is not shown.  Should be calculated based on num_kin_shots and
        # CCD_size only!  Default value = 1024/2 = 512
    
    num_kin_shots = Enum(2,4,5,6,8,10,
    #3, avoid 3-shot kinetics - nasty arithmatic!
    # 5 and 6 added Apr.08 for dua-species imaging. Nasty arithmatic!
        desc='the number of kinetics mode sub-images',
        label='Sub-imgs.')
    
    binning_X = Enum(1,2,4,8,
        desc="Amount of binning on the horizontal direction",
        label="Bins:  horiz.")
    
    binning_Y = Enum(1,2,4,8,
        desc="Amount of binning on the vertical direction",
        label="vert.")
    
    # Keep track of 'previous' kinetics mode settings.
    old_num_kin_shots = Int(2)
    old_binning_X = Int(1)
    old_binning_Y = Int(1)

    def _num_kin_shots_changed(self, old, new):
        """ Executed by traits magic wth num_kin_shots changes.
            Updates kin_win_size.
        """
        # Divide, round, then truncate to integer.
        # NB. Binning is taken care of in acquisition job (pycamera.py)
        self._kin_win_size = ( round(\
            (self._CCD_size/self.num_kin_shots), 0) ).__int__()
        # Update 'old' attribute
        self.old_num_kin_shots = old
        #print "old_num_kin_shots = %d" % self.old_num_kin_shots
        
    def _binning_X_changed(self, old, new):
        """ Executed by traits magic when binning_X changes."""
        self.old_binning_X = old
        #print "old_binning_X = %d" % self.old_binning_X
    
    def _binning_Y_changed(self, old, new):
        """ Executed by traits magic when binning_Y changes."""
        self.old_binning_Y = old
        #print "old_binning_Y = %d" % self.old_binning_Y
    
    # Use this to hide num_kin_shots, binning_X and binning_Y during 
    # acquisition to prevent nasty threading / acquisition errors.
    _kin_visible = true
    
    # Linear size of the Pixis image data array.  Updated by in method 'acquire'.
    # Size in BYTES is computed and checked in 'lib_pixis.cpp'.
    img_array_size = CInt
   
    # Other camera and aquisition settings
    CCD_temp = Int

    view = View(Group(
            Group(  Item('exposure', width=-100), 
                    Item('gain', width=-100),
                    Item('selected_ADC_speed',
                        editor=EnumEditor(name='_ADC_list'), width=-100), 
                    Item('row_shift_time', width=-100,)
                 ), 
            spring,
            Group(  Item('num_kin_shots', width=-100),
                    Group( 
                        'binning_X',
                        spring,
                        'binning_Y',
                        orientation='horizontal',
                        ),
                    visible_when='_kin_visible',
                    show_border=True, label="Kinetics"),
            orientation='vertical'),
            )
    
    def open(self):
        """ Call to initialiase the camera.
        """
        pixis_interface.open()
        #return self.get_buffer_size()

    def close(self):
        """ Call to close the camera.
        """
        self.close_shutter() # redundant!
        return pixis_interface.close()

    def open_shutter(self):
        """ Opens the shutter right away.
        """
        return pixis_interface.open_shutter()

    def close_shutter(self):
        """ Close the shutter right away.
        """
        return pixis_interface.close_shutter()

    def acquire(self):
        """ Captures the sequence. This is a blocking call, it does not return
            until the images are captured and transfered. This uses the 
            attributs of the object to derive the parameters of the acquisition.
        """
        # Compute image array size, i.e. number of PIXELS, based on number
        # of CCD pixels and binning settings.
        
        # You better type this as integer, or when you pass it to C, 
        # you're in for surprises !
        self.img_array_size = (self._CCD_size/self.binning_X)*\
                            (self._CCD_size/self.binning_Y)
        #print "pixis_camera.acquire() -> image size = %d pixels" % self.img_array_size  
       
        image_data = zeros(self.img_array_size,dtype='i')
        # Calculate Y_size (the number of ROWS) before the call to acquire, 
        # so that if someone changes self.binning_X through the GUI during 
        # the acquisition will stil have the value with wich the acquisition 
        # was started.
        Y_size = float(self._CCD_size/self.binning_Y)
        not_error = pixis_interface.acquire(self.exposure, self.number_images, 
                            self.binning_X, self.binning_Y, self.gain, 
                            self._kin_win_size,\
                            self.ADC_speed[self.selected_ADC_speed],\
                            # convert float [microsec] into int [nanosec.]
                            (1000*self.row_shift_time).__int__(), image_data)
        # Now convert this to floats, or you're in for other surprises !
        return image_data.reshape((Y_size,-1)).astype('f')
        # 'image_data' now has 'Y_size' rows and 'buffer_size/Y_size' columns.

#    def get_buffer_size(self):
#        """ Sets the buffer size by querying lib_pixis. It is necessary to do
#            This before an acquisition, as python needs to provide the right
#            amount of memory to the C routines.
#        """
#        # We need an array with an int in it to be able to pass an int by
#        # reference to C.
#        buffer_size = array([0])
#        not_error = pixis_interface.get_buffer_size(self.number_images, \
#                        self.binning_X, self.binning_Y, buffer_size)
#        self.buffer_size = buffer_size[0]

    def triggered_open_shutter(self):
        """ Opens the shutter at the next trigger.
        """
        return pixis_interface.triggered_open_shutter()

    def triggered_close_shutter(self):
        """ Close the shutter at the next trigger.
        """
        return pixis_interface.triggered_close_shutter()


if __name__ == "__main__":
    # This does not get run when the file is used as a module
    from numpy import *
    from pylab import imshow, draw, show
    import pdb

    camera = PixisCamera()
    #pdb.set_trace()
    camera.open()
    data = camera.acquire()
    imshow(data)
    draw()
    #show()
    camera.close_shutter()
    camera.close()

# :vim:nocindent:
