###############################################################################
# lib_pixis.py                                  Gael Varoquaux, 10/06/06
# ############################################################################
# Python file creating the interface between lib_pixis.dll, coded in C, and 
# python. This uses ctypes with numpy. See http://scipy.org/Cookbook/Ctypes2
# for a minimal exemple
#
# The basic idea of the interfacing is that the input and output type of each
# C function are defined, then a C wrapper function is created.
#
# I wouldn't recommend modifying this if you don't know what you are doing.
#
###############################################################################

import numpy
import ctypes

# Load the C library
_lib_pixis = numpy.ctypeslib.load_library('lib_pixis', '.')
#_lib_pixis = numpy.ctypeslib.load_library('lib_pixis', '.\\C_interface\\Debug')

## Define the bar function
#_lib_pixis.bar.restype = ctypes.c_int
#_lib_pixis.bar.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int]
#def bar(x):
#    return _lib_pixis.bar(x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(x))
#x = numpy.random.randn(10)
#n = bar(x)

_lib_pixis.init.restype = ctypes.c_int
def init():
    """ Calls lib_pixis to initialise the camera.
    """
    error = _lib_pixis.init()
    if error:
        print "The camera initialisation returned an error, are you sure the camera is on ?"
        print "Error number was: %d" % error
    return error

_lib_pixis.close.restype = ctypes.c_int
def close():
    """ Calls lib_pixis to close the camera.
    """
    error = _lib_pixis.close()
    if error:
        print "The camera closing returned an error, are you sure the camera is initialized ?"
        print "Error number was: %d" % error
    return error

_lib_pixis.triggered_open_shutter.restype = ctypes.c_int
def triggered_open_shutter():
    """ Calls lib_pixis to open the shutter at next trigger.
    """
    error = _lib_pixis.triggered_open_shutter()
    if error:
        print "Couldn't open the shutter"
        print "Error number was: %d" % error
    return error

_lib_pixis.triggered_close_shutter.restype = ctypes.c_int
def triggered_close_shutter():
    """ Calls lib_pixis to close the shutter at next trigger.
    """
    error = _lib_pixis.triggered_close_shutter()
    if error:
        print "Couldn't close the shutter"
        print "Error number was: %d" % error
    return error

if __name__ == '__main__':
    print "A few tests to see if the camera is working\n"
    print "testing init()"
    init()
    print "testing triggered_open_shutter()"
    triggered_open_shutter()
    print "testing triggered_close_shutter()"
    triggered_close_shutter()
    print "testing close()"
    close()
    raw_input("Press a key to end")


# :vim:nocindent:
