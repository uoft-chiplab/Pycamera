from __future__ import division
from mock_camera import PixisCamera as Camera
from numpy import *
from datanal import gfit1D, analyse

camera = Camera()

camera.get_buffer_size()
data = camera.acquire()

# As we are using Kinetics mode, the image actually is two images.
Xsize=int(data.shape[0]/2)
data = {'image1':data[:Xsize,:], 'image2':data[Xsize:,:]}
data['absorption'] = log(data['image1']/data['image2'])

params, h_fit, v_fit, h_xdata, h_gdata ,v_xdata, v_gdata = analyse(data)

