#!/usr/bin/env python
"""
An example of how to interact with the plotting canvas by connecting
to move and click events
"""
import sys
from pylab import *
#import pdb

#t = arange(0.0, 1.0, 0.01)
#s = sin(2*pi*t)
[X,Y] = meshgrid(arange(-20,20),arange(-20,20))
s = 2*exp( -(X-2)**2 /(2*3**2) - (Y+2)**2 /(2*4**2) )

ax = subplot(111)
#ax.plot(t,s)
ax.imshow(s,interpolation='nearest')

def on_move(event):
    # get the x and y pixel coords
    x, y = event.x, event.y

    if event.inaxes:
        ax = event.inaxes  # the axes instance
    #    print 'data coords', event.xdata, event.ydata

def on_click(event):
    #pdb.set_trace()
    # get the x and y coords, flip y from top to bottom
    x, y = event.x, event.y
    if event.button==1:
        print 'button press'
        if event.inaxes is not None:
            print 'data coords', event.xdata, event.ydata
            print 'matrix value = ', s[event.xdata,event.ydata]

binding_id = connect('motion_notify_event', on_move)
connect('button_press_event', on_click)

if "test_disconnect" in sys.argv:
    print "disconnecting console coordinate printout..."
    disconnect(binding_id)

show()

