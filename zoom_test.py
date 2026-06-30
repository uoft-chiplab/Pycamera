# test some zoom crap
from numpy import *
from scipy import *
from matplotlib.figure import Figure
# Generate some 2D data
g = rand(100,100)

#fig = Figure()
#ax = fig.add_axes()
#ax.set_xlim((20,40))
#ax.set_ylim((30,60))
fig = gcf()
fig.clf()
#plot(g[:,0])
imshow(g,interpolation='nearest')
ax = gca()
ax.set_xlim((0,40))
ax.set_ylim((0,45))
draw()
show()
