"""
This is a file used to create the documentation, I also use it to put
project wide documentation.

Flow control in pycamera
=========================

pycamera is a multithreaded program. It runs an `event-loop`_ that
triggers the execution of GUI driven events, as in most of the GUI
programs. Along side with the handling of the GUI we want pycamera to run
an acquisition loop: be listening to information coming from the camera in
order to collect the data and drive the camera. Finally we want pycamera
to process the data presented by the camera.

The program has three different "timelines" to keep: be responsive to the
GUI, be responsive to the camera, and complete its processing job. The
responsivness to the GUI event is taken care of by the event-loop. The two
other tasks are implemented in separate `threads`_. A thread is a "piece
of program" that lives in its own timeline. The acquisition loop is
implemented in a thread that is spawned by a button on the GUI. It loops
untils told to stop by the GUI. Each time it receives a set of images from
the camera it spawns a new thread that processes the data. If the
computation done by this analysis thread taken to long it will not be
finished when a new set of data arrives. The program then complains and
does not spawn a new thread, to avoid overloading the computer.

Talking to the camera: the C interface
======================================

The camera is accessed through PVCAM, the C library. A small wrapper has
been written in C. The goal of this wrapper is to expose simple operations
that can be perfomed on the camera to the python layer. It is interfaced
with python using the ctypes modules [#]_. A first file (
`pixis_interface.py`) is created that provides a mapping of the C
functions into python. This help the debugging of the C code. A second
file ( `pixis_camera.py`) uses the functions accessible throught this
interface to create a "pythonic" object that provides high-level
operations to the rest of the program, and a representation of itself for
the GUI.

Design philosophy
=================

Maintenance of a software system represents a significant portion of its
life-cycle. Developpers spend far more trying to understand and modify
programming artifacts than creating them in the first place 
(Goldberg 87 [#]_ ). Designing modular code helps localizing in the notions
necessary to understand the program: by looking at the code by itself the
developper can understand it, and does not have to encompass the whole
application.

Pycamera makes heavy use of `object-oriented programming`_ (OOP), and
tries to achieve as much modularity as possible. The use of OOP helps for
this (see this `reference`_), but in the design of the differentlayers
and objects it is important to separate the hadrware interface (camera
specific) from the GUI layer, the flow control layer, and most important
the physics. The use of an interface file to the low-level library for
the hardware enforce the separation of the camera related code. The heavy
use of automatically generated graphical interfaces created by
`TraitsUI`_ allows the programmer not to worry about the user interface,
and limits the number of "callbacks" to create manually.

**It is paramount that the program keeps this code separation**, as it
allows to modify the program for a new experiment or to add a new camera
without having to understand it completely.

Kown bugs
==========

* I have not been able to control the shutter of the camera. The camera
  is not shut down cleanly, the event-loop often closes the camera while
  the acquisition thread is still running and trying to talk to the
  camera. This causes un ugly exit (probably a crash on exit), but should
  not raise any problems as long as no other operations should be
  performed on exit.

* The time-series of the data acquired fill in a (purposly) limited
  buffer. When this buffer is full the GUI does not report an error but
  the graphs freeze. There is no obvious way of cleaning the buffer. An
  event (represented as a button) should be added to the acquisition
  object.

* Controling the gain does not work. A certain amount of work in C is
  required to make this work. See page 8 of the manual.

* I do not think I was able to reach very fast exposures. I think the
  current unit of exposure time is milliseconds, and it is an integer !
  There is probably a way to tell the camera to switch units. More work
  to be done in C !

.. [#] Good and simple examples of to use ctypes with numerical python
       can be found on the `scipy wiki`_

.. [#] A. Goldberg: Programmer as reader. *IEEE Software*, **4** (5): 62-70,
       September 1987

.. _`event-loop` : http://en.wikipedia.org/wiki/Event_loop

.. _`threads` : http://en.wikipedia.org/wiki/Thread_(computer_science)

.. _`object-oriented programming` : http://en.wikipedia.org/wiki/Object-oriented_programming

.. _`reference` : http://www.webopedia.com/TERM/O/object_oriented_programming_OOP.html

.. _`TraitsUI` : http://code.enthought.com/traits/

.. _`scipy wiki` : http://www.scipy.org

"""

from enthought import endo
from os import system, sep
from glob import glob

file_list = " ".join(glob("*.py"))
endo_path = endo.__path__[0] + sep + "scripts" + sep + "endo.py "
system(endo_path  + " --rst "+ file_list)

