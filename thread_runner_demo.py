"""Demo file to test the ThreadRunner object.

"""
# Author: Gael Varoquaux <gael_dot_varoquaux_at_normalesup_dot_org>
# Copyright (c) 2006, Gael Varoquaux
# License: BSD Style.
from enthought.traits.api import push_exception_handler
from enthought.debug.fbi import fbi

push_exception_handler ( handler = lambda o,t,ov,nv: fbi(),
                        main = True,
                        locked = True )

from thread_runner import ThreadRunner
from enthought.traits import *
from enthought.pyface.api import SplitApplicationWindow, GUI
from enthought.traits.ui import CodeEditor, View, Item

##########################################################################
# A sample application.
##########################################################################

class TextDisplay( HasTraits ):
    string =  String()
    view= View( Item('string',show_label=False, springy=True, style='custom' ))
    #view= View( Item('string',show_label=False, springy=True, style='readonly' ))

display = TextDisplay()

def set_display(value):
    global display
    GUI.invoke_later(setattr, display, 'string', value)

import time
def job(thread):
    set_display("\nJob started\n" + display.string)
    for i in xrange(0,20):
        if thread.wants_abort:
            return
        time.sleep(0.1)
        set_display("." + display.string)
    set_display("Job ended\n" + display.string)

thread_runner = ThreadRunner(job = job)

class MainWindow(SplitApplicationWindow):

    traits_panel = Instance(ThreadRunner)
    display = Instance(TextDisplay)

    def _create_rhs(self, parent):
        self.traits_panel = thread_runner
        return self.traits_panel.edit_traits(parent = parent,
                            kind="subpanel").control

    def _create_lhs(self, parent):
        self.display = display
        return self.display.edit_traits(parent = parent,
                            kind="subpanel").control

if __name__ == '__main__':
    gui = GUI()
    window = MainWindow()
    window.open()
    gui.start_event_loop()

