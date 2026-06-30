""" 
ThreadRunner object class: traited wrapper to a simple job that allows to
execute it in a seperate thread. Provides a traitsUI view of the thread.

"""
# Author: Gael Varoquaux <gael_dot_varoquaux_at_normalesup_dot_org>
# Copyright (c) 2006, Gael Varoquaux
# License: BSD Style.

from threading import Thread
from enthought.traits import *
from enthought.traits.ui import View, Item, ButtonEditor, spring, Group
try:
    from enthought.pyface.api import GUI
except ImportError:
    from enthought.pyface import GUI

class ThreadRunner( HasTraits, ):
    """ ThreadRunner(job=None)
        
        Class that run a job in a separate thread and provides a graphical
        representation to start or stop it.
    """
    # boolean giving the whether the thread is running or not.
    # should not be set by a method external to the object.
    _running = false

    # event that stops the thread if running, and vice versa.
    toggle = Event(desc="start/abort the job")

    button_go_label = Str('Start')

    button_abort_label = Str('Abort')

    _button_label = Str('Start')

    thread = Instance(Thread)

    _status = Enum('idle','running','aborting',desc="status of the job")

    view = View(Group(
                Item( 'toggle' , show_label=False,
                        # Not yet in current "production", too bad !
                        #editor = ButtonEditor(label_value = '_button_label')),
                        editor = ButtonEditor(), width=-55),
                #spring,
                Item( '_status', show_label=False, style='readonly', 
                        width=-100),
                orientation = 'horizontal', #layout = 'flow',
                 ), )


    # job(thread): callable that is executed by the thread, to be set
    #                 at run time, or at initialisation. It should accept one
    #                 argument: thread object in which it will live. To be able
    #                 to abort in the middle of the job, it should return when
    #                 thread.wants_abort is true.
    job = Callable(allow_none=False)

    def start(self):
        """ Call to start the thread if not running."""
        self.thread = Thread()
        self.thread.wants_abort = False
        def my_job():
            self.job(self.thread)
            self._running = False
        self.thread.run = my_job
        self._running = True
        self.thread.start()

    def abort(self):
        """ Call to abort the thread if running. """
        if self._running:
            self._status = 'aborting'
            self.thread.wants_abort = True

    def _job_default(self):
        def do_nothing(*args, **kwargs):
            pass
        return do_nothing

    def __running_changed(self, new_running):
        if new_running:
            GUI.invoke_later(setattr, self, '_button_label', 
                                            self.button_abort_label)
            self._status = 'running'
        else:
            GUI.invoke_later(setattr, self, '_button_label', 
                                            self.button_go_label)
            self._status = 'idle'

    def _toggle_fired(self):
        if self._running:
            self.abort()
        else:
            self.start()

