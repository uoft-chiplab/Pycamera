# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

PyCamera is a GUI application for **absorption imaging of ultracold atoms** with a
Princeton Instruments **Pixis 1024BR** CCD camera. It runs the camera in *kinetics
mode* (a single CCD exposure is sliced into several stacked sub-images), builds an
optical-density (O.D.) image from signal/reference/dark frames, fits the atom cloud,
and computes physical quantities (atom number, temperature, phase-space density, RF
splitting fractions) live, shot by shot. Data is saved as MATLAB `.mat` files.

This is **Python 2** code (Python 2.5-era). It uses `print` statements, old-style
`except Exception, inst:` syntax, `enthought.traits`/`enthought.pyface` (legacy
Enthought Traits, *not* the modern `traits`/`traitsui` packages), and **wxPython**
as the GUI backend. It will not run under Python 3 or a modern scientific stack.

## Running

There is no build system, package manifest, or test runner. The app is launched
directly:

```
python pycamera.py                 # main single-species absorption-imaging app
python pycamera-dualSpecies.py     # variant for dual-species (Rb + K) imaging
```

The real camera is only reachable on the Windows lab machine via `lib_pixis.dll`
(loaded through ctypes). To work **without hardware**, swap the camera import in
`experiment2.py`:

```python
from pixis_camera import PixisCamera as Camera   # real hardware (default)
#from mock_camera import PixisCamera as Camera    # software emulation
```

`mock_camera.py` generates synthetic gaussian clouds so the full GUI/analysis
pipeline can run on a dev machine (the rest of the app is hardware-agnostic).

The native side lives in `C_interface/` (`lib_pixis.cpp`, `lib_pixis.sln`) and is
built with Visual Studio against the PVCAM SDK (`pvcam.h`, `master.h`, `Pvcam32.lib`)
to produce `lib_pixis.dll`. There is no automated test suite; `test_analysis.py`
imports a `datanal` module that does not exist in this tree and is stale.

## Architecture

The app is built on the **Traits MVC pattern**: domain objects subclass `HasTraits`
and declare a `view = View(...)`; TraitsUI auto-renders the GUI from those views, and
assigning to a trait automatically fires `_<name>_changed` / `_<name>_fired` handlers.
Understanding the data flow means tracing these trait-change handlers, not explicit
call sites.

### Singleton objects (created in `pycamera.py`)
- `camera` (`pixis_camera.PixisCamera`) — Pixis settings (exposure, gain, ADC speed,
  binning, `num_kin_shots`) and the blocking `acquire()` call.
- `experiment` (`experiment2.Experiment`) — physical parameters (detuning, TOF, trap
  frequencies, species). Pure parameter holder.
- `splitting` (`experiment2.Splitting`) — pixel-box boundaries for RF-splitting
  analysis of two cloud regions.
- `analysis` (`experiment2.Analysis`) — the processing engine; holds the rolling
  `stored_data` "history" buffer and pointers to all the above plus the plot widgets.

### Threading model (see `documentation.py` for the original design notes)
Three concurrent "timelines":
1. **GUI event loop** — wxPython, started by `gui.start_event_loop()`.
2. **Acquisition thread** — `acquisition_job()` in `pycamera.py`, wrapped by
   `ThreadRunner` and started/aborted by the GUI "Start/Abort" button. It loops:
   `camera.acquire()` → slice frame into sub-images by `num_kin_shots` → `make_OD()`
   → optionally `savemat()` → assign to `analysis.current_data`.
3. **Analysis thread** — spawned fresh each shot by `Analysis._current_data_changed`
   (the trait handler that fires when `current_data` is set). Runs `do_analysis()`.
   If analysis is still running when the next frame arrives, the new thread is
   skipped ("skipped a beat") to avoid overload.

Cross-thread GUI updates **must** go through `GUI.invoke_later(...)` (Traits/wx are
not thread-safe). This pattern is used throughout for setting displayed text and
redrawing canvases.

### Kinetics-mode frame slicing (critical, fragile)
The single CCD readout is reshaped and sliced into `num_kin_shots` sub-images inside
`acquisition_job`. Supported values are 2, 4, 5, 6, 8, 10 — each is a separate
hard-coded branch with its own O.D. construction and choice of which sub-images to
save. The 5- and 6-shot branches use **hard-coded row boundaries** (1024 doesn't
divide evenly), so edits there must keep the row ranges summing to 1024.

### Analysis pipeline (`data_analysis2.py`)
`do_analysis()` orchestrates: `fit_prep`/`fit_prep2` (crop region + build 1D
horizontal/vertical profiles) → `seed_gaussian1D` → `gfit1D` (scipy `leastsq`
gaussian fit) → `crunch_params`/`crunch_split_params` (convert fit widths/areas into
physical units using hard-coded constants: pixel size `2.72e-6` m, magnification,
species linewidth `gamma`, ħ, k_B). `make_OD()` builds the optical-density image from
signal and reference (and their backgrounds).

### View / display layer
- `mplwidget.py` (`MPLWidget`) — wraps a matplotlib `Figure` on a wxAgg canvas; the
  base plotting widget.
- `matrix_view.py` (`MatrixView`, `StaticCursor`) — upper-left false-color image
  display with draggable axes (the drag limits define the analysis crop region) and a
  crosshair cursor.
- `pycamera.py` assembles the window: `MainWindow` (vertical split) → `LeftPanel`
  (image view over analysis plots) + `MyPanel` (tabbed TraitsUI controls for
  acquisition / experiment / camera / splitting).

### Hardware interface layers (two parallel copies — beware)
`pixis_interface.py` (used by `pixis_camera.py`) and `lib_pixis.py` are both ctypes
wrappers around `lib_pixis.dll`, with **different** DLL search paths and inverted
error-check logic. `pixis_camera.py` imports `pixis_interface`. Confirm which wrapper
you're touching before editing.

## Conventions and gotchas

- `os.environ['NUMERIX']='numpy'` is set at the top of most modules so legacy Traits
  uses numpy arrays — keep it when adding new modules.
- Data is written to a Windows network share; the base path is hard-coded as
  `base_dir` in `experiment2.py`. Auto-path generation walks dated subdirectories and
  will break if an unexpected folder appears (noted in-code).
- The repo contains many stale artifacts that are **not** the source of truth: `.pyc`,
  `.py~`, `.bak` files, `New Folder (2)/`, Windows `.lnk` shortcuts, and `doc/`
  (generated HTML API docs). Edit the plain `.py` files only.
- `experiment2.py` is the live module (despite the `2`); plain `experiment.py` does
  not exist here. Likewise `data_analysis2.py` is current.
