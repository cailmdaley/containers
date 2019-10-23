----------
Quick Look
----------

Often when developing or diagnosing problems, it is useful to be able to interact with the G3Frames as they are being processed. The SPT3G software provides a utility for exploring and manipulating frames which can be used standalone or injected directly into a pipeline.

.. contents:: Contents

G3Viewer
========

Overview
--------

The G3Viewer is a curses application which can be injected directly into a pipeline to aid in development and diagnostics. It can also be called from standalone scripts if necessary. The G3Viewer can load any type of object, but provides some additional features when viewing a frame queue.

Usage
-----

The G3Viewer is simple to invoke. It can either be called standalone with

.. code-block:: python

    G3Viewer(object_to_explore)

or within a pipeline:

.. code-block:: python

	pipe.Add(G3Viewer)

This will start a curses window for each frame that passes through the pipeline. Hitting "q" will exit the curses window and the pipeline will continue until the next frame. A help menu can be accessed at any time by hitting '?'.

When using G3Viewer as a pipeline module, you can specify a subset of frame types to look at with the `frame_types` argument.

Navigation
----------

When an object is opened in the G3Viewer, the contents of the object will be shown on the left-most panel. If a list of G3Frames is passed, you will see the list of frames in the bottom window, labeled by the frame's type code, and only the first frame will be opened in the main viewer. To navigate between frames, hit shift+navigation (e.g. shift+right or shift+l).

The current selected object will be highlighted. If a member of the object can be explored further, you will see an arrow next to the object's description. Hitting a right-navigation key (right/l/d) or enter will open that object's contents in the next window. The right-most window is reserver for displaying the type of each member of the active window.

Object manipulation
-------------------

Sometimes it is useful to be able to actively manipulate the object you are exploring. When you hit 'o', you will be dropped into an interactive iPython console. The highlighted object will be stored as the variable 'selected', and the selected frame (if applicable) will be stored as 'selected_frame'. Exiting out of the iPython console will bring you back to the G3Viewer window.

InjectDebug
===========

Overview
--------

InjectDebug is a simple pipeline module which will initiate a pdb instance at the insertion point in the pipeline. The frame object passed to InjectDebug will be stored as `frame`.
