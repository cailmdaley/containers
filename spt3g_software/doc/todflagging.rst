------------------
Timstream Flagging
------------------

Flags are stored in a G3MapVectorString that has the form: 

.. code-block:: python

   {bolometer_id: ['BadCalibrator', 'Glitchy', 'ActuallyAKitten']}

The bolometer id maps to a list of reasons why the bolometer is flagged.  

When developing code to flag timestreams we are trying to not force 
flagging decisions on people using intermediate data products.
By that I mean, if we need to compute something in order to flag the timestreams, 
we should compute the thing, store it in the frame, and then use a separate, 
later, Module to actually set the flags. 

As an example, with the glitch finder, we have a module that finds the number of 
glitches above a threshold, and a seperate module for flagging based off of those stored values.

Remember that data products in the frames are supposed to be immutable.  
Do *NOT* do something like 

.. code-block:: python

   frame['Flags'][bolo_id].append("NoisyTimestream").  

*BAD BAD BAD, NO. SHAME. BAD*

You will need to copy the flags, append the flag you want, 
delete the previous flags in the frame and  then add the new flags.  
Because having you do that every time would be insane please use the add_flag
convenience function.
