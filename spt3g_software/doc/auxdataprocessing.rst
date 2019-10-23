Aux Data Processing
===================

For SPT3G, the aux data observations are separated into their own files.  For most of the aux data processing scripts there are no depedencies.  In that case the auxdata processing should just be a function that accepts a filename and some other generic arguments and saves the data in some file.  The actual format of the file it saves it in doesn't matter as long as it can transfered between machines.  Pickle, JSON, HDF5, text, and fits are all fine.  Binary files with data dumped to it is not fine.  If you store it in XML, I hate you.  You are encouraged to save it in JSON, txt, or pickle files.  The output file should have the relevant information from the processing.  Once that processing step has been written we just need to write a module that jams the information into the frame.  Nathan and Nick have a lot of opinions about how the information ends up in the frame so please talk to them about it.

The only aux data processing with dependencies is the source observation processing and Nathan is doing that so he is in charge.

Coding Style
------------

1) Do not hard code keys in the frame, they are likely to change.  Store them as keyword argument to the function.  ``frame["yourkey"]`` BAD.
#) Try to limit library dependencies.  We are very likely to run this code on clusters in the future and fewer dependencies are better.
#) Do not have any dependencies on sptpol_software.  If you need code from that library import it into spt3g.
#) If you have any questions about what to do, talk to Nathan or Nick.  We care a lot about making the code easy to maintain, work with and use.  If code gets committed that is going to cause problems in the future, we will ask you to change it.

