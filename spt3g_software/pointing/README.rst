--------
Pointing
--------

The pointing directory contains files related to calculating detector pointing (boresight, etc.).

Mapmaking procedures in standard processing (``std_processing``) used some of these functions, which were moved to ``python/pointing_for_mapmaking.py`` in May 2019.

Many files (in ``python``, two in ``scripts``) use sptpol software (and some hardcoded directories) and should be reviewed by pointing people before they are either: restructured and merged, archived, or deleted altogether. **In their current state they are not usable by 3G.**

In ``mapmaker/mapmakersutils.py`` there is a ``CalculatePointing`` module that was *not* moved to this pointing directory. This was because it involves some C++ code and would be tricky to move.

*The authors of this restructuring are Wei Q. and Melanie A. -- if you have questions or suggestions about the structuring here, or cannot find something, please contact us. May 2019*
