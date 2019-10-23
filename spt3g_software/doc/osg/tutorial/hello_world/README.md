# Hello World Tutorial

This tutorial is meant to be a simple first step for submitting jobs. It does not depend on anything that is not provided by any OSG site, i.e. a shell. To use the tutorial just `condor_submit hello_world.submit` and take a look at the output files that will appear in the tutorial directory. 

As the next step, we want to output the log files to `/scratch/<username>` rather than the tutorial directory. This will require changing the lines starting with `Error`, `Output`, and `Log` accordingly. 