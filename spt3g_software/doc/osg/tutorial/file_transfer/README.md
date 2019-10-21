# File Transfer Tutorial

This tutorial is meant to be familiarize you with getting the SPT 3G software dependencies via CVMFS and transferring input/output files. Before submitting anything, you will need to have your grid certificate in place, i.e. a `usercert.pem` and `userkey.pem` file in `$HOME/.globus/`. With your grid certificate in place, please execute

`grid-proxy-init -valid 168:00` 

This will generate a X509 user proxy for you. It is a temoprary authentication file that will allow you top access `/spt/` from a remote site. 

Now that you have the proxy in place, we can to take a quick look at `file_transfer.sh`. There are two lines that are crucial:

```
# Getting the SPT software dependencies
eval `/cvmfs/spt.opensciencegrid.org/py2-v1/setup.sh`

# Copying the file to the local directory
globus-url-copy -vb $1 file://$PWD/
```

We can submit the job:

`condor_submit file_transfer.submit`

As the next step, we want to output the log files to `/scratch/<username>` rather than the tutorial directory. This will require changing the lines starting with `Error`, `Output`, and `Log` accordingly. 