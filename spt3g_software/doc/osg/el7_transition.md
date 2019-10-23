# Transition Guide to Scientific Linux 7 (SL7)

## Why transition to SL7?

The transition to SL7 is necessary for three reasons: Push towards RHEL7-only 
infrastructure at MWT2, Ceph optimizations, and having same OS at Pole as at
UChicago. MWT2 is pushing towards deploying RHEL7 derivates to all 
infrastructure that support it. In this case this means SL7 for bare-metal 
machines and CentOS7 for VMs. 

The worker nodes will transition to SL7 later this year. We haven't 
determined whether it will be a rolling upgrade or there will be as mass 
migration. That being said, the transition will require that SPT-3G code 
is both SL6 and SL7 compatible.

The last major release of Ceph we deployed in Dec. 2017 added new Ceph 
kernel optimizations. These optimizations require a C++1x-compatible 
system compiler. RHEL6 and its derivatives do not provide a 
C++1x-compatible system compiler. RHEL7 or a derivative do provide 
said compiler. With the optimizations turned on we should see better 
and more stable performance on `/spt` on scott and amundsen.

Having the same OS on both pole and UChicago machines will facilitate
overall development of SPT 3G software. A uniform environment will allow
analyzers to test possible code in a similar envionment and reveal
possible issues. 

## EL6 vs. EL7

Like any OS upgrade, this transition will bring major changes with it.
Overall the changes can be summarized as:

* Replace init-system with systemd`
* Python version 2.7
* GCC 4.8.5
* Kernel 3.10 by default (we typically run newer kernels on scott and amundsen)
* Newer glibc, libstdc++, etc.

## How to compile code using Singularity

Across the OSG there will both RHEL/CentOS/SL6 and RHEL/CentOS/SL7 machines.
There is a gradual transition to RHEL/CentOS/SL7 occuring. We expect most of
the sites to SL7 by middle of 2019. To get the most out of OSG, it is good practice to
allow for both SL6 and SL7 machines. 

To allow you to compile your code for both RHEL/CentOS/SL6 and RHEL/CentOS/SL7, we 
recommend using [Singularity](https://singularity.lbl.gov/#). Singularity is a 
containerization solution focused on scientific computing, especially HPC. 
Containers change the system environment to mimic another Linux 
distribution, while isolating your workload from the system.

To invoke a simple Singularity environment:

```
singularity shell docker://ubuntu:16.04
```

You will be met with the following shell

```
[briedel@login03 ~]$ singularity shell docker://ubuntu:16.04
Docker image path: index.docker.io/library/ubuntu:16.04
Cache folder set to /home/briedel/.singularity/docker
[5/5] |===================================| 100.0%
Creating container runtime...
Singularity: Invoking an interactive shell within container...

Singularity ubuntu:16.04:~>
```

When we check the release file we see

```
Singularity ubuntu:16.04:~> cat /etc/os-release
NAME="Ubuntu"
VERSION="16.04.4 LTS (Xenial Xerus)"
ID=ubuntu
ID_LIKE=debian
PRETTY_NAME="Ubuntu 16.04.4 LTS"
VERSION_ID="16.04"
HOME_URL="http://www.ubuntu.com/"
SUPPORT_URL="http://help.ubuntu.com/"
BUG_REPORT_URL="http://bugs.launchpad.net/ubuntu/"
VERSION_CODENAME=xenial
UBUNTU_CODENAME=xenial
```

It _looks_ and _feels_ like Ubuntu, but it really isn't. Looking at the 
reported kernel version, we see 

```
Singularity ubuntu:16.04:~> uname -r
3.10.0-693.5.2.el7.x86_64
```

which is the kernel of the underlying system (SL7), rather than the Ubuntu 
(Kernel version 4.4).

To generate an environment that looks  familiar pieces we recommend using

```
singularity shell --bind /cvmfs --bind /spt /cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el7:latest/ --login
```

For a good overview of the Singularity CLI, please take a look at the 
[Singularity Quickstart Guide](https://singularity.lbl.gov/quickstart).

## Create your own Container

In case you need to install specialized software that is not provided 
in CVMFS (e.g. newest version of astropy) or require some special
environment (e.g. FORTRAN 2018 compiler), building your own container
is the best option. The easiest way to "build" your own container is
to use a Dockerfile. A Dockerfile is a description of what
commands should be run as a container is "build". A build container is 
a tarball that includes an OS as well as your software.

Sample Dockerfile is available at in the [OSG GitHub repo](https://github.com/opensciencegrid/osgvo-el6/blob/master/Dockerfile). 
In the sample Dockerfilei, the first line

```
FROM centos:6
```

is the equivalent of an `import` statement in Python or `include` 
statement in C/C++. Docker will take image referenced after the `FROM` 
statement and import them as a base image. In this case we are using 
the CentOS 6 base image. One can also use Ubuntu 18.04 LTS with

```
FROM ubuntu:18:04
```

To install your software you will use typically use the respective Linux 
version's  package manager (e.g. `yum` for RHEL-derivaties or `apt` for
Debian derivaties). The basic BASH commands, such as `yum -y install python`, 
are run using the `RUN` keyword, such that

```
RUN yum -y install python-numpy
```

will install the numpy package from the respective OS's repo. In case the
standard repos do not have the package, you can add repos as well. 
For full documentation of the Dockerfile commands see [Dockerfile docs](https://docs.docker.com/engine/reference/builder/).

Once the Dockerfile is complete, you can test it in multiple ways. One is to 
link a repo in GitHub to your DockerHub account. This connection will allow
DockerHub automatically build a new version of the image for every push
to the repo. To allow other SPT users to see your image you can also add it 
to the `southpoletelescope` group on DockerHub. Please let us know if you
wanted to be added to the group. Once image has been added to the group, we
will add it to the OSG Singularity CVMFS repository and it will be distributed 
across OSG.

There is also an option to install software inside a container and then commit it
to DockerHub, but this is somewhat complicated compared to using a Dockerfile.
In case this is necessary please contact us for support.

