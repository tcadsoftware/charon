# Charon | tcad-charon

[Charon](https://charon.sandia.gov) is an open-source semiconductor device
modeling code, widely referred to as a TCAD (technology computer-aided design)
code, developed at [Sandia National Laboratories](https://www.sandia.gov).  It
is written in C++ and relies on another Sandia open-source project,
[Trilinos](https://github.com/trilinos/Trilinos), for supporting code, such as
nonlinear and linear solvers, finite-element and finite-volume libraries, I/O,
etc.  In addition to running on most Linux-based computers, Charon also
supports simulation of extremely large problems on massively parallel computing
systems that support the MPI standard.

## Contents

1. [Getting Up and Running with Charon](#getting-up-and-running-with-charon)
   1. [Git Instruction](#git-instruction)
   1. [Git LFS](#git-lfs)
   1. [Set Up your Compiler Toolchain](#set-up-your-compiler-toolchain)
   1. [Clone the Repositories](#clone-the-repositories)
   1. [Configure, Build, and Test](#configure-build-and-test)
   1. [Updating Charon](#updating-charon)
1. [Where to Ask Questions](#where-to-ask-questions)
1. [Contributing to Charon](#contributing-to-charon)
1. [License](#license)



## Getting Up and Running with Charon

### Git Instruction

Charon uses the git version control system.  If you are unfamiliar with how to
use git and GitLab, the following resources will likely be beneficial:
*  [Our Git Cheat Sheet](https://cee-gitlab.sandia.gov/Charon/tcad-charon/wikis/Git-Cheat-Sheet)
*  [Introduction to Version Control and Collaboration with Git](https://sems-atlassian-son.sandia.gov/confluence/display/GIT1)
*  [Intermediate Git](https://sems-atlassian-son.sandia.gov/confluence/display/GIT2)

[↑ Contents](#contents)

### Git LFS

Charon also uses git LFS in many of its repositories.

> *You must have git LFS installed and configured on your machine before
> you'll be able to work with Charon.*

Visit [this wiki page](https://cee-gitlab.sandia.gov/Charon/tcad-charon/wikis/git-lfs)
to read up on what it is, how to install it, and how to use it.

[↑ Contents](#contents)

### Set Up your Compiler Toolchain

Center 1400's Software Engineering Maintenance and Support group
([SEMS](https://sems.sandia.gov)) provides almost all the compilers and
third-party libraries (TPLs) that Charon requires.  These are now widely
available on Sandia machines, including the SON, SRN, CEE LAN and HPC.  This is
our preferred standardized source for TPLs and they are consistent with our
testing environment.  See [this
page](https://sems-atlassian-son.sandia.gov/confluence/display/SEMSKB/SEMS+NFS+TPL+System)
for instructions on how to gain access to the SEMS modules, and if you need
help, [contact the SEMS team for
assistance](https://sems.sandia.gov/content/submit-service-request).

Once you have access to the modules, you'll need to load the following:
*  a compiler (gcc, intel, clang)
*  OpenMPI
*  CMake
*  Python
*  Boost
*  HDF5
*  NetCDF

[↑ Contents](#contents)

### Clone the Repositories

Getting all the Charon repositories set up in the right places is mostly
automated.  First, make sure you have SSH keys setup on cee-gitlab by following [these instructions](https://cee-gitlab.sandia.gov/help/ssh/README.md#adding-an-ssh-key-to-your-gitlab-account)

Next clone this repository with
```bash
cd <someBaseDirectory>
git clone git@cee-gitlab.sandia.gov:Charon/tcad-charon
```

and then run the `setup` script with
```bash
cd tcad-charon
./scripts/setup
```

Follow the prompts as the script sets everything up for you.  If you accept the
default build directory location, you should have the following directory
structure:
```
<someBaseDirectory>/ # a.k.a. ${WORKSPACE}
  build/
    configure
    tcad-charon/
  tcad-charon/
    ...
    docs/
    ...
    src/
    test/
      heavyTests/
      heavyTestsOUO/
      nightlyTests/
      nightlyTestsOUO/
      unitTests/
    ...
    Trilinos/
```

[↑ Contents](#contents)

### Configure, Build, and Test

To configure Charon `cd` to a `build` directory, either the one that
was created above or one of your own. It is generally a good idea to
create another subdirectory utilizing a name representative of the
build you will be performing. For example, if you're using the GNU
compilers and building with debug support
```
cd ../build
mkdir gnu.dbg
cd gnu.dbg
```
You can then configure with
```
../../tcad-charon/scripts/build/all/build_charon.py
```
If the configure is successful, you can build with
```
make -j <numProcs>
```
Assuming the build is successful, you can run the test suite with
```
ctest -j <numProcs> -L nightly
```
additionally, if you build a debug version you'll want to run ctest like
```
ctest -j <numProcs> -L nightly -LE debugexclude
```

For futher details on the charon build script, `build_charon.py`, see
its associated [wiki page](https://cee-gitlab.sandia.gov/Charon/tcad-charon/wikis/Using-the-python-build-script).


[↑ Contents](#contents)

### Updating Charon

Any time you need to grab the latest from all the Charon repositories:
```bash
cd path/to/tcad-charon
gitdist checkout develop
gitdist pull
```

The [`gitdist`](https://tribits.org/doc/TribitsDevelopersGuide.html#gitdist-dist-help-all)
command will run those git commands across all of Charon's repositories.  For
more details on what it is and how to use it, see [this wiki
page](https://cee-gitlab.sandia.gov/Charon/tcad-charon/wikis/Using-gitdist).

[↑ Contents](#contents)



## Where to Ask Questions

If you need help with Charon, feel free to ask questions by [creating an
issue](https://cee-gitlab.sandia.gov/Charon/tcad-charon/issues/new).  Select
the ~Question issue template, and that will pre-populate the *Description*
field, giving you instructions on submitting the issue.

[↑ Contents](#contents)



## Contributing to Charon

If you're interested in contributing to Charon, we welcome your collaboration.
Please read [our contributing
guidelines](https://cee-gitlab.sandia.gov/Charon/tcad-charon/blob/develop/CONTRIBUTING.md)
for details on our workflow, submitting merge requests, etc.

[↑ Contents](#contents)



## License

See [Copyright.txt](https://cee-gitlab.sandia.gov/Charon/tcad-charon/blob/develop/Copyright.txt).

[↑ Contents](#contents)
