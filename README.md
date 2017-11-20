# SPIDER
*Simulating Planetary Interior Dynamics with Extreme Rheology*

## Usage
This code is not (yet) open source; you have been granted access by Dan J. Bower.
Until released under an open source license, the original authors retain all copyright.

You are free to create your own private fork of the code.
Please follow the standard development practice of forking and pull requests.  For example:

https://confluence.atlassian.com/bitbucket/forking-a-repository-221449527.html

https://blog.scottlowe.org/2015/01/27/using-fork-branch-git-workflow/

https://help.github.com/articles/fork-a-repo/

Whilst you are encouraged to develop the code for your own purposes, we expect that useful new feature(s) 
will be merged into the mainline following publication of your work.  You can issue a pull request for this.  
In the meantime, try and follow the structure and format of the existing code to ensure changes can be relatively easily merged.

And of course, please issue pull requests for bugs immediately.

## Citation

### 1. SPIDER code
Bower, D.J., P. Sanan, and A.S. Wolf (2017), Numerical solution of a non-linear conservation law applicable to the interior dynamics of partially molten planets, Phys. Earth Planet. Inter., doi: 10.1016/j.pepi.2017.11.004, accepted.

### 2. Use of MgSiO3 melt data tables (RTpress) within SPIDER
Wolf, A.S. and D.J. Bower (2017), An equation of state for high pressure-temperature liquids (RTpress) with application to MgSiO3 melt, Phys. Earth Planet. Inter.

## Building

1. Provide a working C compiler
2. Build Dependencies, either in (2a) double or (2b) quadruple precision
3. Build C Code

### 1. C Compiler

First, ensure that you have a suitable, working C compiler.

#### Aside for Mac OS X]
A simple way to proceed is to install GCC compilers to the default location from:
    http://hpc.sourceforge.net
    
This site contains precompiled binaries for versions of OS X, and will install everything into:
    `/usr/local`
    
Note, however, that binaries for High Sierra are currently not available.

#### Otherwise

If you use MacPorts, Homebrew, or apt, these can also be used to easily install GCC.

You can test that you have a working `gcc` C compiler by running the following from the command line

        echo '#include<stdio.h>' > t.c && echo 'int main(){printf("It seems to work!\n");}' >> t.c && gcc t.c && ./a.out && rm -f t.c a.out

If you plan to use quadruple precision, you need to ensure that you have an actual version
of GCC installed. This is harder than it should be, because Apple irresponsibly installs
something called "`gcc`" which is actually a wrapper for their own compiler. You can check by running

        gcc --version

If you see a message about "Apple LLVM", then you have the wrong compiler! You can also test with this command:

        echo '#include<stdio.h>' > t.c && echo '#include<quadmath.h>' >> t.c && echo 'int main(){printf("It seems to work!\n");}' >> t.c && gcc t.c && ./a.out && rm -f t.c a.out

### 2. Build Dependencies
Before you begin, comment out any existing references to `PETSC_DIR` and `PETSC_ARCH` in your
`.profile`, `.bash_profile`,`.bashrc`, etc. and clear these variables.

### 2a. Build Dependencies: Double Precision

In this case you must use @psanan's hacked version of
PETSc because this version uses the dense direct solver that is required. However,
you can allow PETSc to automatically download and install SUNDIALS
(with double support) so there is no need to install SUNDIALS separately.

1. Get the hacked version of PETSc

        cd /somewhere/to/install
        git clone https://bitbucket.org/psanan/petsc --depth=1 -b psanan/ts-sundials-quad-hack petsc-double-direct

2. Change directories

        cd petsc-double-direct/

3. Configure PETSc using the following command:[Aside] For a debug build, amend command above to use --with-debugging=1

        ./configure --with-debugging=0 --with-fc=0 --with-cxx=0 --with-cc=gcc --download-sundials --download-mpich --download-mpich --COPTFLAGS="-g -O3" --CXXOPTFLAGS="-g -O3"

4. Make.  PETSc's configure process, if successful, will end by printing out a command which you can copy and paste, e.g.

        make PETSC_DIR=/Users/dan/Programs/petsc/petsc-double-direct PETSC_ARCH=arch-darwin-c-opt all

5. Run tests. PETSc's make process will print this out for you to copy and paste e.g.

        make PETSC_DIR=/Users/dan/Programs/petsc/petsc-double-direct PETSC_ARCH=arch-darwin-c-opt test

6. Note the value of PETSC_ARCH for later (will look something like arch-xxx-yyy, arch-darwin-c-opt in the example above)

### 2b. Build Dependencies: Quadruple Precision


#### 2b.i Build SUNDIALS with quadruple precision

1. Clone our hacked version of SUNDIALS from the git repository.

       cd /somewhere/to/install
       git clone https://bitbucket.org/psanan/sundials-quad

2. Make sure that you have CMake available by typing `cmake --version`.  If this fails, then install CMake from your package manager (homebrew, macports, apt,..) or by following the instructions at cmake.org/download.

       sudo port install cmake           # MacPorts
       brew install cmake                # Homebrew
       sudo apt-get install cmake        # apt

3. Configure, build, and install SUNDIALS.

        cd sundials-quad
        mkdir install build
        cd build
        cmake ..
        ccmake .            # with apt you may need to install this separately

4.  Use the ccmake interface to set values similar to those below.
Make sure you type "c" to configure once you have entered these values, then "g" to generate and exit.
Note: specify the same C compiler you used to install PETSc (probably "gcc")

        CMAKE_C_COMPILER: gcc
        CMAKE_INSTALL_PREFIX: ../install
        EXAMPLES_INSTALL_PATH: ../install/examples
        SUNDIALS_PRECISION: quadruple

5. Make and install

        make && make install

#### 2b.ii Build PETSc with quadruple precision


1. Get the hacked version of PETSc and change directories:

        cd /somewhere/to/install
        git clone https://bitbucket.org/psanan/petsc --depth=1 -b psanan/ts-sundials-quad-hack petsc-quad-direct
        cd petsc-quad-direct

2. Configure PETSc using the following command.  Crucially, in the next step we point PETSc to the quadruple precision installation of SUNDIALS that we just created (change /somewhere/to/install to the place that you installed SUNDIALS). [Aside] For a debug build, amend the command to use `--with-debugging=1`

        ./configure --with-debugging=0 --with-fc=0 --with-cxx=0 --with-cc=gcc --with-precision=__float128 --with-sundials=1 --with-sundials-dir=/somewhere/to/install/sundials-quad/install --download-mpich --download-f2cblaslapack --COPTFLAGS="-g -O3" --CXXOPTFLAGS="-g -O3"

3. Make by copying the line (the following all assume an optimised build, i.e, 3a above) e.g.

        make PETSC_DIR=/Users/dan/Programs/petsc/petsc-quad-direct PETSC_ARCH=arch-darwin-c-opt all

4. Run tests by copying the provided line, e.g.

        make PETSC_DIR=/Users/dan/Programs/petsc/petsc-quad-direct PETSC_ARCH=arch-darwin-c-opt test

5. Note the value of `PETSC_ARCH` (will look something like `arch-xxx-yyy`) for later

### 3. C Code

1. In your environment, set `PETSC_DIR` and `PETSC_ARCH` to the PETSc installation that you wish to use (either the quad or double precision, as above)

        export PETSC_DIR=/somewhere/to/install/petsc-double-direct # or /somewhere/to/install/petsc-quad-direct
        export PETSC_ARCH=arch-xxx-yyy

2. Make (from this directory)

        make clean
        make -j

3. Test (TODO WIP)

       make test

You should now be ready to use the code!

## Plotting and data processing

There is a basic python script 'plot\_spider.py' that plots model output, although some parameters remain hard-coded.  When you run a model, output data is stored in an output directory in petsc binary format.  You can see in the source code the order of the vector data.  This order is then mimicked in the python script to access the data.  The script itself needs more clean up, but it should be fairly straightforward to reverse engineer to understand how the data is accessed.
