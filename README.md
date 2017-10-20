# SPIDER

## Building

1. Provide a working C compiler
2. Build Dependencies, either in (2a) double or (2b) quadruple precision
3. Build C Code

### 1. C Compiler

First, ensure that you have a suitable, working C compiler.

A simple way to proceed is install GCC compilers to the default location from:
         http://hpc.sourceforge.net

This will install everything into `/usr/local`

If you use MacPorts, Homebrew, or apt, these can also be used to easily install GCC.

You can test that you have a working `gcc` C compiler by running the following from the command line

        echo '#include<stdio.h>' > t.c && echo 'int main(){printf("It seems to work!\n");}' >> t.c && gcc t.c && ./a.out && rm -f t.c a.out

If you plan to use quadruple precision, you need to ensure that you have an actual version
of GCC installed. This is harder than it should be, because Apple irresponsibly installs
something called "`gcc`" which actually a wrapper for their own compiler. You can check by running

        gcc --version

If you see a message about "Apple LLVM", then you have the wrong compiler!

### 2. Build Dependencies
Before you begin comment out any existing references to `PETSC_DIR` and `PETSC_ARCH` in your
`.profile`, `.bash_profile`,`.bashrc`,etc. and re-source (clear these variables)
profile or bash_profile and resource (clear these variables)

### 2a. Build Dependencies: Double Precision

In this case you must use Patrick's hacked version of
PETSc because this version uses the dense direct solver that is required. However,
you can allow PETSc to automatically download and install SUNDIALS
(with double support) so there is no need to install SUNDIALS separately.

1. get the hacked version of PETSc (Patrick must give you read access to his
bitbucket - email patrick.sanan@gmail.com with your bitbucket username):

        cd /somewhere/to/install
        git clone https://bitbucket.org/psanan/petscfork -b psanan/ts-sundials-quad-hack petsc-double-direct

2. change directories

        cd petsc-double-direct/

3. configure PETSc using the following command:[Aside] For a debug build, amend command above to use --with-debugging=1

        ./configure --with-debugging=0 --with-fc=0 --download-sundials --download-mpich --with-cc=gcc --with-cxx=g++ --download-mpich --COPTFLAGS="-g -O3" --CXXOPTFLAGS="-g -O3"

4. make (the following all assume an optimised build, i.e, 3a above).
PETSc's configure process, if successful, will end by printing out a command which you can copy and paste, e.g.

        make PETSC_DIR=/Users/dan/Programs/petsc/petsc-double-direct PETSC_ARCH=arch-darwin-c-opt all

5. run tests. PETSc's make process will print this out for you to copy and paste e.g.:

        make PETSC_DIR=/Users/dan/Programs/petsc/petsc-double-direct PETSC_ARCH=arch-darwin-c-opt test

6. Note the value of PETSC_ARCH (will look something like arch-xxx-yyy) for later

### 2b. Build Dependencies: Quadruple Precision


#### 2b.i Build SUNDIALS with quadruple precision

1. Patrick needs to give you read access to his bitbucket to obtain
   the quad-supported version of SUNDIALS

2. Clone from git repository:

       cd /somewhere/to/install
       git clone https://bitbucket.org/psanan/sundials-quad.git sundials-quad

3. Make sure that you have CMake available by typing `cmake --version`.  If this fails, then install CMake from your package manager (homebrew, macports, apt,..) or by following the instructions at cmake.org/download

       sudo port install cmake           # MacPorts
       brew install cmake                # Homebrew
       sudo apt-get install cmake        # apt

4. Configure, build, and install SUNDIALS

        cd sundials-quad
        mkdir install build
        cd build
        cmake ..
        ccmake .            # with apt you may need to install this separately

5.  Use the ccmake interface to set the values below.
Make sure you type "c" to configure once you have entered these values, then [q]uit.  
Note: specify the same C compiler you used to install PETSc (probably "gcc")
Note: you may also directly edit `CMakeCache.txt` to edit these values, followed by "cmake .." again

        CMAKE_C_COMPILER: /usr/local/bin/gcc
        CMAKE_C_FLAGS: -O3
        CMAKE_INSTALL_PREFIX: ../install
        EXAMPLES_INSTALL_PATH: ../install/examples
        SUNDIALS_PRECISION: quadruple

6. Make and install:

        make && make install

(Note: you will see many warnings. TODO fix this)

#### 2b.ii Build PETSc with quadruple precision


1. get the hacked version of PETSc (@psanan must give you read access to the
bitbucket repository) and change directories:

        cd /somewhere/to/install
        git clone https://bitbucket.org/psanan/petscfork -b psanan/ts-sundials-quad-hack petsc-quad-direct
        cd petsc-quad-direct

2. configure PETSc using the following command.  Crucially, in the next step we point PETSc to the quadruple precision installation of SUNDIALS that we just created (change /somewhere/to/install to the place that you installed SUNDIALS). [Aside] For a debug build, amend command above to use `--with-debugging=1`

        ./configure --with-debugging=0 --with-fc=0 --with-precision=__float128 --with-sundials=1 --with-sundials-dir=/somewhere/to/install/sundials-quad/install --download-mpich --download-f2cblaslapack --with-cc=gcc --with-cxx=g++ --COPTFLAGS="-g -O3" --CXXOPTFLAGS="-g -O3"

3. make by copying the line (the following all assume an optimised build, i.e, 3a above) e.g.

        make PETSC_DIR=/Users/dan/Programs/petsc/petsc-quad-direct PETSC_ARCH=arch-darwin-c-opt all

4. run tests by copying the provided line, e.g.

        make PETSC_DIR=/Users/dan/Programs/petsc/petsc-quad-direct PETSC_ARCH=arch-darwin-c-opt test

5. Note the value of `PETSC_ARCH` (will look something like `arch-xxx-yyy`) for later

### 3. C Code

1. In your environment, set `PETSC_DIR` and `PETSC_ARCH` to the PETSc installation that you wish to use (either the quad or double precision, as above)

        export PETSC_DIR=/somewhere/to/install/petsc-double-direct # or /somewhere/to/install/petsc-quad-direct
        export PETSC_ARCH=arch-xxx-yyy

2. make (from this directory)

        make clean
        make -j

3. test (TODO WIP)

       make test

You should now be ready to use the code!
