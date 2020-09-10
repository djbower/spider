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

Whilst you are encouraged to develop the code for your own purposes, we expect that useful new feature(s) will be merged into the mainline prior to or following publication of your work.  You can issue a pull request for this.  In the meantime, try and follow the structure and format of the existing code to ensure changes can be relatively easily merged.  And of course, please issue pull requests for bugs immediately.

## Citation

### 1. SPIDER code
Bower, D.J., P. Sanan, and A.S. Wolf (2018), Numerical solution of a non-linear conservation law applicable to the interior dynamics of partially molten planets, Phys. Earth Planet. Inter., 274, 49-62, doi: 10.1016/j.pepi.2017.11.004, arXiv: https://arxiv.org/abs/1711.07303, EarthArXiv: https://eartharxiv.org/k6tgf

### 2. Use of MgSiO3 melt data tables (RTpress) within SPIDER
Wolf, A.S. and D.J. Bower (2018), An equation of state for high pressure-temperature liquids (RTpress) with application to MgSiO3 melt, Phys. Earth Planet. Inter., 278, 59-74, doi: 10.1016/j.pepi.2018.02.004, EarthArXiv: https://eartharxiv.org/4c2s5

### 3. Use of volatiles / atmosphere coupling
Bower, D.J., Kitzmann, D., Wolf, A.S., Sanan, P., Dorn, C., Oza, A.V. (2019), Linking the evolution of terrestrial interiors and an early outgassed atmosphere to astrophysical observations, Astron. Astrophys., 631, A103, doi: 10.1051/0004-6361/201935710, arXiv: https://arxiv.org/abs/1904.08300

## Building

1. Provide a working C compiler
2. Build Dependencies, either in (2a) double or (2b) quadruple precision
3. Build C Code

### 1. C Compiler

You need a "pure" GCC compiler for quadruple precision calculations.
Unfortunately on a Mac, "`gcc`" is a wrapper for the default Apple compiler ("`clang`") which is not a pure GCC compiler package (and hence does not support quadruple precision math).
However, "`gcc`" ("`clang`") should work fine for double precision calculations.  Nevertheless, particularly for Mac OSX you may choose to install GCC to support both double and quadruple precision calculations using a single compiler.

A basic test to ensure you have a working compiler is (note you can swap out "`gcc`" for another compiler binary name):

        echo '#include<stdio.h>' > t.c && echo 'int main(){printf("It seems to work!\n");}' >> t.c && gcc t.c && ./a.out && rm -f t.c a.out

#### 1a. GCC compiler for both double and quadruple precision 

If you use MacPorts, Homebrew, or apt, these can be used to easily install GCC.  For example, using MacPorts:

        sudo port install gcc8 

[DB: add instructions for Homebrew, but since I don't use Homebrew I can't test the commands].  This will install a set of GCC binaries, typically in "`/opt/local/bin`".  The C compiler (for GCC8) will be called "`gcc-mp-8"` where the mp is obviously clarifying that it was installed by MacPorts.
Now on a Mac, you will usually access the default Apple compiler ("`clang`") using "`gcc`", but you can easily access the MacPorts GCC you just installed by using "`gcc-mp-8`" instead.
So when you are installing the software in the next sections, just use "`gcc-mp-8`" instead of "`gcc`" when you are asked to specify the C compiler.

If you plan to use quadruple precision, you must have an actual version of gcc installed (e.g., using MacPorts as above; the Apple compiler will not work). You can check by running (for example):

        gcc --version
        gcc-mp-8 --version

If you see a message about "Apple LLVM", then you have the wrong compiler (i.e., quadruple precision is not supported.)

You can also test with this command (replace `gcc` by `gcc-mp-8` to test the MacPorts compiler):

        echo '#include<stdio.h>' > t.c && echo '#include<quadmath.h>' >> t.c && echo 'int main(){printf("It seems to work!\n");}' >> t.c && gcc t.c && ./a.out && rm -f t.c a.out

In the following sections, it is assumed that you have installed "`gcc-mp-8`", but you can swap this out for your preferred compiler. 

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

3. Configure PETSc using the following command. For a debug build, use --with-debugging=1 instead.

        ./configure --with-debugging=0 --with-fc=0 --with-cxx=0 --with-cc=gcc-mp-8 --download-sundials --download-mpich --COPTFLAGS="-g -O3" --CXXOPTFLAGS="-g -O3"

4. Make.  PETSc's configure process, if successful, will end by printing out a command which you can copy and paste, e.g.

        make PETSC_DIR=/Users/dan/Programs/petsc/petsc-double-direct PETSC_ARCH=arch-darwin-c-opt all

5. Run tests. PETSc's make process will print this out for you to copy and paste e.g.

        make PETSC_DIR=/Users/dan/Programs/petsc/petsc-double-direct PETSC_ARCH=arch-darwin-c-opt test

6. Note the value of PETSC_ARCH for later (will look something like arch-xxx-yyy, arch-darwin-c-opt in the example above)

### 2b. Build Dependencies: Quadruple Precision


#### 2b.i Build SUNDIALS with quadruple precision

1. Clone our hacked version of SUNDIALS from the git repository.

        cd /somewhere/to/install
        mkdir src
        git clone https://bitbucket.org/psanan/sundials-quad src

2. Make sure that you have CMake available by typing `cmake --version`.  If this fails, then install CMake from your package manager (homebrew, macports, apt,..) or by following the instructions at cmake.org/download.  Note: it is necessary to comment out the cmake_policy(SET CMP0042 NEW) in src/CMakeLists.txt if you are using an old version of cmake.

        sudo port install cmake           # MacPorts
        brew install cmake                # Homebrew
        sudo apt-get install cmake        # apt

3. Configure, build, and install SUNDIALS.

        mkdir install build
        cd build
        cmake ../src
        ccmake .            # with apt you may need to install this separately

4.  Use the ccmake interface to set values similar to those below.
Make sure you type "c" to configure once you have entered these values, then "g" to generate and exit.
Note: specify the same C compiler you used to install PETSc (probably "gcc")

        CMAKE_C_COMPILER: gcc-mp-8
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

2. Configure PETSc using the following command.  Crucially, in the next step we point PETSc to the quadruple precision installation of SUNDIALS that we just created (change /somewhere/to/install to the place that you installed SUNDIALS). For a debug build, amend the command to use `--with-debugging=1`

        ./configure --with-debugging=0 --with-fc=0 --with-cxx=0 --with-cc=gcc-mp-8 --with-precision=__float128 --with-sundials=1 --with-sundials-dir=/somewhere/to/install/sundials-quad/install --download-mpich --download-f2cblaslapack --COPTFLAGS="-g -O3" --CXXOPTFLAGS="-g -O3"

3. Make by copying the line (the following all assume an optimised build, i.e, 3a above) e.g.

        make PETSC_DIR=/Users/dan/Programs/petsc/petsc-quad-direct PETSC_ARCH=arch-darwin-c-opt all

4. Run tests by copying the provided line, e.g.

        make PETSC_DIR=/Users/dan/Programs/petsc/petsc-quad-direct PETSC_ARCH=arch-darwin-c-opt test

5. Note the value of `PETSC_ARCH` (will look something like `arch-xxx-yyy`) for later

### 3. SPIDER C Code

1. In your environment, set `PETSC_DIR` and `PETSC_ARCH` to the PETSc installation that you wish to use (either the quad or double precision, as above)

        export PETSC_DIR=/somewhere/to/install/petsc-double-direct # or /somewhere/to/install/petsc-quad-direct
        export PETSC_ARCH=arch-xxx-yyy

2. Get the SPIDER C code:

        cd /somewhere/to/install
        git clone https://djbower@bitbucket.org/djbower/spider-dev.git spider-dev
        cd spider-dev

2. Make (from this directory)

        make clean
        make -j

3. Test basic functionality:

        make test

4. Test more functionality:

        make test_all

You should now be ready to use the code!

### 4. Make your own working branch

You should prototype/develop your models using your own branch (not master).  So once you have verified that the code passes the inbuilt tests (above), you should immediately create and switch to your own working branch.  You can find information about how to create your own branch here:

https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging

The best naming strategy is to set up a branch called [username]/[feature] where username is your bitbucket account username and feature is the feature you are working on.  So for example, I could use:

        git checkout -b djbower/alieninvasion

Note that you create your branch locally, but then need to push the branch to the remote repository.  The guide above, and helpful onscreen messages, will help you with this process (if you are unsure, ask Dan).  In essence, you should only commit changes to your own branch.  Once you are happy that your branch is working and wish to contribute your feature to the mainline of the code, please issue a pull request (or ask Dan for advice).

### 5. SPIDER manual

There is a manual located in ```manual/```, and you are strongly encouraged to update this manual with new features and notes as you use and/or develop the code.  The document is very much a work in progress, so if needs be just create a new section or subsection for your work and consolidate your writing there.  Over time, I will continue to merge the material into consistent chapters.

### 6. Code development tips

1. Ensure all previous tests pass (i.e., you haven't broken a part of the code that was previously working)
1. Use valgrind to confirm that your code is free of memory leaks (or alternatively, we will do this before merging your branch with master)
1. Most importantly, comment comment comment!  And try to predict potential bugs as you code and add a commment about those as well!

   - Comments are no substitute for well-structured and clear code (code can't lie).
   - Do not include any (new) "TODO", "FIXME", or other personal comments in master. Fix problems or leave a clear comment or error message.

## Running the code

Add the installation directory to your `$PATH` so you can call the `spider` binary without requiring the absolute path.  Typically, you will set up a project folder elsewhere to store your models, i.e., outside the installation directory where the SPIDER source code resides.  You can run `spider` without an argument and standard parameters will be used for a simple magma ocean cooling model.  But in general, you will add an argument, i.e. `spider -options_file input.opts` to use the parameters specified in `input.opts` (this input filename is an example).  There are example input files in the `examples/` directory.

Some of the following tutorial examples have been moved.  Instead, look in the tests/opts directory to find example model input (opts) files.

## Tutorial example 1: Simple magma ocean cooling

1. Add the installation directory to your `$PATH`
1. Add the `py3` directory to your `$PATH` and `$PYTHONPATH`.   Note that the plotting scripts require a Python 3 installation.  I recommend using anaconda to install a python 3 environment that is completely separate from any python versions that may have been installed by your operating system.
1. `cd examples/bower_2019/blackbody`
1. `mkdir output`
1. `spider -options_file bu_input.opts`
1. `plot_spider.py -t 0,100,200,400,804,1200,1500`

where `-t` is the argument to specify output times (in years) that exist in the output directory (`output/`) of the model in JSON format.  Note that the output JSON format can be read by a simple text (ascii) reader.

## Tutorial example 2: Coupled interior-atmosphere evolution of a magma ocean

This example is a prototype, because the transition from melt to solid at the surface needs careful treatment (work in progress).

1. `cd examples/bower_2019/atmosphere`
1. `mkdir output`
1. `spider -options_file bu_input.opts`
1. `plot_spider.py -t 0,201313,500030,1000122`
1. `plot_atmosphere.py`

## Note on Python scripts

Python scripts are in `py3/` and can be used as templates for building your own analysis and plotting scripts.  Do not modify the scripts that are already in `py3/` unless you are fixing bugs or retaining backwards compatibility.  Instead, copy one of the existing scripts, rename, and then modify as you require.  You can of course check your plotting scripts into the repository, if you wish.
