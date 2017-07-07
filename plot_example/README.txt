1. The C-code outputs files in a directory called 'output' with files named dSdr_b_aug.TIMESTEP.m  These are the files containing the augmented data (value of dS/dr at the basic nodes and S0 at the uppermost staggered node).

2. globals_def.h essentially contains all the parameters of the model run

3. at the moment, you must recompile the code everytime you change the parameters in global_defs.h, since this header file is compiled as part of the code.  We will eventually change this to be an input file that is read in, much like the python script currently does.

4. after you have run the c-code, you need to copy the entire 'output' directory, which contains the result files, to your 'model' directory.  In the model directory there will be a magma_ocean.cfg file (the python configuration file) with EXACTLY THE SAME PARAMETERS as you used in your c model run.  Note that the python script may still contain non-dimensional values, whereas the c-code is now all dimensional.  This is OK, but the dimensional and non-dimensional entries must be equivalent between the magma_ocean.cfg and global_defs.h

5. you can then run the plotting script as before.  Basically, the script checks to see if the directory called 'output' exists.  If it does, it processes the dSdr_b_aug.TIMESTEP.m files and outputs the data into an 'out' directory.  Previously, an 'out' directory was created by the output of the data from the python script itself.

6. Assuming that magma_ocean.cfg has the same parameters as what you used to run the model (global_defs.h), then the plotting should be correct.

7. I recommend that you run a couple of models that you have already run using the python code, to confirm that the results are the same -- i.e., there are no bugs in the workflow when running with the C-code and then plotting using the workflow described above.
