#/user/bin/env bash

# Template file for running SPIDER within the framework of a coupled interior/atmosphere/escape
# scheme

# ======================
# Description of options
# ======================
# -options_file bu_input.opts:     reads the options file to provide the parameters for SPIDER
#                                  subsequent options on the command line over-ride these defaults
# -IC_INTERIOR 2:                  read in interior ic from file (prescribed by -ic_interior_filename, see below)
# -ic_interior_filename 50008.json:JSON file to read in interior initial condition
# -IC_ATMOSPHERE 2:                read in atmosphere ic from file (prescribed by -ic_atmosphere_filename, see below)
# -ic_atmosphere_filename 50008.json:JSON file to read in atmosphere initial condition
# -SURFACE_BC 4:                   constant heat flux boundary condition
# -surface_bc_value 1.0E4:         prescribed surface heat flux (e.g., 10^4 W/m^2)
#                                  this value will be determined by SOCRATES and passed in
# -activate_rollback:              revert to previous timestep when an event has been triggered
#                                  (see H2O_poststep_change and CO2_poststep_change below)
# -activate_poststep:              run a poststep function to determine if an event has occurred
# -H2O_poststep_change 0.05:       fractional change in H2O partial pressure that triggers an event
#                                  (default here is 5% change, where the reference value is read in from the
#                                  ic file)
# -CO2_poststep_change 0.05:       as above for CO2, also a 5% change to trigger an event
# now, the options below effectively enable you to accommodate a maximum timestep before the code terminates, since
# you may want to update atmospheric escape (or another quantity) within a given time frame even if an event
# has not been triggered due to the change in volatile concentration.  Hence we are using a single time step
# to provide another 'event' that halts the code.  In essence, dtmacro prescribes the maximum time (in years)
# to advance before the code again halts even if an event has not occurred.
# -nstepsmacro 1:                  make a single timestep
# -dtmacro 30000:                  delta time to advance by, here is 30,000 years
#
# -tsurf_poststep_change:          maximum absolute temperature change in Kelvin

# ==================
# Overview of result
# ==================
# TODO: Dan has not updated the following in light of the tsurf_poststep_change criteria!
# data is read from 50008.json, and the time-stepper runs until the H2O partial pressure tration in the mantle exceed a
# 5% change.  Two files are output in the output/ directory.  The first, output/50008.json is the first output
# that reflects the new state of the system computed from the restart data. The second file output is 56710.json.
# This output is triggered by the fact that H2O partial pressure now exceeds a 5% change, relative to the value
# that was read in from 50008.json (you can easily see this by manually checking the values in the JSONs).

# Now, the coupler could either search for the last time output, and use this to restart, or we could generate
# a generically named restart file which is a clone of the last output.  Basically, this would allow us to have
# one file named 'restart.json' that contains the data of the last time step before an event occurred.  My
# general preference is to let the coupler / external python scripts duplicate and copy files rather than asking
# SPIDER to do it.

# generate the restart file from the initial condition
# this will output 50008.json in output/
#spider -options_file bu_input.opts

# main restart example is here:
spider -options_file bu_input.opts -IC_INTERIOR 2 -ic_interior_filename 50008.json -SURFACE_BC 4 -surface_bc_value 1.0E4 -IC_ATMOSPHERE 2 -ic_atmosphere_filename 50008.json -activate_rollback -activate_poststep -H2O_poststep_change 0.05 -CO2_poststep_change 0.05 -nstepsmacro 1 -dtmacro 30000
