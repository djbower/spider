#/user/bin/env bash

# updated version of ../socrates example for restarting the atmosphere from
# known partial pressures of volatile species

# Template file for running SPIDER within the framework of a coupled interior/atmosphere/escape
# scheme

# ======================
# Description of options
# ======================
# -options_file bu_input.opts:     reads the options file to provide the parameters for SPIDER
#                                  subsequent options on the command line over-ride these defaults
# -IC_INTERIOR 2:                  read in interior ic from file (prescribed by -ic_interior_filename, see below)
# -ic_interior filename 50003.json:JSON file to read in interior initial condition
# -SURFACE_BC 4:                   constant heat flux boundary condition
# -surface_bc_value 1.0E4:         prescribed surface heat flux (e.g., 10^4 W/m^2)
#                                  this value will be determined by SOCRATES and passed in
# -activate_rollback:              revert to previous timestep when an event has been triggered
#                                  (see H2O_poststep_change and CO2_poststep_change below)
# -activate_poststep:              run a poststep function to determine if an event has occurred
# -H2O_poststep_change 0.05:       fractional change in H2O concentration in the melt phase that triggers an event
#                                  (default here is 5% change, where the reference value is read in from the
#                                  ic file)
# -CO2_poststep_change 0.05:       as above for CO2, also a 5% change to trigger an event
# -IC_ATMOSPHERE 3:                read in partial pressures of volatiles in Pa
# -CO2_initial_atmos_pressure      partial pressure of CO2 in Pa
# -H2O_initial_atmos_pressure      partial pressure of H2O in Pa
# now, the options below effectively enable you to accommodate a maximum timestep before the code terminates, since
# you may want to update atmospheric escape (or another quantity) within a given time frame even if an event
# has not been triggered due to the change in volatile concentration.  Hence we are using a single time step
# to provide another 'event' that halts the code.  In essence, dtmacro prescribes the maximum time (in years)
# to advance before the code again halts even if an event has not occurred.
# -nstepsmacro 1:                  make a single timestep
# -dtmacro 30000:                  delta time to advance by, here is 30,000 years

# main restart example is here:
spider -options_file bu_input.opts -IC_INTERIOR 2 -ic_interior_filename 50003.json -SURFACE_BC 4 -surface_bc_value 1.0E4 -activate_rollback -activate_poststep -H2O_poststep_change 0.05 -CO2_poststep_change 0.05 -IC_ATMOSPHERE 3 -H2O_initial_atmos_pressure 101325 -CO2_initial_atmos_pressure 101325 -nstepsmacro 1 -dtmacro 30000