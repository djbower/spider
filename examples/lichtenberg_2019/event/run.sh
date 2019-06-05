#/user/bin/env bash

spider -options_file bu_input.opts -initial_condition 2 -ic_filename 0.json -SURFACE_BC 4 -surface_bc_value 1.0E4 -SOLVE_FOR_VOLATILES 1 -activate_rollback -activate_poststep -H2O_poststep_change 0.05 -CO2_poststep_change 0.05
