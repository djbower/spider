#!/usr/bin/env bash

# sledgehammer
#spider -options_file reactions.opt  -atmosic_snes_verbose_monitor -atmosic_snes_view -atmosic_ksp_rtol 1e-8 -atmosic_ksp_atol 1e-8  -atmosic_snes_converged_reason -atmosic_snes_linesearch_damping 0.01 -atmosic_snes_max_it 10000 -snes_linesearch_max_it 1000

# Dan's playing around with the below
spider -options_file reactions.opt  -atmosic_snes_verbose_monitor -atmosic_snes_view  -atmosic_snes_converged_reason -atmosic_snes_atol 1.0E-6 -atmosic_snes_atol 1.0E-6 -atmosic_ksp_atol 1.0E-6 -atmosic_ksp_atol 1.0E-6 -atmosic_snes_max_it 10000 -snes_linesearch_max_it 1000
