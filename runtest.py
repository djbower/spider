#!/usr/bin/env python

import os


cmd_l = ['make clean','make -j CFLAGS_EXTRA="-O0"','lldb -- spider -options_file tests/opts/blackbody_bc.opts']
#cmd_l = ['make clean','make -j CFLAGS_EXTRA="-O0"','lldb -- spider -options_file tests/opts/reaction.opts -nstepsmacro 0']

for cmd in cmd_l:
    print(cmd)
    os.system(cmd)

