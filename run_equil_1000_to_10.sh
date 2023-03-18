#!/bin/bash

# This is a script to run equilibritation on an md system, slowly relating the position restriction over a number of runs.

gmx grompp -f equil_memb_1000.mdp -c 13_GlyT2_POPC_CLR_water_ions_EM2.gro -r 13_GlyT2_POPC_CLR_water_ions_EM2.gro -n GlyT2_POPC_CLR.ndx -p GlyT2_POPC_CLR.top -o 14_GlyT2_POPC_CLR_EQ_1000.tpr  -maxwarn 2  2>&1 | tee 14_grompp.txt
gmx mdrun -v -deffnm 14_GlyT2_POPC_CLR_EQ_1000

gmx grompp -f equil_memb_500.mdp -c 14_GlyT2_POPC_CLR_EQ_1000.gro -r 14_GlyT2_POPC_CLR_EQ_1000.gro -n GlyT2_POPC_CLR.ndx -p GlyT2_POPC_CLR.top -o 15_GlyT2_POPC_CLR_EQ_500.tpr  -maxwarn 2  2>&1 | tee 15_grompp.txt
gmx mdrun -v -deffnm 15_GlyT2_POPC_CLR_EQ_500

