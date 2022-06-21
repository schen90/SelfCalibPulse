#!/bin/bash
#./run.sh 1

./macros/MakeData_AddPS rootfiles/noPS/G4SimData$(printf "%04d" $1).root rootfiles/PS/G4SimData$(printf "%04d" $1).root
