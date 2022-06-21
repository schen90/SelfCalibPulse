#!/bin/bash
#./run.sh 1

./macros/MakeData_AddPS rootfiles/noPS_back/G4SimData$(printf "%04d" $1).root rootfiles/PS_back/G4SimData$(printf "%04d" $1).root
