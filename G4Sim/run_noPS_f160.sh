#!/bin/bash
#./run_noPS_f160.sh 1

#obj=/home/userfs/s/sc2950/server/agata/trunk/build-gcc9.3.0/Agata
obj=/home/userfs/s/sc2950/server/agata/trunk/build-gcc7.4.0/Agata
mac=G4mac/sim1gamma_f160.mac

#($obj -Path ./trunk/ -b $mac -run $1)
($obj -seed -Path ./trunk/ -b $mac -run $1)

./macros/MakeDataG4_noPS trunk/GammaEvents.$(printf "%04d" $1) rootfiles/noPS_f160/G4SimData$(printf "%04d" $1).root

rm -f trunk/GammaEvents.$(printf "%04d" $1)
