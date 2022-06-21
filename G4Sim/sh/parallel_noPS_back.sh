#!/bin/bash
#./parallel_noPS_back.sh 0 10

Nproc=16  # <----------- max parallel procs

function PushQue {  # push PID into Que
    Que="$Que $1"
    Nrun=$(($Nrun+1))
}

function GenQue {  # update Que
    OldQue=$Que
    Que=""; Nrun=0
    for PID in $OldQue
    do
	if [[ -d /proc/$PID ]]; then
	    PushQue $PID
	fi
    done
}

function ChkQue {  # check Que
    OldQue=$Que
    for PID in $OldQue
    do
	if [[ ! -d /proc/$PID ]]; then
	    GenQue; break
	fi
    done
}

# loop all jobs
for i in `seq $1 $2`
do
    ./run_noPS_back.sh $i &  # <---------- CMD
    PID=$!
    PushQue $PID
    sleep 5
    while [[ $Nrun -ge $Nproc ]]
    do
	ChkQue
	sleep 1
    done
done
wait
