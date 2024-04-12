#!/bin/bash -x

narr=512
nt=100
tsc=2.494102
np=40
host=cas114
kernel=tl_f90_cg_calc_w

#mpicc -O2 -Wall -o ${kernel}_tsc_tf.x ./${kernel}_tf.c  -DTIMING -DUSE_TSC -DNTEST=$nt -DNPASS=1 -lgsl -lopenblas
#mpicc -O2 -Wall -o ${kernel}_cgt_tf.x ./${kernel}_tf.c  -DTIMING -DUSE_CGT -DNTEST=$nt -DNPASS=1 -lgsl -lopenblas
#mpicc -O2 -Wall -o ${kernel}_papi_tf.x ./${kernel}_tf.c  -DTIMING -DUSE_PAPI -DNTEST=$nt -DNPASS=1 -lgsl -lopenblas -lpapi                               
#mpicc -O2 -Wall -o ${kernel}_papix6_tf.x ./${kernel}_tf.c  -DTIMING -DUSE_PAPIX6 -DNTEST=$nt -DNPASS=1 -lgsl -lopenblas -lpapi                           
#mpicc -O2 -Wall -o ${kernel}_likwid_tf.x ./${kernel}_tf.c  -DTIMING -DUSE_LIKWID -DLIKWID_PERFMON -DNTEST=$nt -DNPASS=1 -lgsl -lopenblas -llikwid

mpicc -O2 -Wall -o ${kernel}_tsc.x ./${kernel}.c  -DTIMING -DUSE_TSC -DNTEST=$nt -DNPASS=1 -lgsl -lopenblas
mpicc -O2 -Wall -o ${kernel}_cgt.x ./${kernel}.c  -DTIMING -DUSE_CGT -DNTEST=$nt -DNPASS=1 -lgsl -lopenblas
mpicc -O2 -Wall -o ${kernel}_papi.x ./${kernel}.c  -DTIMING -DUSE_PAPI -DNTEST=$nt -DNPASS=1 -lgsl -lopenblas -lpapi
mpicc -O2 -Wall -o ${kernel}_papix6.x ./${kernel}.c  -DTIMING -DUSE_PAPIX6 -DNTEST=$nt -DNPASS=1 -lgsl -lopenblas -lpapi                        
mpicc -O2 -Wall -o ${kernel}_likwid.x ./${kernel}.c  -DTIMING -DUSE_LIKWID -DLIKWID_PERFMON -DNTEST=$nt -DNPASS=1 -lgsl -lopenblas -llikwid  

for m in tsc cgt papi papix6
do
    mpirun --map-by core --bind-to core -np ${np} ./${kernel}_${m}.x $narr $tsc
    rm  -r ${kernel}${narr}n${nt}t_${m}_${host}
    mkdir ${kernel}${narr}n${nt}t_${m}_${host}
    mv ./*.csv ${kernel}${narr}n${nt}t_${m}_${host}
done

likwid-mpirun -np ${np} -g L3 -m ./${kernel}_likwid.x $narr $tsc
rm  -r ${kernel}${narr}n${nt}t_likwid_${host}
mkdir ${kernel}${narr}n${nt}t_likwid_${host}
mv ./*.csv ${kernel}${narr}n${nt}t_likwid_${host}

