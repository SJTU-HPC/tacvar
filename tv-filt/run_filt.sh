#!/bin/bash -x 
nt=100
np=40
host=cas114
filt_dir=~/03-Project/_Research_Perf_var/10-Accuracy/codes/FilTAndVKern/filt
#kernel=tl_f90_cg_calc_w
for kernel in jacobi2d5p tl_f90_cg_calc_w
do
    for m in tsc cgt papi papix6
    do
        for iarr in 60 64 100 120 128 200 250 256 400 510 512 800 1000 1020 1024 2000 2048 3000 3072 4000 4096
        do
            res_dir=~/03-Project/_Research_Perf_var/10-Accuracy/experiments/ipdps-stencil/${kernel}${iarr}n${nt}t_${m}_${host}_filt
            cd $filt_dir
            cp ./filt.x $res_dir/filt.x
            cd $res_dir
            echo 1.0 > er.out
            er0=`cat er.out`
            for iw in 10 50 100 150 200
            do
                ./filt.x -w $iw -n 100000 -l 0.01 -z 0.01
                er=`cat er.out`
                result=$(echo "$er - $er0" | bc)
                if (( $(echo "$result < 0" | bc -l) )); then
                    mv sim_cdf.csv sim_cdf_0.csv
                    mv tr_hist.csv tr_hist_0.csv
                    mv er.out er_0.out
                    mv ep.out ep_0.out
                    mv wd.out wd_0.out
                    er0=$er
                fi
            done
            mv sim_cdf_0.csv sim_cdf.csv
            mv tr_hist_0.csv tr_hist.csv
            mv er_0.out er.out
            mv ep_0.out ep.out
            mv wd_0.out wd.out
            rm filt.x
        done
    done
done
