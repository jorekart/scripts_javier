#!/bin/bash

# Instructions
#   1. Set all input params below correctly
#   2. Compute everything and plot (./plot_seps.sh -all)
#   3. Just plot (./plot_seps.sh )

# Input params 
postproc="../jorek2_postproc_norho_n13_Fnoll_dZ"
sep_times_ex="/home/ITER/artolaj/scripts_hub/scripts_javier/set_sep_times.py"
n_seps=10             # number of separatrices
input="input"
zeroDfile="get_traces/0D.dat"


# Code

python $sep_times_ex $zeroDfile $n_seps 

# Read times to plot separatrices
declare -a tim=()
declare -a steps=()

readarray -t dat < ./times_restarts.txt
for line in "${dat[@]}"
do
  time_restart=($line)
  tim+=(${time_restart[0]})
  steps+=(${time_restart[1]})
done

len=${#steps[@]}
len2=$((len-1))


####### Call postproc ###############################

if [ "$1" == "-all" ]; then
  # create postproc input file
  for i in $(seq 0 $len2)
  do
    echo "namelist ${input}" >> pinp_sep
    echo "for step ${steps[i]} do" >> pinp_sep
    echo "separatrix" >> pinp_sep
    echo "done" >> pinp_sep
    $postproc < pinp_sep
    rm pinp_sep
  done
  
  for i in $(seq 0 $len2)
  do
     SEPfile="postproc/separatrix_s${steps[i]}.dat"
     cp $SEPfile "toplot2/sep$((i+1)).dat"
  done
fi

rm gnu_script2

# Do gnuplot script
echo "reset" >> gnu_script2 
echo "set term qt" >> gnu_script2 
echo " "
#echo "set size 1,1" >> gnu_script2 
echo "set linetype 1 lc rgb '#800000' lw 2" >> gnu_script2 
echo "set linetype 2 lc rgb '#ff0000' lw 2" >> gnu_script2 
echo "set linetype 3 lc rgb '#ff4500' lw 2" >> gnu_script2 
echo "set linetype 4 lc rgb '#ffa500' lw 2" >> gnu_script2 
echo "set linetype 5 lc rgb '#006400' lw 2" >> gnu_script2 
echo "set linetype 6 lc rgb '#0000ff' lw 2" >> gnu_script2 
echo "set linetype 7 lc rgb '#9400d3' lw 2" >> gnu_script2 
echo "set grid" >> gnu_script2 
echo "" >> gnu_script2 
echo "nplot=${len}" >> gnu_script2 
echo "" >> gnu_script2 

echo "">> gnu_script2 
echo "# separatrix ">> gnu_script2 
echo "">> gnu_script2 
echo "set size ratio -1">> gnu_script2 
echo "set key outside">> gnu_script2 
echo "set xtics">> gnu_script2 
echo "set ylabel 'Z (m)' font 'Times-Roman,$label_font'">> gnu_script2 
echo "set xlabel 'R (m)' font 'Times-Roman,$label_font'">> gnu_script2 

echo -n "p 'toplot2/sep1.dat' u 1:2 w l t 't=${tim[0]} ms', ">> gnu_script2 
for i in $(seq 1 $len2)
do
echo -n "'toplot2/sep$((i+1)).dat' u 1:2 w l t 't=${tim[i]} ms',  ">> gnu_script2 
done



echo "">> gnu_script2 
echo "">> gnu_script2 

gnuplot -p -c gnu_script2 

