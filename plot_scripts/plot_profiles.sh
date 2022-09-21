#!/bin/bash

# Purpose computes and plots JOREK profiles
# Usage, use -c option to re-compute the data with postrpoc

# Steps to plot
steps=(00180 00260 00340 00440 00540 00640 00740)
tim=(3.01 6.09 11.0 17.2 23.3 29.4 35.4)
#tim=("${steps[@]}")
postproc="../jorek/jorek2_postproc"
input="input_CQ"


len=${#steps[@]}
len2=$((len-1))

if [ "$1" == "-c" ]; then
  # create postproc input file
  for i in $(seq 0 $len2)
  do
    echo "namelist ${input}" >> pinp_profiles
    echo "si-units" >> pinp_profiles
    echo "for step ${steps[i]} do" >> pinp_profiles
    echo "set surfaces 200" >> pinp_profiles
    echo "qprofile" >> pinp_profiles
    echo "separatrix" >> pinp_profiles
    echo "expressions Psi_N Te ne currdens" >> pinp_profiles
    echo "mark_coords 1">> pinp_profiles
    echo "average">> pinp_profiles
    echo "done" >> pinp_profiles
    $postproc < pinp_profiles
    rm pinp_profiles
  done
  
  rm -r toplot
  mkdir toplot
  
  for i in $(seq 0 $len2)
  do
     AVfile="postproc/exprs_averaged_s${steps[i]}.dat"
     SEPfile="postproc/separatrix_s${steps[i]}.dat"
     qfile="postproc/qprofile_s${steps[i]}.dat"
     cp $AVfile "toplot/average$((i+1)).dat"
     cp $qfile "toplot/qprof$((i+1)).dat"
     cp $SEPfile "toplot/sep$((i+1)).dat"
  done
fi

rm gnu_script2

# Do gnuplot script
echo "reset" >> gnu_script2 
echo " "
echo "set size 1,1" >> gnu_script2 
echo "set multiplot " >> gnu_script2  
echo "set lmargin at screen 0.1" >> gnu_script2 
echo "unset key" >> gnu_script2 
echo "set linetype 1 lc rgb '#800000' lw 2" >> gnu_script2 
echo "set linetype 2 lc rgb '#ff0000' lw 2" >> gnu_script2 
echo "set linetype 3 lc rgb '#ff4500' lw 2" >> gnu_script2 
echo "set linetype 4 lc rgb '#ffa500' lw 2" >> gnu_script2 
echo "set linetype 5 lc rgb '#006400' lw 2" >> gnu_script2 
echo "set linetype 6 lc rgb '#0000ff' lw 2" >> gnu_script2 
echo "set linetype 7 lc rgb '#9400d3' lw 2" >> gnu_script2 
echo "set grid" >> gnu_script2 
echo "" >> gnu_script2 
echo "nplot=7" >> gnu_script2 
echo "" >> gnu_script2 
echo "# ne" >> gnu_script2 
echo "set size 0.5, 0.25" >> gnu_script2 
echo "set origin 0., 0.75" >> gnu_script2 
echo "set format x ''" >> gnu_script2 
echo "set ylabel 'n_e (1e20 m-3)'" >> gnu_script2 
echo "plot for [i=1:nplot] 'toplot/average'.i.'.dat' u 1:(\$3*1e-20) w l" >> gnu_script2 
echo "" >> gnu_script2 
echo "" >> gnu_script2 
echo "# ne" >> gnu_script2 
echo "set size 0.5, 0.25" >> gnu_script2 
echo "set origin 0., 0.54" >> gnu_script2 
echo "set format x ''" >> gnu_script2 
echo "set ylabel 'T_e (eV)'" >> gnu_script2 
echo "plot for [i=1:nplot] 'toplot/average'.i.'.dat' u 1:2 w l" >> gnu_script2 
echo "" >> gnu_script2 
echo "# Jphi " >> gnu_script2 
echo "set size 0.5, 0.25" >> gnu_script2 
echo "set origin 0., 0.33" >> gnu_script2 
echo "set format x ''" >> gnu_script2 
echo "set ylabel 'Jphi (MA/m2)'" >> gnu_script2 
echo "plot for [i=1:nplot] 'toplot/average'.i.'.dat' u 1:(-\$4*1e-6) w l" >> gnu_script2 
echo "" >> gnu_script2 
echo "# q " >> gnu_script2 
echo "set size 0.5, 0.37" >> gnu_script2 
echo "set origin 0., 0.0 " >> gnu_script2 
echo "unset format x" >> gnu_script2 
echo "set xtics" >> gnu_script2 
echo "set ylabel 'q' " >> gnu_script2 
echo "set xlabel 'Psi_N'" >> gnu_script2 
echo "plot for [i=1:nplot] 'toplot/qprof'.i.'.dat' u 1:(-\$2) w l">> gnu_script2 

echo "">> gnu_script2 
echo "# separatrix ">> gnu_script2 
echo "unset lmargin">> gnu_script2 
echo "set rmargin 15">> gnu_script2 
echo "set key at screen 1, graph 1">> gnu_script2 
echo "">> gnu_script2 
echo "set size 0.5, 1.0">> gnu_script2 
echo "set size ratio -1">> gnu_script2 
echo "set origin 0.5, 0.0 ">> gnu_script2 
echo "unset format x">> gnu_script2 
echo "set xtics">> gnu_script2 
echo "set ylabel 'Z (m)'">> gnu_script2 
echo "set xlabel 'R (m)'">> gnu_script2 

echo -n "p 'toplot/sep1.dat' u 1:2 w l t 't=${tim[0]} ms', ">> gnu_script2 
for i in $(seq 1 $len2)
do
echo -n "'toplot/sep$((i+1)).dat' u 1:2 w l t 't=${tim[i]} ms',  ">> gnu_script2 
done



echo "">> gnu_script2 
echo "">> gnu_script2 
echo "unset multiplot">> gnu_script2 

gnuplot -p -c gnu_script2 

