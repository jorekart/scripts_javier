#!/bin/bash

# Purpose computes and plots JOREK time traces
# Use -c0 to re-compute 0D quantities
# Use -cs to re-compute separatrices
# Use -s to plot steps instead to time
# Use -t to obtain restart indexes from equally spaced times
# slow font initialization in some systems, load twice gnu_script2

# Instructions
#   1. Set all input params below correctly
#   2. Compute the 0D files (./plot_traces.sh -c0)
#   3. Set separatrices times (./plot_traces.sh -t)
#   4. Compute separatrices (./plot_traces.sh -cs)
#   5. Plot all (./plot_traces.sh )

# Input params 
postproc="./jorek2_postproc"
sep_times_ex="/home/artolaj/scripts_javier/set_sep_times.py"
n_seps=5              # number of separatrices
input="input_fake"

# Default settings
#zeroDfile="0D.dat"
zeroDfile="postproc/zeroD_quantities_s00000..99999.dat"
Ip_col=65
Ihalo_col=67
Zaxis_col=5
q95_col=76
li_col=17 #66

label_font=14
li_lab="W_{th} (MJ)"  #"li_3"
fact_li=1.e-6

# Calculate indices from space times (from initial 0D file)
if [ "$1" == "-t" ]; then
  python $sep_times_ex $zeroDfile $n_seps 
fi

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
if [ "$1" == "-c0" ]; then

   echo "namelist ${input}" >> pinp_sep
   echo "si-units" >> pinp_sep
   echo "for step 0 to 99999 by 100 do" >> pinp_sep
   echo "zeroD_quantities" >> pinp_sep
   echo "done" >> pinp_sep
   $postproc < pinp_sep
   rm pinp_sep
   
   rm -r toplot2
   mkdir toplot2
   
   for i in $(seq 0 $len2)
   do
      SEPfile="postproc/separatrix_s${steps[i]}.dat"
      cp $SEPfile "toplot2/sep$((i+1)).dat"
   done
   
fi

if [ "$1" == "-cs" ]; then

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

cp $zeroDfile "toplot2/0D.dat" 

################# PLOT #############################


rm gnu_script2

# Do gnuplot script
echo "reset" >> gnu_script2 
echo "set term qt" >> gnu_script2 
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
echo "nplot=${len}" >> gnu_script2 
echo "" >> gnu_script2 
echo "# Ip" >> gnu_script2 
echo "set size 0.5, 0.25" >> gnu_script2 
echo "set origin 0., 0.75" >> gnu_script2 
echo "set format x ''" >> gnu_script2 
echo "set ylabel 'I (MA)' font 'Times-Roman,$label_font'" >> gnu_script2
echo "set key" >> gnu_script2

if [ "$1" == "-s" ]; then
  echo "p 'toplot2/0D.dat' u 2:(\$ $Ip_col*1e-6) lc 'black' t 'Ip', '' u 2:(\$ $Ihalo_col*1e-6) lc 'red' t 'Ihalo'" >> gnu_script2 
else
  echo "p 'toplot2/0D.dat' u (\$ 1*1e3):(\$ $Ip_col*1e-6) w l lc 'black' t 'Ip', '' u (\$ 1*1e3):(\$ $Ihalo_col*1e-6) w l lc 'red' t 'Ihalo'" >> gnu_script2 
fi

echo "" >> gnu_script2 
echo "" >> gnu_script2 
echo "# Zaxis" >> gnu_script2 
echo "set size 0.5, 0.25" >> gnu_script2 
echo "unset key" >> gnu_script2
echo "set origin 0., 0.54" >> gnu_script2 
echo "set format x ''" >> gnu_script2 
echo "set ylabel 'Z_{axis} (m)' font 'Times-Roman,$label_font'" >> gnu_script2 

if [ "$1" == "-s" ]; then
  echo "p 'toplot2/0D.dat' u 2:(\$ $Zaxis_col) lc 'black'" >> gnu_script2 
else
  echo "p 'toplot2/0D.dat' u (\$ 1*1e3):(\$ $Zaxis_col) w l lc 'black'" >> gnu_script2 
fi

echo "" >> gnu_script2 
echo "# q95 " >> gnu_script2 
echo "set size 0.5, 0.25" >> gnu_script2 
echo "set origin 0., 0.33" >> gnu_script2 
echo "set format x ''" >> gnu_script2 
echo "set ylabel 'q_{95}' font 'Times-Roman,$label_font' " >> gnu_script2 

if [ "$1" == "-s" ]; then
  echo "p 'toplot2/0D.dat' u 2:(\$ $q95_col) lc 'black'" >> gnu_script2 
else
  echo "p 'toplot2/0D.dat' u (\$ 1*1e3):(\$ $q95_col) w l lc 'black'" >> gnu_script2 
fi


echo "" >> gnu_script2 
echo "# li3 " >> gnu_script2 
echo "set size 0.5, 0.37" >> gnu_script2 
echo "set origin 0., 0.0 " >> gnu_script2 
echo "unset format x" >> gnu_script2 
echo "set xtics" >> gnu_script2 
echo "set ylabel '$li_lab' font 'Times-Roman,$label_font'" >> gnu_script2 

if [ "$1" == "-s" ]; then
  echo "set xlabel 'Step' font 'Times-Roman,$label_font'" >> gnu_script2 
  echo "p 'toplot2/0D.dat' u 2:(\$ $li_col * $fact_li) lc 'black'" >> gnu_script2 
else
  echo "set xlabel 'Time (ms)' font 'Times-Roman,$label_font'" >> gnu_script2 
  echo "p 'toplot2/0D.dat' u (\$ 1*1e3):(\$ $li_col * $fact_li) w l lc 'black'" >> gnu_script2 
fi

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
echo "set ylabel 'Z (m)' font 'Times-Roman,$label_font'">> gnu_script2 
echo "set xlabel 'R (m)' font 'Times-Roman,$label_font'">> gnu_script2 

echo -n "p 'toplot2/sep1.dat' u 1:2 w l t 't=${tim[0]} ms', ">> gnu_script2 
for i in $(seq 1 $len2)
do
echo -n "'toplot2/sep$((i+1)).dat' u 1:2 w l t 't=${tim[i]} ms',  ">> gnu_script2 
done



echo "">> gnu_script2 
echo "">> gnu_script2 
echo "unset multiplot">> gnu_script2 

gnuplot -p -c gnu_script2 

