import numpy as np
import os
import re

######## README ##############################################
#
#  This code does the following
#
#       - Computes JOREK profiles via postproc at selected times
#       - Generates a gnuplot script to visualize them
#
#  Specify the input parameters requested below
#  and run it with the following command
#
#         python plot_profiles.py && gnuplot -p gnu_script 
#
##############################################################


######## Input parameters ####################################

# Add or remove the profiles that you want (the names must
# match the jorek2_postproc names, except q-profile)
profs = ["ne", "nimp", "Te", "currdens", "Psi_N"]

# Specify the requested steps to plot 
steps = [3000, 3180, 3600, 4160, 4880]

# Provide location of postproc executable and JOREK input file
postproc = "./jorek2_postproc"
input    = "input_fake"
midplane = True  # Put to true to get midplane profiles instead of average

# Put to 1 to compute the profiles and to 0 just to tune the plotting
compute  = 1
terminal = 'qt'  #or x11

#################################################################

# Internal parameters
ntimes  = len(steps)
if midplane: 
    out_name = "midplane"
else:
    out_name = "average"
plot_folder  = "toplot"


if (compute == 1):
    # create postproc input file
    post_inp = "pinp_profiles"

    # expression list
    expr = "expressions t Psi_N" 
    if midplane:
        expr = "expressions t R" 

    for prf in profs:
        expr = expr + " " + prf

    # Create directory to save the data
    try:
        os.system("rm -r "+plot_folder)
    except:
        pass
    
    os.mkdir(plot_folder)

    for i in range(0,ntimes):
        f2 = open(post_inp, "w")
        f2.write( "namelist "+ input +"\n")
        f2.write( "si-units" +"\n")
        f2.write( "for step "+str(steps[i])+ " do" +"\n")
        f2.write( "set surfaces 200" +"\n")
        f2.write( "qprofile" +"\n")
        f2.write( expr +"\n")
        f2.write( "mark_coords 1"+"\n")
        if midplane:
            f2.write( "midplane"+"\n")
        else:
            f2.write( "average"+"\n")
        f2.write( "done"+"\n" )
        f2.close()  

        if midplane: 
            os.system("rm postproc/exprs_midplane_s*")
        else: 
            os.system("rm postproc/exprs_averaged_s*")
            os.system("rm postproc/qprofile_s*")

        # Run postproc
        command = postproc +" < "+ post_inp 
        os.system(command)

        # Copy generated files to plotting folder
        if midplane: 
            os.system("cp "+ "postproc/exprs_midplane_s*" +" toplot/midplane"+str((i+1))+".dat"  )
        else:
            os.system("cp "+ "postproc/exprs_averaged_s*" +" toplot/average"+str((i+1))+".dat"  )
            os.system("cp "+ "postproc/qprofile_s*"+" toplot/qprof"+str((i+1))+".dat"  )

nplot = len(profs)

class profile:
    label="ne (1e20 m^-3)"
    postproc="ne"
    fact=1
    log=False

# Profile settings (add of modify profiles here)
prof_list = []
for prof in profs:
    if (prof=="ne"): 
        density = profile()
        density.label    = "n_e (1e20 m^-3)"
        density.fact     = 1e-20
        density.postproc = "ne"
        density.log      = False
        prof_list.append(density)

    elif (prof=="Te"): 
        temperature = profile()
        temperature.label    = "T_e (eV)"
        temperature.fact     = 1
        temperature.postproc = "Te"
        temperature.log      = True
        prof_list.append(temperature)

    elif (prof=="currdens"): 
        curr = profile()
        curr.label    = "J_{\phi} (MA/m^2)"
        curr.fact     = -1e-6
        curr.postproc = "currdens"
        curr.log      = False
        prof_list.append(curr)

    elif (prof=="nimp"): 
        densimp = profile()
        densimp.label    = "n_{imp} (1e20 m^-3)"
        densimp.fact     = 1e-20
        densimp.postproc = "nimp"
        densimp.log      = False
        prof_list.append(densimp)

    elif (prof=="Psi_N"): 
        densimp = profile()
        densimp.label    = "Psi_N"
        densimp.fact     = 1
        densimp.postproc = "Psi_N"
        densimp.log      = False
        prof_list.append(densimp)

    elif (prof=="q-profile"): 
        qprof = profile()
        qprof.label    = "q-profile"
        qprof.fact     = 1
        qprof.postproc = "None"
        qprof.log      = False
        prof_list.append(qprof)

    else:
        print("invalid profile!")
        break

# Generate gnuplot file
f = open("gnu_script", "w")
f.write("reset \n")
f.write( "set terminal "+terminal+ "  size 600,900\n")  
f.write( "set lmargin at screen 0.25\n")  
f.write( "set multiplot\n")  
f.write( "unset key\n")  
f.write( "set linetype 1 lc 'black'       lw 2\n") 
f.write( "set linetype 2 lc 'grey'        lw 2\n") 
f.write( "set linetype 3 lc rgb '#800000' lw 2\n") 
f.write( "set linetype 4 lc rgb '#ff0000' lw 2\n") 
f.write( "set linetype 5 lc rgb '#ff4500' lw 2\n") 
f.write( "set linetype 6 lc rgb '#ffa500' lw 2\n") 
f.write( "set linetype 7 lc rgb '#006400' lw 2\n") 
f.write( "set linetype 8 lc rgb '#0000ff' lw 2\n") 
f.write( "set linetype 9 lc rgb '#9400d3' lw 2\n") 
f.write( "set linetype 10 lc 'pink'       lw 2\n") 
f.write( "set grid\n") 
f.write("\n")

xsize = 0.7
ysize = 1.0/ float(nplot)

f.write("ntimes=" + str(ntimes) + "\n")

for i in range(0,nplot):

    f.write("set size " +str(xsize)+", " + str(ysize) + "\n")
    f.write("set origin 0.,"+str(1-ysize*(i+1))+ "\n")
    f.write("set format x ''" + "\n")
    f.write("unset logscale" + "\n")
    f.write("set ylabel '"+prof_list[i].label+"'\n")

    if (prof_list[i].log):
        f.write("set logscale y \n")

    if ( i == nplot-1):
        f.write("unset format x" + "\n")
        f.write("set xtics" + "\n")
        if (not midplane): 
            f.write("set xlabel 'Psi_N'" + "\n")
        else:
            f.write("set xlabel 'R (m)'" + "\n")

    if (prof_list[i].label == 'q-profile'):        
        f.write("plot for [i=1:ntimes] 'toplot/qprof'.i.'.dat' u 2:($"+str(3)+"*"+str(prof_list[i].fact)+") w l\n")
    else:
        if midplane: 
            f.write("plot for [i=1:ntimes] 'toplot/midplane'.i.'.dat' u 2:($"+str(i+3)+"*"+str(prof_list[i].fact)+") w l\n")
        else:
            f.write("plot for [i=1:ntimes] 'toplot/average'.i.'.dat' u 2:($"+str(i+3)+"*"+str(prof_list[i].fact)+") w l\n")

    f.write("\n")


# plot legends
f.write("set key right center" + "\n")
f.write("set origin 0.7,0.4" + "\n")
f.write("set size " +str(1-xsize)+", " + str(ysize) + "\n")
f.write("set border 0" + "\n")
f.write("set yrange [-1e40:-1e39]" + "\n")
f.write("unset tics" + "\n")
f.write("unset xlabel" + "\n")
f.write("unset ylabel" + "\n")

times = np.zeros(ntimes)
for i in range(ntimes):
    try:
        if midplane: 
            dat  = np.loadtxt('toplot/midplane'+str(i+1)+'.dat')
        else:
            dat  = np.loadtxt('toplot/average'+str(i+1)+'.dat')
        times[i] = dat[0,0]*1e3
    except: 
        times[i] = 0.0

legends = "p 2 t '" + "{:.2f}".format(times[0]) + " ms'"
for i in range(1,ntimes):
    legends = legends + ", 2 t '" + "{:.2f}".format(times[i]) + " ms'"      
f.write(legends + "\n")

f.write( "unset multiplot\n") 
f.close()
