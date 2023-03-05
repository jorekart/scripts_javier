import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from plot_raw import * 

 
lines = [];           # Lines array (don't touch)      

zeroD1   = 'get_traces/0D.dat'
mag_ener = 'get_traces/magnetic_energies.dat'

ncolor = 1
colors = plt.cm.jet(np.linspace(0,1,ncolor))

c0 = "red";  c1 = "blue";  c2="black"; c3 = "green"

class qtty:
  def __init__(self, name, column, scale, shift): 
    self.name = name;  self.column = column; self.scale=scale;  self.shift=shift;

Ip   = qtty(r"$I_p$ [MA]", 62, 1e-6, 0)  # legend, column, scaling factor, shift
Ih_p = qtty(r"$I_{halo,\theta}$ [MA]", 75, 1, 0)
Zc   = qtty(r"$Z_{curr}$ [m]", 6, 1, 0)
q95  = qtty(r"$q_{95}$ ", 73, 1, 0)
q0   = qtty(r"$q_{0}$ ", 72, 1, 0)
TPF  = qtty(r"TPF ", 76, 1, 0)
Fnoll= qtty(r"$F_{Noll}$ [MN] ", 77, 1, 0)
Wth  = qtty(r"$W_{th}$ [MJ] ", 16, 1e-6, 0)
Wm_n0= qtty(r"$n=0$ ", 1, 1, 0)
Wm_n1= qtty(r"$n=1$ ", 2, 1, 0)
Wm_n2= qtty(r"$n=2$ ", 3, 1, 0)
Wm_n3= qtty(r"$n=3$ ", 4, 1, 0)
 
tscale = 1.e3;  tshift = 0.0
tscale2= 1
p1 = []; p2 = []; p3 = []; p4 = []; p5 = []; 
#############################################################################
### Create here the lines and specify their properties (color,type, etc)  ###
#############################################################################

# Ip and halo    
lines.append( line("", zeroD1, 0, Ip.column,     Ip.name,  c0, "-",  tscale,  Ip.scale,       Ip.shift,  tshift) );   p1.append( len(lines)-1 )
lines.append( line("", zeroD1, 0, Ih_p.column, Ih_p.name,  c1, "-",  tscale,  Ih_p.scale,   Ih_p.shift,  tshift) );   p1.append( len(lines)-1 )

# Zcurr
lines.append( line("", zeroD1, 0, Zc.column,     Zc.name,  c0, "-",  tscale,  Zc.scale,       Zc.shift,  tshift) );   p2.append( len(lines)-1 ) 

# q
lines.append( line("", zeroD1, 0, q95.column,     q95.name,  c1, "-",  tscale,  q95.scale,       q95.shift,  tshift) );   p2.append( len(lines)-1 ) 
#lines.append( line("", zeroD1, 0,  q0.column,      q0.name,  c1, "--",  tscale,   q0.scale,       q0.shift,  tshift) );   p2.append( len(lines)-1 ) 

# TPF
lines.append( line("", zeroD1, 0,  TPF.column,      TPF.name,  c2, "-",  tscale,   TPF.scale,       TPF.shift,  tshift) );   p2.append( len(lines)-1 ) 

# F_Noll
lines.append( line("", zeroD1, 0, Fnoll.column,   Fnoll.name,  c1, "-",  tscale,  Fnoll.scale,      Fnoll.shift,  tshift) );   p3.append( len(lines)-1 ) 

# W_th
lines.append( line("", zeroD1, 0, Wth.column,   Wth.name,  c1, "-",  tscale,  Wth.scale,      Wth.shift,  tshift) );   p4.append( len(lines)-1 ) 

#Wmag
lines.append( line("", mag_ener, 0, Wm_n0.column,   Wm_n0.name,  c0, "-",  tscale2,  Wm_n0.scale,      Wm_n0.shift,  tshift) );   p5.append( len(lines)-1 ) 
lines.append( line("", mag_ener, 0, Wm_n1.column,   Wm_n1.name,  c1, "-",  tscale2,  Wm_n1.scale,      Wm_n1.shift,  tshift) );   p5.append( len(lines)-1 ) 
lines.append( line("", mag_ener, 0, Wm_n2.column,   Wm_n2.name,  c2, "-",  tscale2,  Wm_n2.scale,      Wm_n2.shift,  tshift) );   p5.append( len(lines)-1 ) 
lines.append( line("", mag_ener, 0, Wm_n3.column,   Wm_n3.name,  c3, "-",  tscale2,  Wm_n3.scale,      Wm_n3.shift,  tshift) );   p5.append( len(lines)-1 ) 

# Load and create new lines (with add trace)  
x = [];  y = []; sub_plts   = []; load_data(lines, x, y);  # (look but don't touch)    
    
#############################################################################
###                      Mother plot parameters                           ###      
#############################################################################    
    
scale       = 2.0;                 frame_ratio =0.7;                    
x_label     = r"Time (ms)";       x_lab_font = 10;                      
x_range_log = True;             x_range     = [x[0][0], x[0][-1]*1.5];     
grid_log    = True;              grid_thick  = 1;     
legend_log  = False;             legend_font = 18;             legend_loc   = [0.2,0.5];      
leg_lab     = "";                leg_col    = 1;    
ticks_font  = 7;                 ticks_x_log = True;          ticks_x      = [0,60,5];                
text_log    = False;             text_loc    = [0.25,0.5];     text         = r"";      text_font = 12;    
sup_tit_lg  = False;             sup_tit     = r"Wall density scan, $\gamma_{sheath}=8$";        sup_tit_font = 18;       
save_plot   = True;             p_name      = "overview.pdf";     
tight_save  = True;
xlog_scale  = False;
#############################################################################
###                      Create subplots                                  ###      
#############################################################################

#           curves selection,      ylabel, ylab_font, yr_log,    yrange,  yticks_log,        yticks, leglog, legfont,      legloc    
sub_plts.append( axe2( p1 , "",         8,  False,     [0,1],       False, [0.,10,1],   True,   7, [0.65,0.7], False  ))    
sub_plts.append( axe2( p2 , "",   8,  False,     [0,1],       True, [0.,4,0.5],   True,       7, [0.7,0.7], False  ))    
sub_plts.append( axe2( p3 , "",   8,  False,     [0,1],       False, [0.,10,1],   True,       7, [0.6,0.5], False  ))    
sub_plts.append( axe2( p4 , "",   8,  False,     [0,1],       False, [0.,10,1],   True,       7, [0.65,0.5], False  ))    
sub_plts.append( axe2( p5 , r"$W_{mag}^n$",   8,  False,     [0,1],       False, [0.,10,1],   True,       7, [0.7,0.9], True  ))    


do_plot_raw4(x,y,scale, frame_ratio, x_label, x_lab_font,      
            x_range, grid_log, grid_thick, legend_log, legend_font, legend_loc, leg_lab, leg_col, ticks_font,    
            x_range_log, ticks_x_log, ticks_x, text_log, text_loc,     
            text, text_font, sup_tit, sup_tit_lg, sup_tit_font, sub_plts, xlog_scale, tight_save, lines,save_plot, p_name) 

