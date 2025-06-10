import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

def add_trace(lines,x,y,expression, ix, leg, colr, ltype):
    ynew = []
    for i in range(0,20):
        expression=expression.replace("x"+str(i)+" ", "x["+str(i)+"][i]")
        expression=expression.replace("y"+str(i)+" ", "y["+str(i)+"][i]")
    
    for i in range(0,len(x[ix])):
        ynew.append(eval(expression))
    x.append(x[ix])
    y.append(ynew)
    lines.append(line("new","unknown", 0, 1, leg, colr, ltype,1,1,0,0))
    
def load_data(lines, x, y):
    Nplot=len(lines);  data=[];
    
    for i in range (0, Nplot):
        data.append( np.loadtxt(lines[i].file_ad) ) 
        x.append(data[i][:,lines[i].line_colx]*lines[i].line_fx-lines[i].line_sx) 
        y.append(data[i][:,lines[i].line_coly]*lines[i].line_fy-lines[i].line_sy) 
    
class line:
  def __init__(self, line_name, file_ad,  line_colx, line_coly, line_leg, 
               line_color, line_mark, line_fx, line_fy, line_sx, line_sy):
    self.file_ad   = file_ad;     self.line_leg = line_leg;   self.line_color  = line_color;
    self.line_colx = line_colx;  self.line_coly = line_coly;    self.line_mark = line_mark; 
    self.line_name = line_name;    self.line_fx = line_fx;        self.line_fy = line_fy;
    self.line_sx = line_sx;        self.line_sy = line_sy;


class axe2:
  def __init__(self, curves, ylabel, ylab_font, yr_log, yrange,yticks_log, yticks, leglog, legfont, legloc, logarit):
    self.curves = curves;  self.ylabel = ylabel;   self.ylab_font  = ylab_font;
    self.yr_log = yr_log;  self.yrange = yrange;   self.yticks_log = yticks_log;  
    self.yticks = yticks;  self.leglog = leglog;   self.legfont    = legfont;       self.legloc= legloc;
    self.logarit = logarit


def do_plot_raw4(x,y,scale, frame_ratio, x_label, x_lab_font, x_range,
            grid_log, grid_thick,legend_log, legend_font, legend_loc, leg_lab, leg_col, ticks_font, x_range_log,
            ticks_x_log, ticks_x, text_log, text_loc, text, text_font,
            sup_tit, sup_tit_lg, sup_tit_font, sub_plts, xlog_scale, tight_save, lines, save_plot, p_name):
    
    Nsub = len(sub_plts)
    fig, axes_tmp = plt.subplots(Nsub, 1, sharex=True,figsize=(frame_ratio*5*scale,5*scale))
    if (Nsub==1):   # Fix indexing problems with 1 subplot
        axes = [0]; axes[0]=axes_tmp
    else:
        axes = axes_tmp
    plt.subplots_adjust(hspace=.0)

    for i in range(0,Nsub):
       
        Nplot = len(sub_plts[i].curves);
        
        for j in range (0, Nplot):
            ic = sub_plts[i].curves[j]
            lc = lines[ic].line_color;     ll = lines[ic].line_leg;    lt=lines[ic].line_mark;
            axes[i].plot(x[ic],y[ic],lt,label=ll,color=lc)  # replace plot by semilogx, semilogy or log 
                                                                  
        if (x_range_log):
            axes[i].set_xlim([x_range[0],x_range[1]])                                   

        if (sub_plts[i].logarit):
            axes[i].set_yscale('log')                             

        if (sub_plts[i].yr_log):
            axes[i].set_ylim([sub_plts[i].yrange[0],sub_plts[i].yrange[1]])                                   
    
        axes[i].set(xlabel=x_label, ylabel=sub_plts[i].ylabel) 
        axes[i].xaxis.label.set_size(x_lab_font*scale)
        axes[i].yaxis.label.set_size(sub_plts[i].ylab_font*scale)
              
        if (xlog_scale):
            axes[i].set_xscale('log')

        axes[i].tick_params(axis="both", which="both", labelsize=ticks_font*scale)   # Tick size and format                                           
             
        if (ticks_x_log):
            axes[i].set_xticks(np.arange(ticks_x[0], ticks_x[1], ticks_x[2]))

        if (sub_plts[i].yticks_log):
            axes[i].set_yticks(np.arange(sub_plts[i].yticks[0], sub_plts[i].yticks[1], sub_plts[i].yticks[2]))

        if (grid_log):
            axes[i].grid(visible=True, which="both", axis="both", linewidth=grid_thick)                                               

        if (sub_plts[i].leglog):
            axes[i].legend(title=leg_lab, fontsize=sub_plts[i].legfont*scale, ncol=leg_col, bbox_to_anchor=(sub_plts[i].legloc[0], sub_plts[i].legloc[1]))             # Legend fontsize and location (0 for best)                              
            axes[i].get_legend().get_title().set_fontsize(legend_font*scale) 

    if (text_log):
        axes[0].text(text_loc[0], text_loc[1], text, fontsize=text_font*scale, transform=plt.gcf().transFigure)    # Annotate something in the text            

    if (sup_tit_lg):    
        fig.suptitle(sup_tit, fontsize=sup_tit_font)

    plt.show() 
    
    if (save_plot):
        if (tight_save):
            fig.savefig(p_name, bbox_inches='tight')
        else:
            fig.savefig(p_name)


