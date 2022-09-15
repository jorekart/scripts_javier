# This program fits an R-Z contour to a given function
# It doesn't work very well, but you can use sliders in matplotlib
# to find the optimal parameters by eye
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import Slider, Button

from scipy.optimize import curve_fit

def RZtht(xx,r0,z0,a0,ellip,tria_u,quad_u,tria_l,quad_l):

    R_bnd = np.zeros(len(xx))
    Z_bnd = np.zeros(len(xx))

    for i, x in enumerate(xx):
        if (x > 0.0 ):
            R_bnd[i] = r0 + a0 * np.cos(x + tria_u*np.sin(x) + quad_u*np.sin(2*x))
        else:
            R_bnd[i] = r0 + a0 * np.cos(x + tria_l*np.sin(x) + quad_l*np.sin(2*x))        
        Z_bnd[i] = z0 + a0 * ellip * np.sin(x)

    return R_bnd, Z_bnd


def model(xx,r0,z0,a0,ellip,tria_u,quad_u,tria_l,quad_l):

    Rmd = RZtht(xx,r0,z0,a0,ellip,tria_u,quad_u,tria_l,quad_l)[0]
    Zmd = RZtht(xx,r0,z0,a0,ellip,tria_u,quad_u,tria_l,quad_l)[1]

    return np.sqrt( (Rmd -  Rmid)**2.0 + (Zmd - Zmid)**2.0 )

# Initial guess
r0 = 2.875;   z0 = 0.2213;   ellip = 1.7264;    a0 = 1.026;   
tria_u = 0.276; quad_u = -0.07;
tria_l = 0.2497; quad_l = 0.1274;  

# This must be close to the center of the boundary that you want
Rmid = 2.87;  Zmid = 0.232;

# Fake model
#x = np.linspace(-np.pi, np.pi,200)
#
#Rmod = 0.; Zmod = 0; e_mod=1; a_mod=1;
#y = model(x,Rmod,Zmod,a_mod,e_mod,0.0,0.0, 0.0, 0.0)
#
#Rtmp = RZtht(x,Rmod,Zmod,a_mod,e_mod,0.0,0.0, 0.0, 0.0)[0]
#Ztmp = RZtht(x,Rmod,Zmod,a_mod,e_mod,0.0,0.0, 0.0, 0.0)[1]
#
#np.savetxt( 'test.txt',     np.transpose( [Rtmp, Ztmp] ) )

# Read coordinates 
data = np.loadtxt('RZ_contour.txt')
#data = np.loadtxt('test.txt')
Rin  = data[:,0]
Zin  = data[:,1]
rmin_in = np.sqrt((Rin - Rmid)**2.0 + (Zin - Zmid)**2.0)
theta   = np.arctan2(Zin-Zmid, Rin-Rmid) 

# Fit
guess      = np.array( [r0,z0,a0,ellip,tria_u,quad_u,tria_l,quad_l])
popt, pcov = curve_fit(model, theta, rmin_in, p0=guess)

print(popt)

x = np.linspace(-np.pi, np.pi,200)
R0f = popt[0]
Z0f = popt[1]
af  = popt[2]
Ef  = popt[3]
TUf = popt[4]
QUf = popt[5]
TLf = popt[6]
QLf = popt[7]


Rf = RZtht(x,R0f,Z0f,af,Ef,TUf,QUf, TLf, QLf)[0]
Zf = RZtht(x,R0f,Z0f,af,Ef,TUf,QUf, TLf, QLf)[1]


# Create the figure and the line that we will manipulate
fig, ax = plt.subplots()
line, = plt.plot(Rf, Zf, lw=2)
ax.plot(Rin, Zin)
ax.set_aspect(aspect='equal') 
ax.set(xlabel='R', ylabel='Z')
ax.grid()


# adjust the main plot to make room for the sliders
plt.subplots_adjust(left=0.05, bottom=0.25)

axQL  = plt.axes([0.25, 0.22, 0.65, 0.02])
QL_sl = Slider(
    ax=axQL,
    label="QL",
    valmin=-0.4,
    valmax=0.4,
    valinit=0.00,
)
axTL  = plt.axes([0.25, 0.19, 0.65, 0.02])
TL_sl = Slider(
    ax=axTL,
    label="TL",
    valmin=-0.4,
    valmax=0.5,
    valinit=0.00,
)

axQU  = plt.axes([0.25, 0.16, 0.65, 0.02])
QU_sl = Slider(
    ax=axQU,
    label="QU",
    valmin=-0.4,
    valmax=0.4,
    valinit=0.00,
)
axR  = plt.axes([0.25, 0.13,0.65, 0.02])
R_sl = Slider(
    ax=axR,
    label='R0',
    valmin=2.5,
    valmax=3.0,
    valinit=2.875,
)

axZ  = plt.axes([0.25, 0.1, 0.65, 0.02])
Z_sl = Slider(
    ax=axZ,
    label='Z0',
    valmin=-0.2,
    valmax=0.5,
    valinit=0.2315,
)

axA  = plt.axes([0.25, 0.07, 0.65, 0.02])
a_sl = Slider(
    ax=axA,
    label='a_min',
    valmin=0.9,
    valmax=1.1,
    valinit=1.026,
)

axE  = plt.axes([0.25, 0.04, 0.65, 0.02])
e_sl = Slider(
    ax=axE,
    label="Ellip",
    valmin=1.4,
    valmax=1.9,
    valinit=1.739,
)

axTU  = plt.axes([0.25, 0.01, 0.65, 0.02])
TU_sl = Slider(
    ax=axTU,
    label="TU",
    valmin=-0.1,
    valmax=0.5,
    valinit=0.26,
)




# The function to be called anytime a slider's value changes
def update(val):
    line.set_xdata( RZtht(x,R_sl.val,Z_sl.val,a_sl.val,e_sl.val,TU_sl.val,QU_sl.val, TL_sl.val, QL_sl.val)[0]  )
    line.set_ydata( RZtht(x,R_sl.val,Z_sl.val,a_sl.val,e_sl.val,TU_sl.val,QU_sl.val, TL_sl.val, QL_sl.val)[1] )
    fig.canvas.draw_idle()


# register the update function with each slider
R_sl.on_changed(update)
Z_sl.on_changed(update)
a_sl.on_changed(update)
e_sl.on_changed(update)
TU_sl.on_changed(update)
QU_sl.on_changed(update)
TL_sl.on_changed(update)
QL_sl.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = plt.axes([0.8, 0.95, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')


def reset(event):
    freq_slider.reset()
    amp_slider.reset()
button.on_clicked(reset)

plt.show()
