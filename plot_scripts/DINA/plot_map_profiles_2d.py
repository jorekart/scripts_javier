import imas
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.widgets import Slider, TextBox
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, TwoSlopeNorm, LinearSegmentedColormap

# ---------------- Input ----------------
# 100187 15MA up. 100200 15MA down
uri = 'imas:mdsplus?user=public;pulse=100200;run=1;database=ITER_DISRUPTIONS;version=3'
ind0 = 100   # initial time index
MODES = ['qpar', 'jtor', 'te', 'qcalc']   # disruption q_par, core_profiles j_total, electron temperature, q_par from Te (sheath formula)
POSITIVE_MODES = ('te', 'qcalc')   # non-negative quantities: scale 0..max, white-to-red colormap
MODE = 'te'   # starting mode; press 't' in the plot to cycle through MODES

# q_par from the sheath formula of plasma_calculator.py: q_par = e * gamma * ne * Te * c_s(Te, Ti)
NE_CALC       = 1e20      # hardcoded electron density [m^-3]
GAMMA_SHEATH  = 7.0       # sheath heat transmission factor (not defined in plasma_calculator.py; assumed value)
E_CHARGE      = 1.60217662e-19
ION_MASS      = 1.67e-27 * 2.014   # deuterium, as in plasma_calculator.py


def q_par_calc(Te):
    """Te in eV (Ti = Te assumed) -> q_par in W/m^2."""
    Te = np.asarray(Te, dtype=float)
    c_s = ((Te + Te) * E_CHARGE / ION_MASS) ** 0.5
    return E_CHARGE * GAMMA_SHEATH * NE_CALC * Te * c_s


fw = np.loadtxt('/home/ITER/artolaj/ITER_components/1st_wall.txt')

# ---------------- Read IDSs ----------------
entry = imas.DBEntry(uri=uri, mode='r', dd_version='3.41.0')  # match on-disk DD version
core_t = entry.get('core_profiles')   # full read for the time base and the j_total range
times = np.asarray(core_t.time)
_j = [-np.asarray(p.j_total) for p in core_t.profiles_1d if np.asarray(p.j_total).size]   # skip empty slices; sign as plotted in get_frame
CLIM = {}
if _j:
    CLIM['jtor'] = (np.nanmin([np.nanmin(j) for j in _j]), np.nanmax([np.nanmax(j) for j in _j]))
else:
    print('WARNING: core_profiles.profiles_1d[:].j_total is empty in all slices; '
          'the j_total option will not work (run check_core_fields.py to see what is filled)')
_te = [np.asarray(p.electrons.temperature) for p in core_t.profiles_1d if np.asarray(p.electrons.temperature).size]
if _te:
    CLIM['te'] = (0.0, np.nanmax([np.nanmax(x) for x in _te]))
    CLIM['qcalc'] = (0.0, float(q_par_calc(CLIM['te'][1])))   # q_par is monotonic in Te
else:
    print('WARNING: core_profiles electrons.temperature is empty in all slices; the te option will not work')
del core_t, _j, _te

n_t = len(times)
ind0 = min(ind0, n_t - 1)

# Fixed colour scale per quantity (same for all times); q_par range is read lazily on first use
norm = Normalize(vmin=0.0, vmax=1.0)
# One rainbow, symmetric about 0: white at zero, red at both extremes (negative or positive)
_rainbow = ['white', 'blue', 'cyan', 'lime', 'yellow', 'red']
_up = [(0.5 + 0.5 * i / (len(_rainbow) - 1), c) for i, c in enumerate(_rainbow)]    # 0.5 .. 1
_down = [(1.0 - p, c) for p, c in _up][::-1]                                         # 0 .. 0.5
cmap_sym = LinearSegmentedColormap.from_list('white_rainbow_sym', _down[:-1] + _up)
# Positive-only quantities (temperature): white at 0 up to red at the maximum
cmap_pos = LinearSegmentedColormap.from_list('white_rainbow', _rainbow)
cmap = cmap_sym


def get_clim():
    if MODE not in CLIM:
        dis_t = entry.get('disruption')
        qm = np.nanmax([np.nanmax(q) for q in (np.asarray(p.power_density_conductive_losses)
                                               for p in dis_t.profiles_1d) if q.size])
        qn = np.nanmin([np.nanmin(q) for q in (np.asarray(p.power_density_conductive_losses)
                                               for p in dis_t.profiles_1d) if q.size])
        CLIM[MODE] = (float(qn), float(qm))
    lo, hi = CLIM[MODE]
    if MODE in POSITIVE_MODES:
        return 0.0, hi
    a = max(abs(lo), abs(hi))   # symmetric limits so that 0 (white) is the centre of the colour bar
    return -a, a


def get_frame(ind):
    """Return everything needed to draw the slice closest to times[ind] (as in L2G)."""
    t = times[ind]
    interp = imas.ids_defs.CLOSEST_INTERP
    dis  = entry.get_slice('disruption', t, interp)
    eq   = entry.get_slice('equilibrium', t, interp)
    summ = entry.get_slice('summary', t, interp)
    core = entry.get_slice('core_profiles', t, interp)
    print(f'slice times: dis={dis.time[0]} eq={eq.time[0]} core={core.time[0]} (requested {t})')

    p2d = eq.time_slice[0].profiles_2d[0]
    R = np.asarray(p2d.grid.dim1)
    Z = np.asarray(p2d.grid.dim2)
    psi2d = np.transpose(np.asarray(p2d.psi))   # (nZ, nR)

    # Profile psi always comes from the core_profiles slice
    psi1d = np.asarray(core.profiles_1d[0].grid.psi)
    if MODE == 'jtor':
        qpar =  -np.asarray(core.profiles_1d[0].j_total)
    elif MODE == 'te':
        qpar = np.asarray(core.profiles_1d[0].electrons.temperature)
    elif MODE == 'qcalc':
        qpar = q_par_calc(core.profiles_1d[0].electrons.temperature)
    else:
        qpar = np.asarray(dis.profiles_1d[0].power_density_conductive_losses)

    psi_ax   = float(core.profiles_1d[0].grid.psi_magnetic_axis)
    psi_sep_arr = np.asarray(summ.local.separatrix.position.psi)
    if psi_sep_arr.size:
        psi_sep = float(psi_sep_arr[0])
    else:   # summary slice came back empty: use the equilibrium boundary psi instead
        psi_sep = float(eq.time_slice[0].global_quantities.psi_boundary)
        print('summary separatrix psi empty in this slice, using eq psi_boundary =', psi_sep)
    psi_halo = float(dis.global_quantities.psi_halo_boundary[0])

    # np.interp needs increasing x; outside the profile range -> NaN
    order = np.argsort(psi1d)
    q2d = np.interp(psi2d, psi1d[order], qpar[order], left=np.nan, right=np.nan)
    if MODE in ('jtor', 'te', 'qcalc'):
        q2d = np.ma.masked_invalid(q2d)   # core quantities: no masking
    else:
        # Mask the core: everything between psi_axis and psi_sep (and inside)
        outside = np.sign(psi_sep - psi_ax)
        q2d = np.ma.masked_invalid(np.where((psi2d - psi_sep) * outside >= 0, q2d, np.nan))
    return R, Z, psi2d, q2d, psi_sep, psi_halo, psi_ax


# ---------------- Figure ----------------
fig, ax = plt.subplots(figsize=(7, 9))
plt.subplots_adjust(bottom=0.2)

sm = ScalarMappable(norm=norm, cmap=cmap)
CUSTOM = {}          # per-mode (vmin, vmax) set with the colour sliders
_cache = {}          # last computed frame, so colour-slider moves do not re-read the IDSs
cb = fig.colorbar(sm, ax=ax, pad=0.02)
cb.set_label(r'$q_\parallel$ (power_density_conductive_losses) [W/m$^2$]')


LABELS = {'qpar': r'$q_\parallel$ (power_density_conductive_losses) [W/m$^2$]',
          'jtor': r'$j_\phi$ (core_profiles j_total) [A/m$^2$]',
          'te':   r'$T_e$ (core_profiles electrons.temperature) [eV]',
          'qcalc': r'$q_\parallel$ = e $\gamma$ n$_e$ T$_e$ c$_s$ ($n_e$ = 1e20, $\gamma$ = 7) [W/m$^2$]'}
NAMES = {'qpar': r'$q_\parallel$', 'jtor': r'$j_\phi$', 'te': r'$T_e$', 'qcalc': r'$q_\parallel$ (calc)'}


def draw(ind):
    ind = int(ind)
    key = (ind, MODE)
    if _cache.get('key') != key:
        _cache['key'], _cache['frame'] = key, get_frame(ind)
    R, Z, psi2d, q2d, psi_sep, psi_halo, psi_ax = _cache['frame']
    view = (ax.get_xlim(), ax.get_ylim()) if ax.has_data() else None   # keep the user's zoom
    ax.clear()

    qmin, qmax = CUSTOM.get(MODE) or get_clim()
    if qmin < 0.0 < qmax:   # white stays at 0 whatever the limits
        cur_norm, cmap = TwoSlopeNorm(vcenter=0.0, vmin=qmin, vmax=qmax), cmap_sym
    else:
        cur_norm, cmap = Normalize(vmin=qmin, vmax=qmax), cmap_pos
    sm.set_norm(cur_norm)
    sm.set_cmap(cmap)
    cb.update_normal(sm)
    levels = np.linspace(qmin, qmax, 61)
    cb.set_label(LABELS[MODE])

    # q_par colormap
    ax.contourf(R, Z, q2d, levels=levels, cmap=cmap, norm=cur_norm, extend='both', zorder=1)
    # psi contours
    ax.contour(R, Z, psi2d, levels=40, colors='0.6', linewidths=0.5, zorder=2)
    # separatrix / LCFS and halo boundary
    ax.contour(R, Z, psi2d, levels=[psi_sep],  colors='black', linewidths=2.0, zorder=3)
    ax.contour(R, Z, psi2d, levels=[psi_halo], colors='lime', linewidths=2.0,
               linestyles='--', zorder=3)
    # walls
    ax.plot(fw[:, 0], fw[:, 1], 'k.-', lw=2.5, ms=6)
    # magnetic axis check: contour of the 2D psi at psi_ax (collapses to a point if consistent)
    # and the location of the 2D psi extremum closest to psi_ax
    ax.contour(R, Z, psi2d, levels=[psi_ax], colors='magenta', linewidths=1.5, zorder=3)
    iz, ir = np.unravel_index(np.argmin(np.abs(psi2d - psi_ax)), psi2d.shape)
    ax.plot(R[ir], Z[iz], 'm*', ms=12, mec='k', zorder=4)
    print(f'ind={ind} psi_ax={psi_ax:.3f} psi2d min/max={psi2d.min():.3f}/{psi2d.max():.3f} '
          f'psi_sep={psi_sep:.3f} psi_halo={psi_halo:.3f}')

    ax.legend(handles=[
        Line2D([0], [0], color='0.6', lw=1, label=r'$\psi$ contours'),
        Line2D([0], [0], color='black', lw=2, label='Separatrix / LCFS'),
        Line2D([0], [0], color='lime', lw=2, ls='--', label='Halo boundary'),
        Line2D([0], [0], color='magenta', marker='*', mec='k', ms=10, lw=1.5, label=r'$\psi_{axis}$'),
    ], loc='upper right')
    ax.set_aspect('equal', adjustable='box')
    if view is not None:
        ax.set_xlim(*view[0])
        ax.set_ylim(*view[1])
    ax.set_xlabel('R [m]')
    ax.set_ylabel('Z [m]')
    name = NAMES[MODE]
    ax.set_title(f'{name} map, t = {times[ind]:.4f} s  (index {ind})')
    fig.canvas.draw_idle()


# ---------------- Time slider ----------------
slider_ax = fig.add_axes([0.2, 0.03, 0.6, 0.03])
slider = Slider(slider_ax, 'time index', 0, n_t - 1, valinit=ind0, valstep=1)
slider.on_changed(draw)

# ---------------- Colour-scale input boxes ----------------
cmin_box = TextBox(fig.add_axes([0.2, 0.09, 0.2, 0.04]), 'colour min ', initial='')
cmax_box = TextBox(fig.add_axes([0.6, 0.09, 0.2, 0.04]), 'colour max ', initial='')


def sync_clim_boxes():
    """Show the limits of the current MODE in the boxes."""
    vmin, vmax = CUSTOM.get(MODE) or get_clim()
    cmin_box.set_val(f'{vmin:.4g}')
    cmax_box.set_val(f'{vmax:.4g}')


def on_clim(_):
    try:
        vmin, vmax = float(cmin_box.text), float(cmax_box.text)
    except ValueError:
        print('colour limits must be numbers, e.g. 1e7')
        return
    if vmax <= vmin:
        print('colour max must be larger than colour min')
        return
    CUSTOM[MODE] = (vmin, vmax)
    draw(slider.val)


cmin_box.on_submit(on_clim)   # press Enter in a box to apply
cmax_box.on_submit(on_clim)

# 't' cycles q_par / j_total / T_e, arrow keys step by one index
def on_key(event):
    global MODE
    if cmin_box.capturing or cmax_box.capturing:   # typing in a colour box
        return
    if event.key == 't':
        MODE = MODES[(MODES.index(MODE) + 1) % len(MODES)]
        sync_clim_boxes()
        draw(slider.val)
    elif event.key == 'right':
        slider.set_val(min(slider.val + 1, n_t - 1))
    elif event.key == 'left':
        slider.set_val(max(slider.val - 1, 0))
fig.canvas.mpl_connect('key_press_event', on_key)

sync_clim_boxes()
draw(ind0)
plt.show()
