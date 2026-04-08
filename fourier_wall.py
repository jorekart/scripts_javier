import numpy as np
from scipy.interpolate import interp1d

vv = np.loadtxt('VV.txt')
outline = vv.copy()

# center
R_geo = np.mean(outline[:, 0])
Z_geo = np.mean(outline[:, 1])

# angle
theta = np.arctan2(outline[:, 1] - Z_geo, outline[:, 0] - R_geo)
theta[theta < 0] += 2*np.pi

# sort by angle
ind = np.argsort(theta)
theta = theta[ind]
outline = outline[ind, :]

# periodic extension for interpolation
theta_ext = np.concatenate(([theta[-1] - 2*np.pi], theta, [theta[0] + 2*np.pi]))
R_ext = np.concatenate(([outline[-1, 0]], outline[:, 0], [outline[0, 0]]))
Z_ext = np.concatenate(([outline[-1, 1]], outline[:, 1], [outline[0, 1]]))

# interpolation
fx = interp1d(theta_ext, R_ext, kind='cubic')
fy = interp1d(theta_ext, Z_ext, kind='cubic')

# uniform theta for FFT
npt = 400
theta_uniform = np.linspace(0.0, 2*np.pi, npt, endpoint=False)

R = fx(theta_uniform)
Z = fy(theta_uniform)

# FFT
coeff_R = np.fft.fft(R) / npt
coeff_Z = np.fft.fft(Z) / npt

coeff_R[1:npt//2] *= 2.0
coeff_Z[1:npt//2] *= 2.0

mmax = 100
m_w = np.arange(mmax + 1)

rc_w = np.real(coeff_R[:mmax+1])
rs_w = -np.imag(coeff_R[:mmax+1])
zc_w = np.real(coeff_Z[:mmax+1])
zs_w = -np.imag(coeff_Z[:mmax+1])

def fmt_mode(arr):
    return ", ".join(f"{int(v):10d}" for v in arr)

def fmt_coef(arr):
    return ", ".join(f"{v: .4E}".replace("E", "D") for v in arr)

with open("harmonics_w.txt", "w") as f:
    f.write(f"  m_w  = {fmt_mode(m_w)},\n")
    f.write(f"  rc_w = {fmt_coef(rc_w)},\n")
    f.write(f"  rs_w = {fmt_coef(rs_w)},\n")
    f.write(f"  zc_w = {fmt_coef(zc_w)},\n")
    f.write(f"  zs_w = {fmt_coef(zs_w)},\n")

# --- reconstruction from Fourier harmonics ---
R_rec = np.zeros_like(theta_uniform)
Z_rec = np.zeros_like(theta_uniform)

R_rec += rc_w[0]
Z_rec += zc_w[0]

for m in range(1, mmax + 1):
    R_rec += rc_w[m] * np.cos(m * theta_uniform) + rs_w[m] * np.sin(m * theta_uniform)
    Z_rec += zc_w[m] * np.cos(m * theta_uniform) + zs_w[m] * np.sin(m * theta_uniform)

# --- plots for comparison ---
import matplotlib.pyplot as plt

plt.figure(figsize=(7, 7))
plt.plot(R, Z, 'k-', lw=2, label='Interpolated contour')
plt.plot(R_rec, Z_rec, 'r--', lw=2, label=f'Fourier reconstruction (mmax={mmax})')
plt.plot(outline[:, 0], outline[:, 1], 'bo', ms=3, label='Input points')
plt.axis('equal')
plt.xlabel('R')
plt.ylabel('Z')
plt.legend()
plt.title('Contour vs Fourier reconstruction')
plt.tight_layout()
plt.show()

print(f"Max reconstruction error = {np.max(err):.6e}")
