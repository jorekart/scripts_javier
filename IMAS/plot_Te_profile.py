import imas
import matplotlib.pyplot as plt

database = 'ITER'
user     = 'public'
shot     = 105033
run      = 1
t_index  = 56

imas_entry = imas.DBEntry(imas.imasdef.MDSPLUS_BACKEND, database, shot, run, user, data_version = '3')
imas_entry.open()

core_profiles = imas_entry.get('core_profiles')

rho = core_profiles.profiles_1d[t_index].grid.rho_tor_norm
Te  = core_profiles.profiles_1d[t_index].electrons.temperature

plt.xlabel(r'$\rho_{norm}$')
plt.ylabel(r'$T_e$')
plt.grid()
plt.plot(rho, Te)
plt.show()


