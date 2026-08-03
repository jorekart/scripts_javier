import imas
import matplotlib.pyplot as plt

database = 'ITER'
user     = 'public'
shot     = 105033
run      = 1
dataversion = '3'
t_index  = 56

imas_entry = imas.DBEntry(imas.imasdef.MDSPLUS_BACKEND, database, shot, run, user, data_version = dataversion)
imas_entry.open()

summary = imas_entry.get('summary')

time = summary.time
Ip   = summary.global_quantities.ip.value

plt.xlabel(r'Time [s]')
plt.ylabel(r'$I_p$ [A]')
plt.grid()
plt.plot(time, Ip)
plt.show()


