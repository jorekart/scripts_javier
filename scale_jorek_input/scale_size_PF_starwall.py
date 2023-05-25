import numpy as np


fR   = 0.5

file1 = open('JET_PF_coils.nml', 'r')
Lines = file1.readlines()

with open('JET_PF_coils.nml', "r") as f:
    lines = f.readlines()

with open("test.txt", "w") as f:
    for line in lines:
        indexR = line.find("R_fila          =")
        indexZ = line.find("Z_fila          =")
        if indexR != -1:
            s=line[indexR+17:]
            arr = np.array(s.strip("[]").split(","), dtype=float)
            arr = arr * fR
            new_array = "{}".format(", ".join(str(x) for x in arr))
            
            new_line = line[:indexR+17] + new_array + "\n"
            f.write(new_line)
        elif indexZ != -1:
            s=line[indexZ+17:]
            arr = np.array(s.strip("[]").split(","), dtype=float)
            arr = arr * fR
            new_array = "{}".format(", ".join(str(x) for x in arr))
            
            new_line = line[:indexZ+17] + new_array + "\n"
            f.write(new_line)

        else:
            f.write(line)

