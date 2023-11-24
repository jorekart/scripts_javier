# G-EQDSK format description
# https://w3.pppl.gov/ntcc/TORAY/G_EQDSK.pdf

import numpy as np

# Amount of columns in the output file
ncol = 5

def writeReal(f, v):
  f.write("% 16.9E"%v)

def write1d(f, v):
  n = len(v)
  for i in range(n):
    writeReal(f, v[i])
    if i>0 and ((i+1)%ncol == 0 or i == n-1):
      f.write("\n")

def write2d(f, v):
  sh = np.shape(v)
  n  = sh[0]
  m  = sh[1]
  k  = 0
  for j in range(m):
    for i in range(n):
      writeReal(f, v[i,j])
#      k = i*m+j
      k += 1
#      if k>0 and ((k+1)%ncol == 0 or k == (n-1)*m+(m-1)):
      if k>ncol-1 or (j==(m-1) and (i==n-1)):
        k = 0
        f.write("\n")


def write(f, data):
  #Write a G-EQDSK file

  if isinstance(f, str):
      # If the input is a string, treat as file name
      with open(f, "w") as fh: # Ensure file is closed
          return write(fh, data) # Call again with file object


  # Write description
  f.write("geqdsk 0 " + str(data['nr']) + " " + str(data['nz']) + "\n")

  writeReal(f, data['rdim'])
  writeReal(f, data['zdim'])
  writeReal(f, data['rcentr'])
  writeReal(f, data['rleft'])
  writeReal(f, data['zmid'])
  f.write("\n")

  writeReal(f, data['rmaxis'])
  writeReal(f, data['zmaxis'])
  writeReal(f, data['simagx'])
  writeReal(f, data['sibdry'])
  writeReal(f, data['bcentr'])
  f.write("\n")

  writeReal(f, data['current'])
  writeReal(f, data['simagx'])
  writeReal(f, 0.0)
  writeReal(f, data['rmaxis'])
  writeReal(f, 0.0)
  f.write("\n")

  writeReal(f, data['zmaxis'])
  writeReal(f, 0.0)
  writeReal(f, data['sibdry'])
  writeReal(f, 0.0)
  writeReal(f, 0.0)
  f.write("\n")

  # Write the arrays
  write1d(f, data['fpol'])
  write1d(f, data['pressure'])
  write1d(f, data['ffprime'])
  write1d(f, data['pprime'])
  write2d(f, data['psirz'])
  write1d(f, data['qpsi'])

  # Boundary and limiter
  f.write(" " + str(data['nbbbs']) + " " + str(data['limitr']) + "\n")
  write1d(f, data['rbbbs'])
  write1d(f, data['zbbbs'])
  write1d(f, data['rlim'])
  write1d(f, data['zlim'])
    
