import numpy as np

# This code creates a toroidal band discretized with triangles
# The band is created from a theta and phi range
# The triangles are stored in STL format for visualization in paraview

def to_cartesian(R,Z,phi):
    x = R*np.cos(phi)
    y = R*np.sin(phi)
    z = Z
    return np.array([x,y,z])*1000 #in mm

def to_RZ(tht, R0, Z0, r):
    R = R0 + r*np.cos(tht)
    Z = Z0 + r*np.sin(tht)
    return np.array([R,Z])

def compute_triangle_normal(p1,p2,p3):
    p21 = p2 - p1 
    p31 = p3 - p1 
    normal = np.cross(p21,p31)
    return normal/np.linalg.norm(normal)

phi_start = 0.0;   phi_end = 2*np.pi

# Data for the ITER DFW!!
#tht_start = -2*np.pi*15/360.0;   tht_end = 2*np.pi*15/360.0
#dR = 0.1 # From R.A. Pitts, the surface is recessed 10 cm
#R0        = 4.0;   Z0      = 0.63;   r_min = 8.4+0.1-R0

tht_start = -2*np.pi*15/360.0;   tht_end = 2*np.pi*15/360.0
dR        = 0.1 
R0        = 4.0;   Z0 = 0.63;   r_min = 8.4+0.1-R0

#nphi = 500;     ntht   = 100
nphi = 60;     ntht   = 40
dphi = (phi_end - phi_start) / float(nphi)
dtht = (tht_end - tht_start) / float(ntht)

f = open("test.stl", "w")

f.write("solid test \n")

for ip in range(nphi):
    for it in range(ntht):

        # 4 points for a quadrilateral
        phi1 = phi_start + float(ip  )*dphi
        phi2 = phi_start + float(ip+1)*dphi
        phi3 = phi2
        phi4 = phi1

        tht1 = tht_start + float(it  )*dtht
        tht2 = tht1
        tht3 = tht_start + float(it+1)*dtht
        tht4 = tht3

        RZ1 = to_RZ(tht1, R0, Z0, r_min)
        RZ2 = to_RZ(tht2, R0, Z0, r_min)
        RZ3 = to_RZ(tht3, R0, Z0, r_min)
        RZ4 = to_RZ(tht4, R0, Z0, r_min)

        r1 = to_cartesian(RZ1[0],RZ1[1],phi1)
        r2 = to_cartesian(RZ2[0],RZ2[1],phi2)
        r3 = to_cartesian(RZ3[0],RZ3[1],phi3)
        r4 = to_cartesian(RZ4[0],RZ4[1],phi4)

        # output first triangle
        normal = compute_triangle_normal(r1,r2,r3)
        f.write("  facet normal %s %s %s \n"%(normal[0],normal[1],normal[2]))
        f.write("    outer loop \n")
        f.write("      vertex %s %s %s \n"%(r1[0],r1[1],r1[2]))
        f.write("      vertex %s %s %s \n"%(r2[0],r2[1],r2[2]))
        f.write("      vertex %s %s %s \n"%(r3[0],r3[1],r3[2]))
        f.write("    endloop \n")
        f.write("  endfacet \n")

        # output second triangle
        normal = compute_triangle_normal(r1,r3,r4)
        f.write("  facet normal %s %s %s \n"%(normal[0],normal[1],normal[2]))
        f.write("    outer loop \n")
        f.write("      vertex %s %s %s \n"%(r1[0],r1[1],r1[2]))
        f.write("      vertex %s %s %s \n"%(r3[0],r3[1],r3[2]))
        f.write("      vertex %s %s %s \n"%(r4[0],r4[1],r4[2]))
        f.write("    endloop \n")
        f.write("  endfacet \n")
        
f.write("endsolid \n")
f.close()
