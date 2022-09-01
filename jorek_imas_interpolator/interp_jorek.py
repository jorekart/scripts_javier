import numpy as np
import multiprocessing
import time

from functools import partial
################################
# Classes for elements and nodes
################################
class element_class:
    def __init__(self):
        self.vertex_ind = np.empty(0)
        self.size       = np.empty(0)

class element_list_class:
    def __init__(self):
        self.n_elm    = 0 
        self.n_vertex = 0 
        self.n_dof    = 0 
        self.elems    = []

class node_class:
    def __init__(self):
        self.x   = np.empty(0)

class node_list_class:
    def __init__(self):
        self.n_nodes  = 0
        self.n_tor    = 0 
        self.n_period = 0 
        self.nodes    = []


################################################
# Reconstruct node and element list as in JOREK
################################################
def reconstruct_element_list(grid):
    # grid comes from ids%grid_ggd(grid_ind)
    elm_list_imas  = np.array( grid.space[0].objects_per_dimension[2].object.array )

    elm_list       = element_list_class()
    elm_list.n_elm = len(elm_list_imas) 

    for elm in elm_list_imas:
        elm_tmp = element_class()
        elm_tmp.size       = elm.geometry_2d.transpose() 
        elm_tmp.vertex_ind = elm.nodes 
        elm_list.elems.append( elm_tmp )

    elm_list.n_vertex = len(elm_tmp.size[:,0]) 
    elm_list.n_dof    = len(elm_tmp.size[0,:]) 


    return elm_list


def reconstruct_node_list(grid):
    # grid comes from ids%grid_ggd(grid_ind)
    node_list_imas    = np.array( grid.space[0].objects_per_dimension[0].object.array )

    node_list         = node_list_class()
    node_list.n_nodes = len(node_list_imas) 

    for node in node_list_imas:
        node_tmp      = node_class()
        node_tmp.x    = node.geometry_2d.transpose() 
        node_list.nodes.append(node_tmp)

    node_list.n_period = grid.space[1].geometry_type.index

    return node_list
 
################################################
# Interpolation routines for space and values  
################################################

#!> This subroutine interpolates a specific position within one element at a given position (s,t)
def interp_RZ(node_list, element_list, i_elm, s, t):

    basis    = basis_functions(s, t)
    sizes    = element_list.elems[i_elm].size

    coeff_R  = np.empty_like(sizes)
    coeff_Z  = np.empty_like(sizes)

    for kv in range(0, element_list.n_vertex):  
        iv             = element_list.elems[i_elm].vertex_ind[kv] - 1
        coeff_R[kv,:]  = node_list.nodes[iv].x[:,0] 
        coeff_Z[kv,:]  = node_list.nodes[iv].x[:,1] 

    # i=vertex, j=dof
    R = np.einsum('ij,ij,ij->', coeff_R, sizes, basis ) 
    Z = np.einsum('ij,ij,ij->', coeff_Z, sizes, basis ) 

    return R, Z


#!> Same as interp_RZ but including derivatices
def interp_RZ_deriv1(node_list, element_list, i_elm, s, t):

    basis    = basis_functions(s, t)
    basis_s  = basis_functions_s(s, t)
    basis_t  = basis_functions_t(s, t)

    sizes    = element_list.elems[i_elm].size

    coeff_R  = np.empty_like(sizes)
    coeff_Z  = np.empty_like(sizes)

    for kv in range(0, element_list.n_vertex):  
        iv             = element_list.elems[i_elm].vertex_ind[kv] - 1
        coeff_R[kv,:]  = node_list.nodes[iv].x[:,0] 
        coeff_Z[kv,:]  = node_list.nodes[iv].x[:,1] 

    # i=vertex, j=dof
    R  = np.einsum('ij,ij,ij->', coeff_R, sizes, basis   ) 
    Z  = np.einsum('ij,ij,ij->', coeff_Z, sizes, basis   ) 

    Rs = np.einsum('ij,ij,ij->', coeff_R, sizes, basis_s ) 
    Zs = np.einsum('ij,ij,ij->', coeff_Z, sizes, basis_s ) 

    Rt = np.einsum('ij,ij,ij->', coeff_R, sizes, basis_t ) 
    Zt = np.einsum('ij,ij,ij->', coeff_Z, sizes, basis_t ) 

    return R, Rs, Rt, Z, Zs, Zt


#!> This subroutine interpolates space a specific position within one element at a given position (s,t)
def interp_RZ_fast(s, t, coeff_R, coeff_Z):

    basis    = basis_functions(s, t)
    basis_s  = basis_functions_s(s, t)
    basis_t  = basis_functions_t(s, t)

    R  = np.sum( coeff_R * basis  )
    Rs = np.sum( coeff_R * basis_s)
    Rt = np.sum( coeff_R * basis_t)
    Z  = np.sum( coeff_Z * basis  )
    Zs = np.sum( coeff_Z * basis_s)
    Zt = np.sum( coeff_Z * basis_t)

    return R, Rs, Rt, Z, Zs, Zt


#!> This subroutine interpolates a value within one element at a given position (s,t, phi)
def interp_val(val_coeff, node_list, element_list, i_elm, s, t, phi):

    n_vertex  = element_list.n_vertex
    n_dof     = element_list.n_dof   
    n_tor     = len( val_coeff[:,0] )

    basis     = basis_functions(s, t)
    tor_basis = toroidal_basis(n_tor, node_list.n_period, phi, False)
    sizes     = element_list.elems[i_elm].size

    coeff     = np.zeros( shape=( n_vertex, n_dof, n_tor  ) )

    for kv in range(0, n_vertex):  

        inode  = element_list.elems[i_elm].vertex_ind[kv] - 1 

        for idof in range(0, n_dof):

            icoeff              =  inode + idof*node_list.n_nodes 
            coeff[kv,idof,:]  =  val_coeff[:,icoeff]

    # i=vertex, j=dof,  k=harmonic
    return np.einsum('ijk,ij,ij,k->', coeff, sizes, basis, tor_basis ) 



def value_at_point(val_coeff, node_list, elm_list, phi, RZ_elm, RZv):

    # Find RZ takes most of the time...
    jor_coords =  find_RZ(node_list, elm_list, RZv[0], RZv[1], RZ_elm) 

    i_elm = jor_coords[0]
    s_out = jor_coords[1]
    t_out = jor_coords[2]
    ifail = jor_coords[3]

    if (ifail == 0):
        val = interp_val(val_coeff, node_list, elm_list, i_elm, s_out, t_out, phi)
    else:
        print('Failed finding point at R, Z' + str(RZv[0])+', '+ str(RZv[1]))
        val = 0.0

    return val



# Gets value at R, Z, phi coordinates given in 1D arrays
def value_at_RZphi(grid, val_coeff, Rp, Zp, phi):
 
    print( '  Reconstruct node/element list') 
    elm_list  = reconstruct_element_list( grid )
    node_list = reconstruct_node_list(    grid )

    # Sanity check
    if ( (len(Rp) != len(Zp)) or (len(Rp) != len(phi)) or (len(Zp) != len(phi)) ):
        print( " The length of the coordinate arrays are different! ")
        return

    # List with element center coordinates, for faster find_RZ
    n_elm    = elm_list.n_elm
    RZ_elm   = np.zeros( (n_elm, 2) )

    for i in range(0,n_elm):
        outRZ       = interp_RZ(node_list, elm_list, i, s=0.5, t=0.5)
        RZ_elm[i,:] = np.array( [ outRZ[0], outRZ[1]] )

    parallel = True   # Activate or nor parallel version

    print( '  Calculating values at given points...') 
    t0  = time.time()

    if (parallel):
        # Auxiliary routine that only takes RZ as input
        aux2 = partial( value_at_point, val_coeff, node_list, elm_list, 0, RZ_elm)

        # Multiprocessing
        #pool = multiprocessing.Pool()
        pool = multiprocessing.Pool(8)  # 8 is faster than default in ITER SDCC
        RZv = zip( Rp, Zp )
        val = pool.map(aux2, RZv)
        pool.close()

    else:
        val = np.zeros_like( Rp )
        for i in range(0, len(Rp)):
            val[i] = value_at_point(val_coeff, node_list, elm_list, phi[i], RZ_elm,np.array( [Rp[i], Zp[i]]))
    
    t1 = time.time()
    print('Elapsed time in expensive point loop = ' + str(t1-t0) )

    return val


####################################################
# Routines to find RZ points in JOREK representation
####################################################

def find_RZ_single(node_list,element_list,i_elm,R_find,Z_find):
#-------------------------------------------------------------------------
#< solves two non-linear equations using Newtons method (from numerical recipes)
#< LU decomposition replaced by explicit solution of 2x2 matrix.
#<
#< finds the crossing of two coordinate lines given as a series of cubics in element
#< i_elm
#-------------------------------------------------------------------------

    ntrial = 20
    tolx = 1.e-8
    tolf = 1.e-15
    
    ##### Save element properties (for speed) #####
    sizes    = element_list.elems[i_elm].size

    coeff_R  = np.empty_like(sizes)
    coeff_Z  = np.empty_like(sizes)

    for kv in range(0, element_list.n_vertex):  
        iv             = element_list.elems[i_elm].vertex_ind[kv] - 1
        coeff_R[kv,:]  = node_list.nodes[iv].x[:,0] * sizes[kv,:] 
        coeff_Z[kv,:]  = node_list.nodes[iv].x[:,1] * sizes[kv,:]
    ################################################################

    x    = np.zeros(2)
    p    = np.zeros(2)
    fvec = np.zeros(2)
    fjac = np.zeros((2,2))

    #  Initialize parameters for output
    ifail = 999;  R_out = 0.0;  Z_out = 0.0; s_out = 0.0; t_out = 0.0;
    
    for istart in range(1,6):
    
        if (istart == 1): 
            x[0] = 0.5
            x[1] = 0.5
        elif (istart == 2):
            x[0] = 0.75
            x[1] = 0.75
        elif (istart == 3):
            x[0] = 0.75
            x[1] = 0.25
        elif (istart == 4):
            x[0] = 0.25
            x[1] = 0.75
        elif (istart == 5):
            x[0] = 0.25
            x[1] = 0.25
    
        for i in range(1,ntrial+1):
    
            out_interp = interp_RZ_fast(x[0],x[1], coeff_R, coeff_Z)

            RRg1     = out_interp[0]
            ZZg1     = out_interp[3]
        
            fvec[0]   = RRg1 - R_find
            fvec[1]   = ZZg1 - Z_find
            fjac[0,0] = out_interp[1]
            fjac[0,1] = out_interp[2]
            fjac[1,0] = out_interp[4]
            fjac[1,1] = out_interp[5]
        
            errf=np.absolute(fvec[0])+np.absolute(fvec[1])
        
#            print(' %s %s %s %s %s %s %s '%(i,errf,x,RRg1,R_find,ZZg1,Z_find))

            if (errf <= tolf):
        
                s_out     = x[0]
                t_out     = x[1]
        
                ielm      = i_elm
                R_out     = RRg1
                Z_out     = ZZg1
        
                ifail = 0
                return R_out,Z_out,i_elm,s_out,t_out,ifail
        
            p    = -fvec
            dis  = fjac[1,1]*fjac[0,0]-fjac[0,1]*fjac[1,0]
            tmp  = p[0]

            if (dis != 0.0):
                p[0] = (fjac[1,1]*p[0]-fjac[0,1]*p[1])/dis
                p[1] = (fjac[0,0]*p[1]-fjac[1,0]*tmp )/dis
            else:
                break
        
            errx=np.absolute(p[0]) + np.absolute(p[1])

            for k in range(0,2):        
                p[k] = np.minimum(p[k], 0.25)
                p[k] = np.maximum(p[k],-0.25)
        
            x += p
        
            for k in range(0,2):        
                x[k] = np.maximum(x[k], 0.0)
                x[k] = np.minimum(x[k], 1.0)
        
            if (errx <= tolx):
        
                s_out     = x[0]
                t_out     = x[1]
        
                R_out     = RRg1
                Z_out     = ZZg1
        
                ifail = 0
                return R_out,Z_out,i_elm,s_out,t_out,ifail
    

        return R_out,Z_out,i_elm,s_out,t_out,ifail



def find_RZ(node_list, element_list, R_find, Z_find, RZ_elm):
    
    # Look at elements which center is closer to specified point
    # Vector with distances to point
    Rvec = np.full_like( RZ_elm[:,0], R_find )
    Zvec = np.full_like( RZ_elm[:,1], Z_find )

    distances = (RZ_elm[:,0]-Rvec[:])**2 + (RZ_elm[:,1]-Zvec[:])**2    

    # Try to find point in the element with minimum distance
    ntry_elm = 10;  

    for itry in range(0,ntry_elm):    
        i_elm = np.argmin(distances)

        out   = find_RZ_single(node_list,element_list,i_elm,R_find,Z_find)

        if (out[5] ==0):  # out[5] is ifail
            break 
        else:
            distances[i_elm] = 1e99

    return i_elm, out[3], out[4], out[5]


"""
Calculate values of the basis functions at positions s and t
Optionally put many values of s and t at once as numpy arrays.
Dimension 0: vertex
Dimension 1: order 
"""
def basis_functions(s,t):
    return np.asarray([
      [ (-1 + s)**2*(1 + 2*s)*(-1 + t)**2*(1 + 2*t),
        3*(-1 + s)**2*s*(-1 + t)**2*(1 + 2*t),
        3*(-1 + s)**2*(1 + 2*s)*(-1 + t)**2*t,
        9*(-1 + s)**2*s*(-1 + t)**2*t],
      [ -(s**2*(-3 + 2*s)*(-1 + t)**2*(1 + 2*t)),
        -3*(-1 + s)*s**2*(-1 + t)**2*(1 + 2*t),
        -3*s**2*(-3 + 2*s)*(-1 + t)**2*t,
        -9*(-1 + s)*s**2*(-1 + t)**2*t],
      [ s**2*(-3 + 2*s)*t**2*(-3 + 2*t),
        3*(-1 + s)*s**2*t**2*(-3 + 2*t),
        3*s**2*(-3 + 2*s)*(-1 + t)*t**2,
        9*(-1 + s)*s**2*(-1 + t)*t**2],
      [ -((-1 + s)**2*(1 + 2*s)*t**2*(-3 + 2*t)),
        -3*(-1 + s)**2*s*t**2*(-3 + 2*t),
        -3*(-1 + s)**2*(1 + 2*s)*(-1 + t)*t**2,
        -9*(-1 + s)**2*s*(-1 + t)*t**2]])

"""
Calculate values of the basis functions derived to s at positions s and t
Optionally put many values of s and t at once as numpy arrays.
Dimension 0: vertex
Dimension 1: order 
optional dimension 2, 3: position s, t
"""
 
def basis_functions_s(s,t):
    return np.asarray([
    [  6*(-1 + s)*s*(-1 + t)**2*(1 + 2*t),
       3*(-1 + s)*(-1 + 3*s)*(-1 + t)**2*(1 + 2*t),
      18*(-1 + s)*s*(-1 + t)**2*t,
       9*(-1 + s)*(-1 + 3*s)*(-1 + t)**2*t],
    [ -6*(-1 + s)*s*(-1 + t)**2*(1 + 2*t),
      -3*s*(-2 + 3*s)*(-1 + t)**2*(1 + 2*t),
      -18*(-1 + s)*s*(-1 + t)**2*t,
      -9*s*(-2 + 3*s)*(-1 + t)**2*t],
    [  6*(-1 + s)*s*t**2*(-3 + 2*t),
       3*s*(-2 + 3*s)*t**2*(-3 + 2*t),
      18*(-1 + s)*s*(-1 + t)*t**2,
      9*s*(-2 + 3*s)*(-1 + t)*t**2],
    [ -6*(-1 + s)*s*t**2*(-3 + 2*t),
       3*(1 - 3*s)*(-1 + s)*t**2*(-3 + 2*t),
     -18*(-1 + s)*s*(-1 + t)*t**2,
      9*(1 - 3*s)*(-1 + s)*(-1 + t)*t**2]])
 

"""
Calculate values of the basis functions derived to t at positions s and t
Optionally put many values of s and t at once as numpy arrays.
Dimension 0: vertex
Dimension 1: order 
optional dimension 2, 3: position s, t
"""
def basis_functions_t(s,t):

    return np.asarray([
    [  6*(-1 + s)**2*(1 + 2*s)*(-1 + t)*t,
       18*(-1 + s)**2*s*(-1 + t)*t,
       3*(-1 + s)**2*(1 + 2*s)*(-1 + t)*(-1 + 3*t),
       9*(-1 + s)**2*s*(-1 + t)*(-1 + 3*t)],
    [  -6*s**2*(-3 + 2*s)*(-1 + t)*t,
       -18*(-1 + s)*s**2*(-1 + t)*t,
       3*s**2*(-3 + 2*s)*(1 - 3*t)*(-1 + t),
       9*(-1 + s)*s**2*(1 - 3*t)*(-1 + t)],
    [  6*s**2*(-3 + 2*s)*(-1 + t)*t,
       18*(-1 + s)*s**2*(-1 + t)*t,
       3*s**2*(-3 + 2*s)*t*(-2 + 3*t),
       9*(-1 + s)*s**2*t*(-2 + 3*t)],
    [  -6*(-1 + s)**2*(1 + 2*s)*(-1 + t)*t,
       -18*(-1 + s)**2*s*(-1 + t)*t,
       -3*(-1 + s)**2*(1 + 2*s)*t*(-2 + 3*t),
       -9*(-1 + s)**2*s*t*(-2 + 3*t)]])


def toroidal_basis(n_tor, n_period, phi, without_n0_mode):
    # Setup toroidal coefficients for each plane and toroidal harmonic
    HZ = np.zeros(n_tor)
    for i in range(n_tor):
        mode = np.floor((i+1)/2)*n_period
        if (i == 0):
            if (not without_n0_mode):
                HZ[i] = 1
        elif (i % 2 == 0):
            HZ[i] = np.sin(mode*phi)
        elif (i % 2 == 1):
            HZ[i] = np.cos(mode*phi)
    return HZ



