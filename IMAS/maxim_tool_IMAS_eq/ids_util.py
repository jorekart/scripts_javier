
import sys
import math
import imas # UAL library
import matplotlib
matplotlib.use('Qt5Agg')
#matplotlib.use('GTK3Agg')

import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches

import numpy as np



def GetGeometryPath(geom):
  verts = []
  codes = []
  
  if (geom.geometry_type == 2):
    # Rectangle
    verts = [
      (geom.rectangle.r - 0.5*geom.rectangle.width, geom.rectangle.z - 0.5*geom.rectangle.height),  # left, bottom
      (geom.rectangle.r - 0.5*geom.rectangle.width, geom.rectangle.z + 0.5*geom.rectangle.height),  # left, top
      (geom.rectangle.r + 0.5*geom.rectangle.width, geom.rectangle.z + 0.5*geom.rectangle.height),  # right, top
      (geom.rectangle.r + 0.5*geom.rectangle.width, geom.rectangle.z - 0.5*geom.rectangle.height),  # right, bottom
      (0., 0.),  # ignored
    ]
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]
    path = Path(verts, codes)
    return path
  
  elif (geom.geometry_type == 3):
    # Oblique
    verts = [
      (geom.oblique.r, geom.oblique.z),  # left, bottom
      (geom.oblique.r + geom.oblique.length_alpha*math.cos(geom.oblique.alpha),
        geom.oblique.z + geom.oblique.length_alpha*math.sin(geom.oblique.alpha)),  # left, top
      (geom.oblique.r + geom.oblique.length_alpha*math.cos(geom.oblique.alpha) -  geom.oblique.length_beta*math.sin(geom.oblique.beta),
        geom.oblique.z + geom.oblique.length_alpha*math.sin(geom.oblique.alpha) +  geom.oblique.length_beta*math.cos(geom.oblique.beta)),  # right, bottom
      (geom.oblique.r - geom.oblique.length_beta*math.sin(geom.oblique.beta),
        geom.oblique.z + geom.oblique.length_beta*math.cos(geom.oblique.beta)),  # right, top
      (0., 0.),  # ignored
    ]
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]
    path = Path(verts, codes)
    return path
  
  elif (geom.geometry_type == 1):
    # Outline
    n = len(geom.outline.r)
    for i in range(n):
      verts.append((geom.outline.r[i], geom.outline.z[i]))
      codes.append(Path.LINETO)
      
    if (n > 0):
      codes[0] = Path.MOVETO
      
      verts.append((0.0, 0.0))
      codes.append(Path.CLOSEPOLY)
    
    path = Path(verts, codes)
    return path
    
  elif (geom.geometry_type == 5):
    # Annulus
    #print(dir(patches))
    #width = geom.annulus.radius_outer - geom.annulus.radius_inner
    circle_out = Path.circle([geom.annulus.r, geom.annulus.z], geom.annulus.radius_outer);
    circle_in = Path.circle([geom.annulus.r, geom.annulus.z], geom.annulus.radius_inner);
    vertices = circle_in.vertices[::-1]
    codes = circle_in.codes
    circle_in = Path(vertices, codes)
    circle = Path.make_compound_path(circle_out, circle_in)
    # Get the path and the affine transformation
    #path = circle.get_path()
    #transform = circle.get_transform()
    
    # Now apply the transform to the path
    #newpath = transform.transform_path(path)
    return circle
  
  elif (geom.geometry_type == 6):
    # Thick line
    verts.append((geom.thick_line.first_point.r, geom.thick_line.first_point.z))
    codes.append(Path.LINETO)
      
    verts.append((geom.thick_line.second_point.r, geom.thick_line.second_point.z))
    codes.append(Path.LINETO)

    codes[0] = Path.MOVETO
      
    verts.append((0.0, 0.0))
    codes.append(Path.CLOSEPOLY)
    
    path = Path(verts, codes)
    return path
  
  path = Path(verts, codes)
  return path


def plot_pf_active(ax, ids, facecolor='orange', edgecolor='blue'):
    for coil in ids.coil:
      for elem in coil.element:
        path = GetGeometryPath(elem.geometry)
        patch = patches.PathPatch(path, facecolor=facecolor, edgecolor=facecolor, lw=1)
        ax.add_patch(patch)
        
        
def plot_pf_passive(ax, ids, facecolor=(0, 0, 0.5), edgecolor=(0, 0, 1)):
    for loop in ids.loop:
      for elem in loop.element:
        path = GetGeometryPath(elem.geometry)
        if elem.turns_with_sign < 0.:
          facecolor = (0, 0.5, 0.5)
        else:
          facecolor = (0.5, 0.5, 0)
        patch = patches.PathPatch(path, facecolor=facecolor, edgecolor=edgecolor)
        ax.add_patch(patch)
        
        
def plot_limiter(ax, ids, color='k-'):
    if (len(ids.description_2d) > 0):
      for unit in ids.description_2d[0].limiter.unit:
        ax.plot(unit.outline.r, unit.outline.z, color, linewidth=1, label='limiter')
    else:
      print("No limiter data in given IDS")
          
          

def GetIDS(meta, idslist):
    imas_entry_init = imas.DBEntry(imas.imasdef.MDSPLUS_BACKEND, meta["database"], meta["shot"], meta["run"], meta["user"], data_version = '3')
    imas_entry_init.open()
    
    if (isinstance(idslist, dict)):
      out = {}
      for key in idslist:
        out[key] = imas_entry_init.get(key)
      
    if (isinstance(idslist, list)):
      out = {}
      for it in idslist:
        out[it] = imas_entry_init.get(it)
      
    elif (isinstance(idslist, str)):
      out = imas_entry_init.get(idslist)
        
    imas_entry_init.close()
    return out

    