from matplotlib.transforms import Bbox, TransformedBbox, \
     blended_transform_factory

from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector,\
     BboxConnectorPatch

from numpy import sqrt

def connect_bbox(bbox1, bbox2,
                 loc1a, loc2a, loc1b, loc2b,
                 prop_lines, prop_patches=None):
    if prop_patches is None:
        prop_patches = prop_lines.copy()
        prop_patches["alpha"] = prop_patches.get("alpha", 1)*0.2

    c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
    c1.set_clip_on(False)
    c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
    c2.set_clip_on(False)

    bbox_patch1 = BboxPatch(bbox1, **prop_patches)
    bbox_patch2 = BboxPatch(bbox2, **prop_patches)

    p = BboxConnectorPatch(bbox1, bbox2,
                           #loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                           loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b,
                           **prop_patches)
    p.set_clip_on(False)

    return c1, c2, bbox_patch1, bbox_patch2, p

def zoom_effect02(ax1, ax2, **kwargs):
    u"""
    ax1 : the main axes
    ax1 : the zoomed axes

    Similar to zoom_effect01.  The xmin & xmax will be taken from the
    ax1.viewLim.
    """

    tt = ax1.transScale + (ax1.transLimits + ax2.transAxes)
    trans = blended_transform_factory(ax2.transData, tt)

    mybbox1 = ax1.bbox
    mybbox2 = TransformedBbox(ax1.viewLim, trans)

    prop_patches=kwargs.copy()
    prop_patches["ec"]="none"
    prop_patches["alpha"]=0.2

    c1, c2, bbox_patch1, bbox_patch2, p = \
        connect_bbox(mybbox1, mybbox2,
                     loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                     prop_lines=kwargs, prop_patches=prop_patches)

    ax1.add_patch(bbox_patch1)
    ax2.add_patch(bbox_patch2)
    ax2.add_patch(c1)
    ax2.add_patch(c2)
    ax2.add_patch(p)

    return c1, c2, bbox_patch1, bbox_patch2, p

def volcalc(j,dj,E,dE,u2,du2,e):
  Units=e**2*1.38e-23/8.85e-12
  V=(E*j)**2*Units/u2
  dj1=E*dj
  dE1=j*dE
  du1=E*j*sqrt(du2/u2)
  S=2*E*j/u2
  dV=S*sqrt(dj1**2+dE1**2+du1**2)*Units
  return V,dV,Units

def dimcalc(j,dj,E,dE,u2,du2,e):
  V,dV,Units=volcalc(j,dj,E,dE,u2,du2,e)
  print V,'+',dV,'m3'
  L=V**(1/3.0)
  dL=V**(-2/3.0)*dV/3.0
  print L,'+',dL,'m'
  print L/16,'+',dL/16,'m'

if __name__ == "__main__":
  import matplotlib.pyplot as plt
  from numpy import arange,sin

  ax1 = plt.subplot(211)
  plt.plot(x,y,x,y**2)
  ax2 = plt.subplot(212)
  plt.plot(x,y,'og-')
  ax1.set_xlim(0.5, 0.7)
  zoom_effect02(ax1, ax2)


  plt.show()
