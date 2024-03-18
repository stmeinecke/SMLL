import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib import gridspec


#plot parameters
DPI = 100
scale = 1.0
scaleCM = 1/2.54
fs = 10*scale
fig = plt.figure(figsize=(5*scale,3.5*scale), dpi=DPI)
plt.rcParams.update({'font.size': fs})



filename = 'out_TS'
data =  np.loadtxt(filename)


CCRTT = 4.0
dt = data[1,0]-data[0,0]
print dt

OSTime = 0
OSSteps = int(OSTime / dt)
TS = data[OSSteps:,1]

SPRT = int(CCRTT/dt)
print SPRT 

Nmax = TS.size
print Nmax

RTs = np.floor(Nmax/SPRT)
print RTs

space_time = np.reshape(TS[0:int(SPRT*RTs)], (int(RTs),int(SPRT))) 

space_time = space_time/space_time.max()



extentarray = [0, CCRTT, 0,RTs]

#plt.imshow(np.flipud(space_time), aspect="auto", cmap="jet", interpolation="bicubic", extent=extentarray, norm=mpl.colors.LogNorm(), vmin = 1E-3)
plt.imshow(np.flipud(space_time), aspect="auto", cmap="jet", interpolation="none", extent=extentarray)
cb0 = plt.colorbar(cmap=space_time, spacing='uniform')
cb0.set_label("Normalized Intensity")
plt.xlabel(r'Time $t$ in $\tau$')
plt.ylabel(r'Roundtrip Number')


plt.tight_layout()
#plt.savefig(filename+"space_time.pdf")
plt.savefig(filename+"_space_time.png", dpi=300)
plt.show()
