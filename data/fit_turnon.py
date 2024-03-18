import math
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib as mpl
#import matplotlib.colors as colors
from matplotlib import gridspec
#from mpl_toolkits import mplot3d

from scipy.optimize import curve_fit


#plot parameters
DPI = 100
scale = 1.0
fs = 8*scale
fig = plt.figure(figsize=(8*scale,8.0*scale))
plt.rcParams.update({'font.size': fs})
gs = gridspec.GridSpec(nrows=1, ncols=1, width_ratios=[1], height_ratios=[1])




TS = np.loadtxt("out_TS")


def fitfunc(t,Gamma, omega, amp0, phi):
  return amp0*np.exp(-t*Gamma)*np.cos(t*omega+phi)



    
ydata = TS[:,1]
ydata = ydata - ydata[-1]
xdata = TS[:,0]

dt = xdata[1]-xdata[0]

Ymax = 0
YmaxPos = 0
for k in range(ydata.size):
  if ydata[k] > Ymax:
    Ymax = ydata[k]
    YmaxPos = k
    
YmaxPos_t = YmaxPos*dt
print Ymax
print YmaxPos_t

FitTimeInt = 10
FitTimeDelta = 2

xdataFit = xdata[YmaxPos+int(FitTimeDelta/dt):YmaxPos+int((FitTimeDelta+FitTimeInt)/dt)]-(YmaxPos_t+FitTimeDelta)
ydataFit = ydata[YmaxPos+int(FitTimeDelta/dt):YmaxPos+int((FitTimeDelta+FitTimeInt)/dt)]


popt, pcov = curve_fit(fitfunc, xdataFit, ydataFit)
print popt

print "RO: " + str(popt[1]/(2.0*np.pi))
print "Damping: " + str(popt[0])

print pcov



##################################################################################################################
##################################################################################################################

ax00 = plt.subplot(gs[0,0])
plt.grid(c='0.5')

#plt.plot(xdata,ydata)
plt.plot(xdataFit,ydataFit, c='r', lw=2)

plt.plot(xdataFit,fitfunc(xdataFit,popt[0],popt[1],popt[2],popt[3]))



plt.tight_layout()
#plt.subplots_adjust(hspace=0.07, wspace=0.07, left=0.10,right=0.96,top=0.98,bottom=0.07)
#plt.savefig("timetraces.pdf", dpi=600)
#plt.savefig("timetraces.png", dpi=600)
plt.show()

