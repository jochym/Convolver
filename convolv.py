#!/usr/bin/python
# -*- coding: utf-8 -*-
# $Id: convolv.py 82 2007-12-14 09:02:08Z Paweł T. Jochym $

ExpXScale=0.2417963
ExpYScale=0.1
#ExpXScale=1
#ExpYScale=1

from pylab import *
import sys,string
from Scientific.Functions.Interpolation import InterpolatingFunction
from Scientific import N
from matplotlib.widgets import Slider, Button, CheckButtons
from scipy.integrate import simps
import scipy

def Gauss(x,hw):
    s=hw/(2*sqrt(log(2)))
    return exp(-(x/s)**2)/(s*sqrt(pi))
    
def Lorenz(x,hw):
    a=hw/2
    return 1/(a*pi*(1+(x/a)**2))
        

def updateG(val):
    global hw, da, cx, cg, xmed, gauss
    hw = fwhmG.val
    gauss=Gauss(cx-xmed,hw)*step
    cg=scipy.convolve(da, gauss, mode=1)
    gpl.set_ydata(cg)
    #print integral(cx,da), integral(cx,cg)
    draw()

def updateL(val):
    global hw, da, cx, cl
    hw = fwhmL.val
    lorenz=Lorenz(cx-xmed,hw)*step
    cl=scipy.convolve(da, lorenz, mode=1)
    lpl.set_ydata(cl)    
    #print integral(cx,da), integral(cx,cl)
    draw()


def integral(dx,da):
    return simps(da,x=dx)

def changeData(label):
    if label=='Gauss' :
        gpl.set_visible(not gpl.get_visible())
    elif label=='Lorenz': 
        lpl.set_visible(not lpl.get_visible())
    elif label=='Theory': 
        dpl.set_visible(not dpl.get_visible())
    elif label=='Experiment': 
        epl.set_visible(not epl.get_visible())
    else :
        pass
    draw()



try:
    data=open(sys.argv[1],'r').readlines()
    data=[[string.atof(x) for x in l.split()] for l in data if len(l.split())>1]
    data=[[d[0] for d in data],[d[1] for d in data]]
except IndexError:
    print "Give data file to convolve as first argument"
    raise



xmax=max(data[0])
xmin=min(data[0])
xrng=xmax-xmin
xmed=(xmin+xmax)/2.

hw=xrng/10

stepd=(xmax-xmin)/max(len(data[0]),100)
steph=hw/10
step=min(stepd,steph)

cx=N.arrayrange(xmin,xmax+step,step)
xmin=cx[0]
xmax=cx[-1]
xmed=(xmin+xmax)/2.

di=InterpolatingFunction((N.array(data[0]),), N.array(data[1]), 0)
da=N.array([di(x) for x in cx])

try:
    edata=open(sys.argv[2],'r').readlines()
    edata=[[string.atof(x) for x in l.split()] for l in edata if len(l.split())>1]
    edata=[[d[0]*ExpXScale for d in edata],[d[1]*ExpYScale for d in edata]]
    edata=array(edata)
    edata[1]=edata[1]*integral(cx,da)/integral(edata[0],edata[1])
    #print edata
except IndexError:
    edata=[[],[]]


#print xmin, xmax, xmed, step, cx, cx-xmed
gauss=step*Gauss(cx-xmed,hw)
lorenz=step*Lorenz(cx-xmed,hw)
cg=scipy.convolve(da, gauss, mode=1)
cl=scipy.convolve(da, lorenz, mode=1)

#tstx=array([-hw/2,hw/2])
#tsty=5*step*Gauss(tstx,hw)
#print tstx, tsty
#plot(tstx+xmed,tsty,'y-')
#print da
#print integral(cx,gauss)/step, integral(cx,lorenz)/step, integral(cx,da), integral(cx,cg), integral(cx,cl)
subplot(111)
dpl,gpl,lpl,epl=plot(data[0], data[1], 'k', cx, cg, 'r', cx, cl, 'b',edata[0] , edata[1], 'go')
grid(True)
title('Linear convolution')
ylabel('Data')


#subplot(212)
#rx=arange(xmin,xmax,0.05*s)
#print xmin, xmax, x0, xrng
#f0=max(Gauss(0),Lorenz(0))
#l0=Lorenz(0)
#gausspl,lorenzpl,hmpl,hm2pl=plot(rx, [Gauss(x-xmed)/f0 for x in rx], 'r', rx, [Lorenz(x-xmed)/f0 for x in rx], 'b', [xmed-hw/2.,xmed+hw/2.],[0.5,0.5], 'g', [xmed-hw/2.,xmed+hw/2.],[0.5*l0/f0,0.5*l0/f0], 'g')
#grid(True)
#ylabel('Response')
#xlabel('Frequency [THz]')

subplots_adjust(bottom=0.25, right=0.95)


axcolor = 'lightgoldenrodyellow'
axfwhmG = axes([0.1, 0.15, 0.6, 0.03], axisbg=axcolor)
fwhmG = Slider(axfwhmG, 'Gauss', 5*step, xrng/2, valinit=hw)
fwhmG.on_changed(updateG)

axfwhmL = axes([0.1, 0.05, 0.6, 0.03], axisbg=axcolor)
fwhmL = Slider(axfwhmL, 'Lorentz', 5*step, xrng/2, valinit=hw)
fwhmL.on_changed(updateL)

active=[True,False,True,True]
lpl.set_visible(False)

cbData=CheckButtons(axes([0.78, 0.05, 0.2, 0.15]),('Gauss','Lorenz','Theory','Experiment'),active)
cbData.on_clicked(changeData)


 
show()

for k in range(len(cx)) :
    print cx[k], cg[k], cl[k]

