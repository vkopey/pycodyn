# encoding: utf-8
"""Simulation of the pumping process by 1-section string
Euler method. 
[s1]---[m1]-+
            |
           [f1]
Copyright © Volodymyr Kopei, 2017, 2019 email: vkopey@gmail.com"""
from pycodyn import *

fs=-34687.0 # section weight
fr=-18499.0 # liquid weight above the plunger          
# components:
s1=SpringDamper(name='s1', c=44650.0, d=2120.7)
m1=Mass(name='m1', m=3961.0)
f1=Force(name='f1')#, f=fs+fr 

# additional equations of the string model, formed by connecting of the components flanges
peqs=s1.pinEqs(1,[m1.pins[0]])
peqs+=m1.pinEqs(1,[f1.pins[0]])
s=System(els=[s1,m1,f1], eqs=peqs) # system

# static problem — the string under the maximum static loads
ics={m1.v:0.0, m1.a:0.0, s1.x1:0.0, s1.x1p:0.0, f1.f:fs+fr}
d=s.solve(ics)
print(d[m1.x])

# for k in sorted(d.keys(), key=lambda k:repr(k)):
#     print(k,d[k])

def motion(t):
    """describes the harmonic motion of the upper point and returns its position at time t"""
    A=2.1/2 # amplitude
    n=6.4/60 # frequency
    return A*math.sin(2*math.pi*n*t) # position
   
def force(v):
    """returns the value of the force on the pump plunger F, depending on the value of its speed v"""
    F=fs # weight of the section
    if v>0: # if upperstroke
        F+=fr # increase the force by value of the fluid weight
    return F*math.tanh(abs(v)/0.01) # smoothing near the point v=0

def fnBC(d, t):
    """boundary conditions at time t for fnBC.vrs components"""
    val = motion(t), force(d[m1.v]) 
    return dict(zip(fnBC.vrs, val))
fnBC.vrs = s.elsd['s1'].x1, s.elsd['f1'].f

# solve the dynamic problem — the upper point has a harmonic motion
T,R=s.solveDyn(d, timeEnd=2*60/6.4, fnBC=fnBC)
R=[r for t,r in zip(T,R) if t>60/6.4] # only last period
plt.plot([d[s1.x1] for d in R], [d[s1.f1]/1000 for d in R]) # wellhead dynamometer card
plt.plot([d[m1.x] for d in R], [(-d[m1.f2]+fs)/1000 for d in R]) # plunger dynamometer card
plt.xlabel('x, m'); plt.ylabel('f, kN')
plt.show()
'''
import csv
with open('res.csv','wb') as f:
    w=csv.writer(f, delimiter=';')
    for r in zip([d[s1.x1] for d in R], [d[s1.f1]/1000 for d in R]):
        w.writerow(r)
'''