# encoding: utf-8
"""Simulation of the variable structure system (breakage of the second section) by event handling.
[s1]---[m1]-+-[s2]---[m2]-+
            |             |
           [f1]          [f2]
Copyright © Volodymyr Kopei, 2017, 2019 email: vkopey@gmail.com"""
from pycodyn import *
        
def event(self, state): # event handler
    # simulation of the breakage of the second section when force>56000
    if state[s1.f1]>56000:
        self.fnBC=fnBC2
        peqs=s1.pinEqs(1,[m1.pins[0]])
        peqs+=m1.pinEqs(1,[f1.pins[0]])
        self.__init__(els=[s1,m1,f1], eqs=peqs) # changing the system             
        self.createCurEqs(fnBC2) # new current equations
System.event=event            
                
fs=(-18494.0, -16193.0) # sections weights
fr=-18499.0 # liquid weight          
# components:
s1=SpringDamper(name='s1', c=114926.0, d=5458.0)
m1=Mass(name='m1',m=2112.0)
f1=Force(name='f1', f=fs[0]) 
s2=SpringDamper(name='s2', c=73021.0, d=3468.0)
m2=Mass(name='m2', m=1850.0)
f2=Force(name='f2')
# additional equations:
peqs=s1.pinEqs(1,[m1.pins[0]])
peqs+=m1.pinEqs(1,[s2.pins[0],f1.pins[0]])
peqs+=s2.pinEqs(1,[m2.pins[0]])
peqs+=m2.pinEqs(1,[f2.pins[0]])
s=System(els=[s1,m1,s2,m2,f1,f2], eqs=peqs) # system

# static problem — the string under the maximum static loads
ics={m1.v:0.0, m1.a:0.0, m2.v:0.0, m2.a:0.0, s1.x1:0.0, s1.x1p:0.0, f2.f:fs[1]+fr}
d=s.solve(ics)
print(d[m2.x])

def motion(t):
    """describes the harmonic motion of the upper point and returns its position at time t"""
    A=2.1/2 # amplitude
    n=6.4/60 # frequency
    return A*math.sin(2*math.pi*n*t) # position
   
def force(v):
    """returns the value of the force on the pump plunger F, depending on the value of its speed v"""
    F=fs[1] # weight of the second section
    if v>0: # if upperstroke
        F+=fr # increase the force by value of the fluid weight
    return F*math.tanh(abs(v)/0.01) # smoothing near the point v=0

def fnBC(d, t):
    """boundary conditions at time t for fnBC.vrs components"""
    val = motion(t), force(d[m2.v]) 
    return dict(zip(fnBC.vrs, val))
fnBC.vrs = s.elsd['s1'].x1, s.elsd['f2'].f

def fnBC2(d, t):
    """boundary conditions 2 at time t for fnBC2.vrs components"""
    val = (motion(t), )
    return dict(zip(fnBC.vrs, val))
fnBC2.vrs = (s.elsd['s1'].x1, )

# solve the dynamic problem — the upper point has a harmonic motion
T,R=s.solveDyn(d, timeEnd=2*60/6.4+10, fnBC=fnBC)
plt.plot([d[s1.x1] for d in R], [d[s1.f1]/1000 for d in R]) # wellhead dynamometer card
plt.xlabel('x, m'); plt.ylabel('f, kN')
plt.show()