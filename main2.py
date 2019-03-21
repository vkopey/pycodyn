# encoding: utf-8
"""Simulation of the pumping process by two-section string  
[s1]---[m1]-+-[s2]---[m2]-+
            |             |
           [f1]          [f2]
© Volodymyr Kopei, 2017, email: vkopey@gmail.com"""
from pycodyn import *

fs=(-14602.0, -14602.0) # sections weights
fr=-16688.0 # liquid weight above the plunger          
# components:
s1=SpringDamper(name='s1', c=79388.0, d=3712.0)
m1=Mass(name='m1',m=1701.0)
f1=Force(name='f1', f=fs[0]) 
s2=SpringDamper(name='s2', c=79388.0, d=3712.0)
m2=Mass(name='m2', m=1701.0)
f2=Force(name='f2', f=fs[1]+fr)
# additional equations of the string model, formed by connecting of the components flanges
peqs=s1.pinEqs(1,[m1.pins[0]])
peqs+=m1.pinEqs(1,[s2.pins[0],f1.pins[0]])
peqs+=s2.pinEqs(1,[m2.pins[0]])
peqs+=m2.pinEqs(1,[f2.pins[0]])
s=System(els=[s1,m1,s2,m2,f1,f2], eqs=peqs) # system
#print s.eqs

# static problem — the string under the maximum static loads
ics={m1.v:0.0, m1.a:0.0, m2.v:0.0, m2.a:0.0} 
ics.update({s1.x1:0.0, s1.x1p:0.0})
d=s.solve(ics)
print d[s1.x1], d[m2.x]

# print d
# for k in sorted(d.keys(), key=lambda k:repr(k)):
#     print k,d[k]

def motion(t):
    """describes the harmonic motion of the upper point and returns its position at time t"""
    A=3.0/2 # amplitude
    n=6.5/60 # frequency
    return A*math.sin(2*math.pi*n*t) # position
   
def force(v):
    """returns the value of the force on the pump plunger F, depending on the value of its speed v"""
    F=fs[1] # weight of the second section
    if v>0: # if upperstroke
        F+=fr # increase the force by value of the fluid weight
    return F*math.tanh(abs(v)/0.01) # smoothing near the point v=0

def fnBC(elsd, d, t):
    """boundary conditions at time t for elsd components"""
    return {elsd['s1'].x1:motion(t), elsd['f2'].f:force(d[m2.v])}

# solve the dynamic problem — the upper point has a harmonic motion
T,R=s.solveDyn(d, timeEnd=2*60/6.5, fnBC=fnBC)
R=[r for t,r in zip(T,R) if t>60/6.5] # only last period
plt.plot([d[s1.x1] for d in R], [d[s1.f1]/1000 for d in R]) # wellhead dynamometer card
plt.plot([d[m2.x] for d in R], [(-d[m2.f2]+fs[1])/1000 for d in R]) # plunger dynamometer card
plt.xlabel('x, m'); plt.ylabel('f, kN')
plt.show()