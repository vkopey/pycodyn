# -*- coding: utf-8 -*-
"""Simulation of the pumping process by two-section string
Assimulo. 
[s1]---[m1]-+-[s2]---[m2]-+
            |             |
           [f1]          [f2]
Copyright © Volodymyr Kopei, 2017, 2019 email: vkopey@gmail.com"""

from pycodynDAE import *

fs=(-18494.0, -16193.0) # sections weights
fr=-18499.0 # liquid weight above the plunger          
# components:
s1=SpringDamper(name='s1', c=114926.0, d=5458.0)
m1=Mass(name='m1', m=2112.0)
f1=Force(name='f1', f=fs[0]) 
s2=SpringDamper(name='s2', c=73021.0, d=3468.0)
m2=Mass(name='m2', m=1850.0)
f2=Force(name='f2') #, f=fs[1]+fr
# additional equations of the string model, formed by connecting of the components flanges
peqs=s1.pinEqs(1,[m1.pins[0]])
peqs+=m1.pinEqs(1,[s2.pins[0],f1.pins[0]])
peqs+=s2.pinEqs(1,[m2.pins[0]])
peqs+=m2.pinEqs(1,[f2.pins[0]])
s=System(els=[s1,m1,s2,m2,f1,f2], eqs=peqs) # system
prnt(s.eqs)

# static problem — the string under the maximum static loads
ics={s1.x1:0.0, s1.Dx1:0.0, f2.f:fs[1]+fr, m1.v:0.0, m1.a:0.0, s1.Dx2:0.0, s2.Dx1:0.0, s2.Dx2:0.0, m2.v:0.0, m2.a:0.0}
state=s.solve(s.eqs,ics)
state.update(ics)

# dynamic problem — the upper point has a harmonic motion
A=2.1/2 # amplitude
n=6.4/60 # frequency
eq=s.eqs.subs({s1.x1: A*sin(2*pi*n*t), s1.Dx1: A*sin(2*pi*n*t).diff(t), m2.f2: Piecewise((fs[1], m2.v<0), (fs[1]+fr*tanh(abs(m2.v)/0.01), m2.v>=0))})
#or add equations:
#eq=s.eqs+Tuple(Eq(s1.x1, A*sin(2*pi*n*t)), Eq(m1.f2, Piecewise((fs, m1.v<0), (fs+fr*tanh(abs(m1.v)/0.01), m1.v>=0))) )

T,Y,Yd=s.solveDAE(eq, state, 2*60/6.4)

import matplotlib.pyplot as plt
index=s.y.index(s1.f1)
Y=[y for t,y in zip(T,Y) if t>60/6.4] # only last period
T=[t for t in T if t>60/6.4]
plt.plot(A*np.sin(2*np.pi*n*np.array(T)), [v[index]/1000 for v in Y])
plt.show()
'''
import csv
with open('res.csv','wb') as f:
    w=csv.writer(f, delimiter=';')
    for r in zip(A*np.sin(2*np.pi*n*np.array(T)), [v[index]/1000 for v in Y]):
        w.writerow(r)
'''