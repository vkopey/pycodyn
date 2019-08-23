# -*- coding: utf-8 -*-
"""Simulation of the pumping process by 1-section string
Assimulo. 
[s1]---[m1]-+
            |
           [f1]
Copyright © Volodymyr Kopei, 2017, 2019 email: vkopey@gmail.com"""

from pycodynDAE import *

fs=-34687.0 # sections weights
fr=-18499.0 # liquid weight above the plunger
# components:
s1=SpringDamper(name='s1', c=44650.0, d=2120.7)
m1=Mass(name='m1', m=3961.0)
peqs=s1.pinEqs(1,[m1.pins[0]]) # list of additional equations
s=System(els=[s1,m1], eqs=peqs) # system
prnt(s.eqs)

# static problem — the string under the maximum static loads
ics={s1.x1:0.0, s1.Dx1:0.0, m1.f2:fs+fr, m1.v:0.0, m1.a:0.0, s1.Dx2:0.0}
state=s.solve(s.eqs,ics)
state.update(ics)

# dynamic problem — the upper point has a harmonic motion
A=2.1/2 # amplitude
n=6.4/60 # frequency
eq=s.eqs.subs({s1.x1: A*sin(2*pi*n*t), s1.Dx1: A*sin(2*pi*n*t).diff(t), m1.f2: Piecewise((fs, m1.v<0), (fs+fr*tanh(abs(m1.v)/0.01), m1.v>=0))})
#or add equations:
#eq=s.eqs+Tuple(Eq(s1.x1, A*sin(2*pi*n*t)), Eq(m1.f2, Piecewise((fs, m1.v<0), (fs+fr*tanh(abs(m1.v)/0.01), m1.v>=0))) )
T,Y,Yd=s.solveDAE(eq, state, 20.0)

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