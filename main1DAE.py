# -*- coding: utf-8 -*-
"""Simulation of free vibrations of the sucker rod string
Assimulo.
Copyright © Volodymyr Kopei, 2017, 2019 email: vkopey@gmail.com"""
from pycodynDAE import *

# create components:
s1=SpringDamper(name='s1', c=44650.0, d=2120.0)
m1=Mass(name='m1',m=3961.0)
peqs=s1.pinEqs(1,[m1.pins[0]]) # list of additional equations
s=System(els=[s1,m1], eqs=peqs) # system
prnt(s.eqs)

bc={s1.x1:0.0, s1.Dx1:0.0}
eq=s.eqs.subs(bc) # constant boundary conditions
prnt(eq)

#static — the column is stretched by 1 m
ics={m1.x:-1.0, m1.v:0.0, m1.a:0.0, s1.Dx2:0.0}
state=s.solve(eq,ics)
state.update(ics)

#dynamic — free vibrations of the string
eq=eq.subs({m1.f2:0.0}) # additional BC
T,Y,Yd=s.solveDAE(eq, state, 10.0)

import matplotlib.pyplot as plt
index=s.y.index(m1.x)
plt.plot(T, [v[index] for v in Y])
plt.show()
'''
import csv
with open('res.csv','wb') as f:
    w=csv.writer(f, delimiter=';')
    for r in zip(T, [v[index] for v in Y]):
        w.writerow(r)
'''