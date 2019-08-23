# -*- coding: utf-8 -*-
"""Simulation of free vibrations of the sucker rod string
Trapezoidal rule.
Copyright © Volodymyr Kopei, 2017, 2019 email: vkopey@gmail.com"""
from pycodyn import *
from trapComponents import SpringDamper,Mass

# create components:
s1=SpringDamper(name='s1', c=44650.0, d=2120.0)
m1=Mass(name='m1', m=3961.0)
peqs=s1.pinEqs(1,[m1.pins[0]]) # list of additional equations
s=System(els=[s1,m1], eqs=peqs) # system
# solve the static problem — the column is stretched by 1 m
ics={m1.x:-1.0, m1.v:0.0, m1.a:0.0, s1.x1:0.0, s1.x1p:0.0, m1.vp:0.0, m1.ap:0.0, s1.v1:0.0, s1.v1p:0.0, s1.v2:0.0, s1.v2p:0.0}
d=s.solve(ics)

def fnBC(d, t):
    """boundary conditions at time t for fnBC.vrs components"""
    val = 0.0, 0.0 
    return dict(zip(fnBC.vrs, val))
fnBC.vrs = s.elsd['s1'].x1, s.elsd['m1'].f2

#solve the dynamic problem — free vibrations of the string
T,R=s.solveDyn(d, timeEnd=10, fnBC=fnBC)
plt.plot(T, [d[m1.x] for d in R])
plt.xlabel('t, s'); plt.ylabel('m1.x, m')
plt.show()
'''
import csv
with open('res.csv','wb') as f:
    w=csv.writer(f, delimiter=';')
    for r in zip(T, [d[m1.x] for d in R]):
        w.writerow(r)
'''