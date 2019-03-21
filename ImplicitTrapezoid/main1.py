# -*- coding: utf-8 -*-
"""Simulation of free vibrations of the sucker rod string (trapezoidal rule)
© Volodymyr Kopei, 2017, email: vkopey@gmail.com"""
from pycodyn import *

s1=SpringDamper(name='s1', c=39694.0, d=1856.0)
m1=Mass(name='m1',m=3402.0)
peqs=s1.pinEqs(1,[m1.pins[0]])
s=System(els=[s1,m1], eqs=peqs)
ics={m1.x:-1.0, m1.v:0.0, m1.a:0.0, s1.x1:0.0, s1.x1p:0.0, m1.vp:0.0, m1.ap:0.0, s1.v1:0.0, s1.v1p:0.0, s1.v2:0.0, s1.v2p:0.0}
d=s.solve(ics)

def fnBC(elsd, d, t):
    return {elsd['s1'].x1:0.0, elsd['s1'].x1p:0.0, elsd['m1'].f2:0.0}

T,R=s.solveDyn(d, timeEnd=10, fnBC=fnBC)
plt.plot(T, [d[m1.x] for d in R])
plt.xlabel('t, s'); plt.ylabel('m1.x, m')
plt.show()