# -*- coding: utf-8 -*-
"""Simulation of free vibrations of the sucker rod string
Symbolic ODE-solver SymPy.
Copyright © Volodymyr Kopei, 2017, 2019 email: vkopey@gmail.com"""

from pycodynDAE import *

# create components:
s1=SpringDamper(name='s1', c=44650.0, d=2120.0)
m1=Mass(name='m1', m=3961.0)
peqs=s1.pinEqs(1,[m1.pins[0]]) # list of additional equations
s=System(els=[s1,m1], eqs=peqs) # system
prnt(s.eqs)
eq=s.eqs.subs({m1.f2:0.0, s1.x1:0.0, s1.Dx1:0.0}) # BC
prnt(eq)

# eliminate variables manually:
eq=eq.subs({s1.f2:-s1.f1, m1.f1:s1.f1, s1.x2:m1.x, s1.Dx2:m1.Dx, s1.vrel:m1.v, m1.a:m1.Dv})
eq=eq.subs(s1.f1, solve(eq[3], s1.f1)[0])
eq=Tuple(*set(eq)-{True})

# or automatically:
# from eliminate import *
# eq=eliminate(eq, keep={m1.Dv,m1.Dx,m1.v,m1.x}, maxEqLen=2)
prnt(eq)

# substitution of functions and derivatives
eq=eq.subs({m1.x:Function('m1_x')(t), m1.v:Function('m1_v')(t)})
eq=eq.subs({m1.Dx:Derivative('m1_x(t)', t), m1.Dv:Derivative('m1_v(t)', t)})
prnt(eq)
eq=dsolve(eq, ics={Function('m1_x')(0):-1.0, Function('m1_v')(0):0.0}) # solve ODE symbolically
print(eq)
plot(eq[1].rhs, xlim=(0,10), ylim=(-1,1))
'''
import numpy as np
T=np.linspace(0,10)
#f=lambdify(t,eq[1].rhs)
import csv
with open('res.csv','wb') as f:
    w=csv.writer(f, delimiter=';')
    for r in [(ti, eq[1].rhs.subs(t,ti)) for ti in T]:
        w.writerow(r)
'''