# -*- coding: utf-8 -*-
"""
Base classes for the easy-to-understand and modify component-oriented acausal hybrid modeling
Differential-algebraic equations with Assimulo solvers.
Copyright © Volodymyr Kopei, 2017, 2019 email: vkopey@gmail.com"""

import numpy as np
from sympy import *

t=Symbol('t')

class Translational1D(object):
    """Base class of mechanical 1D components that have translational motion"""
    def __init__(self, name, args):
        self.name=name # component name
        for k,v in args.items(): # for each key-value pair
            if k in ['name','self']: continue # except name and self
            if v==None: # if value is None
                # create symbolic variable with name name+'_'+k
                self.__dict__[k]=Symbol(name+'_'+k)
            elif type(v) in [float,Float]: # if value is float
                self.__dict__[k]=Number(v) # create constant
        self.eqs=[] # equations list
        self.pins=[] # pins list
        
    def pinEqs(self,pindex,pins):
        eqs=[] # equation list of the flange
        f=Number(0) # sum of forces on flanges of other components
        for pin in pins: # for each flange of the other components
            # add equations describing the equality on the flange:
            eqs.append(Eq(self.pins[pindex]['x'], pin['x'])) # positions
            f+=pin['f'] # add to the sum of forces
        eqs.append(Eq(self.pins[pindex]['f'], -f)) # equality to zero the sum of forces on the flange 
        return eqs
    
class Mass(Translational1D):
    """Mass concentrated at a point, which has translational motion"""
    def __init__(self, name, m=1.0, x=None, v=None, a=None, f1=None, f2=None, Dx=None, Dv=None):
        Translational1D.__init__(self, name, locals()) # base class constructor call
        # system of equations
        self.eqs=[Eq(self.m*self.a, self.f1+self.f2),
                  Eq(self.a, self.Dv),
                  Eq(self.v, self.Dx)]
        self.pins=[dict(x=self.x, f=self.f1),
                   dict(x=self.x, f=self.f2)] # two flanges

class SpringDamper(Translational1D):
    """Translational 1D spring and damper, which are connected in parallel"""
    def __init__(self, name, c=1.0, d=0.1, x1=None, x2=None, vrel=None, f1=None, f2=None, Dx1=None, Dx2=None):
        Translational1D.__init__(self, name, locals())
        # system of equations
        self.eqs=[Eq(self.c*(self.x2-self.x1)+self.d*self.vrel, self.f2),
                  Eq(-self.f2, self.f1),
                  Eq(self.vrel, self.Dx2-self.Dx1)]
        
        self.pins=[dict(x=self.x1, f=self.f1),
                   dict(x=self.x2, f=self.f2)] # two flanges 

class Force(Translational1D):
    """1D force whose application point has translational motion"""
    def __init__(self,name,f=None, x=None):
        Translational1D.__init__(self, name, locals())
        self.pins=[dict(x=self.x, f=-self.f)] # one flange

class System(object):
    """System of components connected by flanges"""
    def __init__(self, els, eqs):
        self.els=els # components list
        self.elsd=dict([(e.name,e) for e in els]) # same, but dict.
        self.eqs=[] # list of system equations
        for e in self.els: # for each component
            self.eqs+=e.eqs # join with component equations
        self.eqs=self.eqs+eqs # join with additional equations
        self.eqs=Tuple(*self.eqs)
        
    def residualArgs(self,eq):
        "returns ordered arguments for residual (functions and derivatives of eq)"
        ss=eq.atoms(Symbol) # set of equation symbols
        ss.discard(t) # without t
        dss=dict([(i.name,i) for i in ss]) # dict name:symbol
        y=set();yd=set() # function; derivative
        for a in ss:
            if 'D' in a.name: yd.add(a)
            else: y.add(a)
        y_=[];yd_=[] # pairs function; derivative    
        for a in yd:
            b_name=a.name.replace('D','') # find pair
            if dss.get(b_name): # if pair
                yd_.append(a)
                y_.append(dss[b_name])
        yyd=set(y_+yd_)        
        y=y_+list(y-yyd)
        yd=yd_+list(yd-yyd)
        return y,yd
        
    def residual(self,t,y,yd): # residuals for Assimulo
        yyd=np.concatenate([[t],y,yd])[:self.nv]
        r=self.lambdfun(*yyd)
        return np.array(r)
            
    def solveDAE(self, eq, state, stopTime=10.0):
        """Solves dynamic task with Assimulo (ODASSL, IDA)
        state - dictionary with initial state"""
        y,yd=self.residualArgs(eq)
        self.y=y
        self.yd=yd
        self.nv=len(y+yd)+1 # number of arguments for lambdfun (with t)
        eq0=[e.rhs-e.lhs for e in eq]
        self.lambdfun=lambdify([t]+y+yd,eq0,'numpy')
        #import inspect
        #print(inspect.getsource(self.lambdfun))
        
        y0=[state[i] for i in y] # initial conditions
        yd0=[state[i] for i in yd]
        #provide the same length y0, yd0 (important for ODASSL):
        dn=len(y0)-len(yd0)
        if dn>0: yd0+=[0.0]*dn
        else: y0+=[0.0]*abs(dn)
        
        from assimulo.problem import Overdetermined_Problem,Implicit_Problem
        from assimulo.solvers import ODASSL,IDA
        
        # model = Overdetermined_Problem(self.residual, y0=y0, yd0=yd0)
        # sim = ODASSL(model)
        
        model = Implicit_Problem(self.residual, y0=y0, yd0=yd0)
        model.algvar = [1]*len(yd)+[0]*dn #[1,1,1,0,0,0,0,0] 
        sim = IDA(model)
        sim.suppress_alg = True
        print(sim.get_options())
        
        T, Y, Yd = sim.simulate(stopTime)
        #sim.plot()
        return T, Y, Yd 
        
    def solve(self,eq,ics): # for static tasks
        eq=eq.subs(ics)
        return solve(eq) # sympy solver

def prnt(eq): # eqations printing        
    print('\nEquations=')
    for i in eq: print(i)
    