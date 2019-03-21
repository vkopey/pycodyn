# -*- coding: utf-8 -*-
"""Base classes for the easy-to-understand and modify component-oriented acausal hybrid modeling
© Volodymyr Kopei, 2017, email: vkopey@gmail.com"""

from sympy import *
import math
import matplotlib
import matplotlib.pyplot as plt

dt=0.1 # time step
#dt=Symbol('dt') # only to obtain equations in a symbolic form

class Translational1D(object):
    """Base class of mechanical 1D components that have translational motion"""
    def __init__(self, name, args):
        self.name=name # component name
        for k,v in args.iteritems(): # for each key-value pair
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
            eqs.append(Eq(self.pins[pindex]['xp'], pin['xp'])) # positions at time t-dt
            f+=pin['f'] # add to the sum of forces
        eqs.append(Eq(self.pins[pindex]['f'], -f)) # equality to zero the sum of forces on the flange 
        return eqs
    
class Mass(Translational1D):
    """Mass concentrated at a point, which has translational motion"""
    def __init__(self, name, m=1.0, x=None, xp=None, v=None, vp=None, a=None, f1=None, f2=None):
        Translational1D.__init__(self, name, locals()) # base class constructor call
        # system of equations
        self.eqs=[Eq(self.m*self.a, self.f1+self.f2),
                  Eq(self.a, (self.v-self.vp)/dt),
                  Eq(self.v, (self.x-self.xp)/dt)]
        self.pins=[dict(x=self.x, xp=self.xp, f=self.f1),
                   dict(x=self.x, xp=self.xp, f=self.f2)] # two flanges

class SpringDamper(Translational1D):
    """Translational 1D spring and damper, which are connected in parallel"""
    def __init__(self, name, c=1.0, d=0.1, x1=None, x2=None, x1p=None, x2p=None, vrel=None, f1=None, f2=None):
        Translational1D.__init__(self, name, locals())
        # system of equations
        self.eqs=[Eq(self.c*(self.x2-self.x1)+self.d*self.vrel, self.f2),
                  Eq(-self.f2, self.f1),
                  Eq(self.vrel, (self.x2-self.x2p)/dt-(self.x1-self.x1p)/dt)]
        
        self.pins=[dict(x=self.x1, xp=self.x1p, f=self.f1),
                   dict(x=self.x2, xp=self.x2p, f=self.f2)] # two flanges 

class Force(Translational1D):
    """1D force whose application point has translational motion"""
    def __init__(self,name,f=None,x=None,xp=None):
        Translational1D.__init__(self, name, locals())
        self.pins=[dict(x=self.x, xp=self.xp, f=-self.f)] # one flange

class System(object):
    """System of components connected by flanges"""
    def __init__(self, els, eqs):
        self.els=els # components list
        self.elsd=dict([(e.name,e) for e in els]) # same, but dict.
        self.eqs=[] # list of system equations
        for e in self.els: # for each component
            self.eqs+=e.eqs # join with component equations
        self.eqs=self.eqs+eqs # join with additional equations
        
    def solveN(self, eqs): # solves the static problem
        import scipy.optimize
        eqs0=[]
        vrs=set()
        for e in eqs:
            eq=e.lhs-e.rhs
            eqs0.append(eq)
            for a in eq.atoms():
                if a.is_Symbol:
                    vrs.add(a)
        vrs=list(vrs)
        f=lambdify([vrs], eqs0, 'numpy')
        goals=[0.0 for i in vrs]
        sol=scipy.optimize.root(f, goals, method='lm') 
        d=dict(zip(vrs,sol.x))
        return d
                                           
    def solve(self, ics): # solves the dynamic problem
        eqs=[e.subs(ics) for e in self.eqs] # substitution of ics
        eqs=[e for e in eqs if e not in (True,False)] # discard all degenerate equations
        # solve the system of equations by:
        #sol=solve(eqs) # SymPy (slow)
        sol=self.solveN(eqs) # SciPy (faster)
        sol.update(ics) # update dictionary by dictionary ics
        return sol    
    
    def solveDyn(self, d, timeEnd, fnBC):
        t=0.0 # time variable
        T=[] # list of time values
        Res=[] # list of results
        ics={} # dictionary with values of variables
        while t<timeEnd: # while t < final time value
            for e in self.els: # for each component
                # save positions and velocities
                if 'x' in e.__dict__:
                    ics.update({e.xp:d[e.x]})
                if 'x1' in e.__dict__:
                    ics.update({e.x1p:d[e.x1]})
                if 'x2' in e.__dict__:
                    ics.update({e.x2p:d[e.x2]})
                if 'v' in e.__dict__:
                    ics.update({e.vp:d[e.v]})
            ics.update(fnBC(self.elsd, d, t)) # update BC
            d=self.solve(ics) # solve the problem
            print t
            T.append(t)
            Res.append(d) # save results
            t+=dt # increase time value
            #if some_condition: # changing the system structure
            #   self.__init__(new_els, new_eqs) 
        return T,Res
