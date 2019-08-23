# -*- coding: utf-8 -*-
"""Base classes for the easy-to-understand and modify component-oriented acausal hybrid modeling.
Difference equations with Euler method.
Copyright © Volodymyr Kopei, 2017, 2019 email: vkopey@gmail.com"""

from sympy import *
import math, time
import matplotlib.pyplot as plt

def byName(d,name): # return value by symbol name
    for k in d:
        if repr(k)==name: return d[k]

dt=0.1 # time step
#dt=Symbol('dt') # only to obtain equations in a symbolic form

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
        
    def solveN(self, eqs): # solve alg. system by scipy
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
                                               
    def solve(self, ics): # solve alg. system at t
        eqs=[e.subs(ics) for e in self.eqs] # substitution of ics
        eqs=[e for e in eqs if e not in (True,False)] # discard all degenerate equations
        # solve the system of equations by:
        #sol=nsolve(eqs) # SymPy (slow) #or solve
        sol=self.solveN(eqs) # SciPy (faster)
        sol.update(ics) # update dictionary by dictionary ics
        return sol
        
    def solv(self, preState): # solve by subs. to sympy expr.
        state=preState.copy()
        for k in self.ceqs: # current equations
            state[k]=self.ceqs[k].subs(preState).evalf()
            #assert type(state[k]) in [float,Float]
        return state
        
    def solvN(self, preState): # same but by lambdafunction
        #use Python 3.7 for fastest execution
        state=preState.copy()
        ls=dict([(repr(a),state[a]) for a in self.vrsp]) # arg. dict
        res=self.ceqsf(**ls) # lambdafunction call
        for a,v in zip([i[0] for i in self.ceqsi], res):
            state[a]=v # update state
        return state
        
    def createCurEqs(self, fnBC): # create current 'fast equations'
        eqs=Tuple(*self.eqs)
        vrs={i for i in eqs.atoms(Symbol) if repr(i)[-1]!='p'} # vars without 'p'
        vrsbc=set(fnBC.vrs)
        vrs=vrs-vrsbc # unknown vars at current step
        self.ceqs=solve(eqs,vrs) # current expressions
        self.vrsp={i for i in eqs.atoms(Symbol) if repr(i)[-1]=='p'} # vars with 'p'
        self.vrsp.update(vrsbc) # known vars at current step
        self.ceqsi=self.ceqs.items() # ordered expressions
        self.ceqsf=lambdify(self.vrsp,[i[1] for i in self.ceqsi],'numpy') # current lambda function
                         
    def solveDyn(self, state, timeEnd, fnBC):
        #solves the dynamic problem 
        t=0.0 # time variable
        T=[] # list of time values
        Res=[] # list of results
        
        self.createCurEqs(fnBC)
        ics={} # for self.solve()
        start = time.time()
        while t<timeEnd:
            for k in state: # previous values "xp=x"...
                if repr(k)[-1]=='p':
                    state[k]=byName(state,repr(k)[:-1])
                    ics[k]=byName(state,repr(k)[:-1]) # for self.solve()
            
            state.update(fnBC(state, t)) # update BC
            #state=self.solv(state) # by sympy expression (slow)
            state=self.solvN(state) # by numpy expression (fast)
            
            # ics.update(fnBC(state, t)) # update BC
            # state=self.solve(ics) # by scipy.optimize.root (slow, for frequent events)
            
            print(t)
            T.append(t)
            Res.append(state) # save results
            t+=dt # increase time value
            
            self.event(state) # event handler
        end = time.time()
        print('simulation time',end-start)
        return T,Res
        
    def event(self, state): # event handler
        pass
