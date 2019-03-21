# -*- coding: utf-8 -*-
"""Base classes for the easy-to-understand and modify component-oriented acausal hybrid modeling (trapezoidal rule)
© Volodymyr Kopei, 2017, email: vkopey@gmail.com"""

from sympy import *
import math
import matplotlib
import matplotlib.pyplot as plt

dt=0.1
#dt=Symbol('dt')

class Translational1D(object):
    def __init__(self, name, args):
        self.name=name
        for k,v in args.iteritems():
            if k in ['name','self']: continue
            if v==None:
                self.__dict__[k]=Symbol(name+'_'+k)
            elif type(v) in [float,Float]:
                self.__dict__[k]=Number(v)
        self.eqs=[]
        self.pins=[]
        
    def pinEqs(self,pindex,pins):
        eqs=[] 
        f=Number(0) 
        for pin in pins: 
            
            eqs.append(Eq(self.pins[pindex]['x'], pin['x'])) 
            eqs.append(Eq(self.pins[pindex]['xp'], pin['xp'])) 
            f+=pin['f'] 
        eqs.append(Eq(self.pins[pindex]['f'], -f)) 
        return eqs
    
class Mass(Translational1D):
    def __init__(self,name,m=1.0,x=None,xp=None,v=None,vp=None,a=None,ap=None,f1=None,f2=None):
        Translational1D.__init__(self, name, locals()) 
        self.eqs=[Eq(self.m*self.a, self.f1+self.f2),
                  Eq((self.a+self.ap)/2, (self.v-self.vp)/dt),
                  Eq((self.v+self.vp)/2, (self.x-self.xp)/dt)] 
        self.pins=[dict(x=self.x, xp=self.xp, f=self.f1),
                   dict(x=self.x, xp=self.xp, f=self.f2)] 

class SpringDamper(Translational1D):
    def __init__(self,name,c=1.0,d=0.1,x1=None,x2=None,x1p=None,x2p=None,vrel=None,f1=None,f2=None,v1=None,v2=None,v1p=None,v2p=None):
        Translational1D.__init__(self, name, locals())
        self.eqs=[Eq(self.c*(self.x2-self.x1)+self.d*self.vrel, self.f2),
                  Eq(-self.f2, self.f1),
                  Eq((self.v1+self.v1p)/2, (self.x1-self.x1p)/dt),
                  Eq((self.v2+self.v2p)/2, (self.x2-self.x2p)/dt),
                  Eq(self.vrel, self.v2-self.v1)] 
        
        self.pins=[dict(x=self.x1, xp=self.x1p, f=self.f1),
                   dict(x=self.x2, xp=self.x2p, f=self.f2)] 

class Force(Translational1D):
    def __init__(self,name,f=None,x=None,xp=None):
        Translational1D.__init__(self, name, locals())
        self.pins=[dict(x=self.x, xp=self.xp, f=-self.f)] 

class System(object):
    def __init__(self, els, eqs):
        self.els=els 
        self.elsd=dict([(e.name,e) for e in els]) 
        self.eqs=[] 
        for e in self.els: 
            self.eqs+=e.eqs 
        self.eqs=self.eqs+eqs 
        
    def solveN(self, eqs):
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
                                           
    def solve(self, ics):
        eqs=[e.subs(ics) for e in self.eqs] 
        eqs=[e for e in eqs if e not in (True,False)] 
        
        sol=self.solveN(eqs) 
        sol.update(ics) 
        return sol    
    
    def solveDyn(self, d, timeEnd, fnBC):
        t=0.0 
        T=[] 
        Res=[] 
        ics={} 
        while t<timeEnd: 
            for e in self.els: 
                
                if 'x' in e.__dict__:
                    ics.update({e.xp:d[e.x]})
                if 'x1' in e.__dict__:
                    ics.update({e.x1p:d[e.x1]})
                if 'x2' in e.__dict__:
                    ics.update({e.x2p:d[e.x2]})
                if 'v' in e.__dict__:
                    ics.update({e.vp:d[e.v]})
                if 'a' in e.__dict__:
                    ics.update({e.ap:d[e.a]})
                if 'v1' in e.__dict__:
                    ics.update({e.v1p:d[e.v1]})
                if 'v2' in e.__dict__:
                    ics.update({e.v2p:d[e.v2]})
            ics.update(fnBC(self.elsd, d, t)) 
            d=self.solve(ics) 
            print t
            T.append(t)
            Res.append(d) 
            t+=dt 
        return T,Res
