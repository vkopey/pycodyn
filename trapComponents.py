# -*- coding: utf-8 -*-
"""Components for Trapezoidal rule.
Copyright © Volodymyr Kopei, 2017, 2019 email: vkopey@gmail.com"""

from pycodyn import *

class Mass(Translational1D):
    def __init__(self,name,m=1.0,x=None,xp=None,v=None,vp=None,a=None,ap=None,f1=None,f2=None):
        Translational1D.__init__(self, name, locals()) # виклик конструктора базового класу
        
        self.eqs=[Eq(self.m*self.a, self.f1+self.f2),
                  Eq((self.a+self.ap)/2, (self.v-self.vp)/dt),
                  Eq((self.v+self.vp)/2, (self.x-self.xp)/dt)] # система рівнянь
        self.pins=[dict(x=self.x, xp=self.xp, f=self.f1),
                   dict(x=self.x, xp=self.xp, f=self.f2)] # два фланця

class SpringDamper(Translational1D):
    def __init__(self,name,c=1.0,d=0.1,x1=None,x2=None,x1p=None,x2p=None,vrel=None,f1=None,f2=None,v1=None,v2=None,v1p=None,v2p=None):
        Translational1D.__init__(self, name, locals())
        
        self.eqs=[Eq(self.c*(self.x2-self.x1)+self.d*self.vrel, self.f2),
                  Eq(-self.f2, self.f1),
                  Eq((self.v1+self.v1p)/2, (self.x1-self.x1p)/dt),
                  Eq((self.v2+self.v2p)/2, (self.x2-self.x2p)/dt),
                  Eq(self.vrel, self.v2-self.v1)] # система рівнянь
        
        self.pins=[dict(x=self.x1, xp=self.x1p, f=self.f1),
                   dict(x=self.x2, xp=self.x2p, f=self.f2)] # два фланця 
