from sympy import *
def eliminate(eq, keep, maxEqLen=2):
    """Eliminate variables from equations eq
    keep - try to keep these variables 
    maxEqLen - max number of equations after elimination"""
    while len(eq)>maxEqLen:
        e=solve(eq, eq.free_symbols-keep, exclude=keep) # for elimination
        e=sorted(e.items(), key=lambda x:len(x[1].atoms())) # first shortest expressions
        eq=eq.subs(e[0][0], e[0][1]) # eliminate
        eq=Tuple(*set(eq)-set([True])) # without dublicates and True
    return eq