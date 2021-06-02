import sympy as sym
from sympy.solvers import solve
from sympy import Eq
import pandas as pd
import os
os.chdir(os.getcwd())

CHW = pd.read_csv("Supernova_CHW1.csv")
D_i = list(CHW["D"])
V_i = list(CHW["V"])

H0 = sym.Symbol('H0')
H1 = sym.Symbol('H1')

X1 = sum(list(map(lambda D, V: D**2 + V**2*(1/H0**2) -2*D*V*(1/H0), D_i, V_i)))
dX1_dH0 = sym.diff(X1, H0)

X2 = sum(list(map(lambda D, V: D**2 + (V**4)*(1/H1**2) + (V**2)*(1/H0**2) - (2*D*(V**2))*(1/H1) - 2*D*V*(1/H0) + (2*(V**3))*(1/(H0*H1)), D_i, V_i)))
dX2_dH0 = sym.diff(X2, H0)
dX2_dH1 = sym.diff(X2, H1)

print(f'X1 penalty function is : {X1}')
print(f'H0 is : {solve(dX1_dH0, H0)}')
print(f'down limit and up limit of H0 is :{solve(X1-1, H0)}')

eq1 = dX2_dH0
eq2 = dX2_dH1

print(f'X2 penalty function is : {X2}')
print(f'H0 & H1 are : {solve((eq1, eq2), (H0, H1))}')