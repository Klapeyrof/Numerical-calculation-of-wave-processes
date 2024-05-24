import numpy as np
import matplotlib.pyplot as plt
from sympy.interactive.printing import init_printing
init_printing(use_unicode=False, wrap_line=False)
from sympy import *

Sx = zeros(5, 5)
rho, lamda, mu = symbols('rho, lambda, mu')
Sx[0, 3] = -(lamda + 2*mu)/mu
Sx[0, 4] = 1
Sx[1, 1] = rho*sqrt(mu/rho)
Sx[1, 3] = 1
Sx[2, 1] = -rho*sqrt(mu/rho)
Sx[2, 3] = 1
Sx[3, 0] = rho*sqrt((lamda + 2*mu)/rho)
Sx[3, 2] = 1
Sx[4, 0] = -rho*sqrt((lamda + 2*mu)/rho)
Sx[4, 2] = 1

Sx = Sx.T

for i in range(5):
    print(Sx[i, :])

print()

Sx_inv = (Sx).inv(method="LU")
for i in range(5):
    print(Sx_inv[i, :])

print()


vx, vy, Gxx, Gyy, Gxy = symbols('vx, vy, Gxx, Gyy, Gxy')

q = zeros(5, 1)
q[0] = vx
q[1] = vy
q[2] = Gxx
q[3] = Gyy
q[4] = Gxy

print(q)

curr_w = Sx_inv*q
print(curr_w)

