from sympy import *

ρ, λ, μ = symbols('ρ λ μ')


Ax = zeros(5, 5)
Ax[0, 2] = -(λ + 2*μ)
Ax[1, 3] = -μ
Ax[1, 4] = -(λ + 2*μ)
Ax[2, 0] = -1/ρ
Ax[3, 1] = -1/ρ

for i in range(5):
    print(Ax[i, :])
print("Собственные значения: ", Ax.eigenvals(), '\n')
print("Собственные векторы: ", Ax.eigenvects())

print('\n'*5)

Ay = zeros(5, 5)
Ay[0, 2] = -(λ + 2*μ)
Ay[0, 3] = -μ
Ay[1, 4] = -(λ + 2*μ)
Ay[3, 0] = -1/ρ
Ay[4, 1] = -1/ρ

for i in range(5):
    print(Ay[i, :])
print("Собственные значения: ", Ay.eigenvals(), '\n')
print("Собственные векторы: ", Ay.eigenvects())


