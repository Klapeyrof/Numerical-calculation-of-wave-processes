from matplotlib import pyplot as plt
import matplotlib.animation as animation
from math import pi, sin, sqrt
import time



# 1. Задание функций (матриц и граничных условий)

def W(v, p, i):
    w1 = -c[i]*rho[i]/2*v + 1/2*p
    w2 = c[i]*rho[i]/2*v + 1/2*p
    return w1, w2

def SolveDiscontinuity(w1, w2, i):
    v = -2.0 * (w1 - w2) / (c[i - 1] * rho[i - 1] + c[i] * rho[i])
    p = 2.0 * (c[i - 1] * rho[i - 1] * w1 + c[i] * rho[i] * w2) / (c[i - 1] * rho[i - 1] + c[i] * rho[i])
    return (v, p)


# Задаёт функцию на левой границе (aP + bV = f_left(t))
def f_left(tay):
    return 0 #sin(tay)

# Задаёт давление/скорость на левой стенке
def Border_left(w1, alpha, beta, tay):
    return (f_left(tay) - w1 * (alpha - beta / (rho[0] * c[0]))) / ((alpha + beta / (rho[0] * c[0])))

# Задаёт функцию на правой границе (aP + bV = f_right(t))
def f_right(tay):
    return 0

# Задаёт давление/скорость на левой стенке
def Border_right(w2, alpha, beta, tay):
    return (f_right(tay) - w2 * (alpha + beta / (rho[N-1] * c[N-1]))) / ((alpha - beta / (rho[N-1] * c[N-1])))



# 2. Задание параметров

N = 300
rho = []
c = []
for i in range(int(N/2)):
    rho.append(1000)
    c.append(1000)
for i in range(int(N/2), N):
    rho.append(10)
    c.append(300)

L = 1000
t = 0.5 * L / max(c)
hx = L/N
Ku = 0.7
dt = Ku*hx/max(c)
steps = int(t/dt)


# 3. Задание массивов данных и начальных условий

x = []
for i in range(N):
    x.append(hx*i)

data0_p = [0] * N
data0_v = [0] * N
curr_v = [0] * N
curr_p = [0] * N
next_v = [0] * N
next_p = [0] * N

for i in range(int(N/3)):
    data0_v[i] = 1/c[0]
    data0_p[i] = 1*rho[0]
    curr_v[i] = data0_v[i]
    curr_p[i] = data0_p[i]


# Весь массив данных
V_All_data = []
P_All_data = []

for k in range(steps + 1):
    V_All_data.append([])
    P_All_data.append([])
    for i in range(N):
        V_All_data[k].append(0)
        P_All_data[k].append(0)

V_All_data[0] = data0_v
P_All_data[0] = data0_p






# 4. Сам код


# w1_next - значение на n шаге в i+1 точке
# w2_prev - значение на n шаге в i-1 точке
# w1_new - значение на n+1 шаге в i-ой точке
# next => R, prev => L. (w1 справа, w2 слева)
for k in range(1, steps):
    for i in range(N):

        # Считываем граничное условие на стенках
        #P_Border_left, V_Border_left = Border_left(j * dt)
        #P_Border_right, V_Border_right = Border_right(j * dt)

        if(i == 0):
            w1_i, _ = W(curr_v[i], curr_p[i], i)                    # Считаем честно из матрицы
            w1_next, _ = W(curr_v[i+1], curr_p[i+1], i)             # Значения, учитывающие материал справа
            w1_new = w1_i - c[i] * dt / hx * (w1_i - w1_next)       # Честно посчитали т.к. справа знаем данные
            w2_new = Border_left(w1_new, 1, 0, k*dt)                # Из граничного условия aP + bV = f_left(t) на каждом шаге

        if(i > 0 and i < N-1):
            w1_i, _ = W(curr_v[i], curr_p[i], i)                    # Считаем честно из матрицы
            _, w2_i = W(curr_v[i], curr_p[i], i - 1)                # Считаем честно из матрицы

            _, w2_prev = W(curr_v[i-1], curr_p[i-1], i-1)           # Значения, учитывающие материал слева
            w1_next, _ = W(curr_v[i+1], curr_p[i+1], i)             # Значения, учитывающие материал справа

            w1_new = w1_i - c[i] * dt / hx * (w1_i - w1_next)
            w2_new = w2_i - c[i-1] * dt / hx * (w2_i - w2_prev)

        if (i == N - 1):
            _, w2_i = W(curr_v[i], curr_p[i], i-1)
            _, w2_prev = W(curr_v[i-1], curr_p[i-1], i-1)
            w2_new = w2_i - c[i] * dt / hx * (w2_i - w2_prev)
            w1_new = Border_right(w2_new, 1, 0, k * dt)              # Из граничного условия aP + bV = f_right(t) на каждом шаге

        vi_next, pi_next = SolveDiscontinuity(w1_new, w2_new, i)
        next_v[i] = vi_next
        next_p[i] = pi_next
    for i in range(N):
        curr_v[i] = next_v[i]
        curr_p[i] = next_p[i]
    for i in range(N):
        V_All_data[k][i] = curr_v[i]
        P_All_data[k][i] = curr_p[i]


print("Интерполяция завершена")





# 5. Отображение анимации

fig_V, ax_V = plt.subplots()
fig_P, ax_P = plt.subplots()

ax_V.scatter(x, V_All_data[0], c='red')
ax_P.scatter(x, P_All_data[0], c='blue')



import numpy as np
def update_V(frame):

    ax_V.clear()
    ax_V.scatter(x, V_All_data[frame], c='red')

    ax_V.set_title("Скорость")
    ax_V.set_xlabel('x')
    ax_V.set_ylabel('V')
    ax_V.grid()
    ax_V.legend()

def update_P(frame):

    ax_P.clear()
    ax_P.scatter(x, P_All_data[frame], c='blue', label='Шаг %.0f'%(frame))

    ax_P.set_title("Давление")
    ax_P.set_xlabel('x')
    ax_P.set_ylabel('P')
    ax_P.grid()
    ax_P.legend()



ani_V = animation.FuncAnimation(fig=fig_V, func=update_V, frames=steps, interval=50)
ani_P = animation.FuncAnimation(fig=fig_P, func=update_P, frames=steps, interval=50)

plt.show()
#ani_V.save('gif/1D/Acoustic_1D_discontinious_V_new.gif', writer='pillow')
#ani_P.save('gif/1D/Acoustic_2D_discontinious_P_new.gif', writer='pillow')