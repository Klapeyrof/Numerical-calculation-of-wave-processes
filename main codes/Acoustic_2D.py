from matplotlib import pyplot as plt
import matplotlib.animation as animation
from math import sin, cos, pi, floor, sqrt, log


# 1. Задание функций (матриц и граничных условий)

def Wx(v, p):
    w1 = -cx*rho/2*v + 1/2*p
    w2 = cx*rho/2*v + 1/2*p
    return w1, w2

def W_1x(w1, w2):
    v = -1/cx/rho*w1 + 1/cx/rho*w2
    p = w1 + w2
    return v, p


def Wy(v, p):
    w1 = -cy*rho/2*v + 1/2*p
    w2 = cy*rho/2*v + 1/2*p
    return w1, w2

def W_1y(w1, w2):
    v = -1/cy/rho*w1 + 1/cy/rho*w2
    p = w1 + w2
    return v, p


# Задаёт функцию на левой границе (aP + bV = f_left(tay,y))
# P = 0 при a = const, b = 0, f_left = 0
def f_left(tay, y):
    #return 0

    # Стоячая волна
    return sin(pi/0.2*(y+1) + 5*pi*tay) + sin(pi/30*(y+1) - 5*pi*tay)


# Задаёт давление/скорость на левой стенке
def Border_x_left(w1, alpha, beta, tay, y):
    return (f_left(tay, y) - w1 * (alpha - beta / (rho * cx))) / ((alpha + beta / (rho * cx)))



def f_right(tay, y):
    return 0 #sin(pi/10*y + pi/5*tay)

def Border_x_right(w1, alpha, beta, tay, y):
    return (f_right(tay, y) - w1 * (alpha - beta / (rho * cx))) / ((alpha + beta / (rho * cx)))



def f_down(tay, x):
    return 0

def Border_y_down(w1, alpha, beta, tay, x):
    return (f_down(tay, x) - w1 * (alpha - beta / (rho * cy))) / ((alpha + beta / (rho * cy)))



def f_up(tay, x):
    return 0

def Border_y_up(w1, alpha, beta, tay, x):
    return (f_left(tay, x) - w1 * (alpha - beta / (rho * cy))) / ((alpha + beta / (rho * cy)))



# 2. Задание параметров

Nx, Ny = 60, 60
rho = 1
cx = 3
cy = 0.001
L = 2
Ku = 1
hx = L/Nx
hy = L/Ny
dt = Ku * min(hx, hy)/max(cx, cy)
steps = 300





# 3. Задание массивов данных и начальных условий

# Сетка
x = []
y = []
for i in range(Ny):
    x.append([])
    y.append([])
    for j in range(Nx):
        x[i].append(j)
        y[i].append(i)

curr_v_x = []               # Данные скоростей по оси х на i-ом шаге
curr_v_y = []               # Данные скоростей по оси y на i-ом шаге
data0_v_x = []              # Начальные данные скоростей по оси х
data0_v_y = []              # Начальные данные скоростей по оси y
curr_p = []                 # Данные давлений на i-ом шаге
data0_p = []                # Начальные данные давлений
next_vx = [0] * Nx          # Данные скоростей на i+1-ом шаге по х
next_px = [0] * Nx          # Данные давлений на i+1-ом шаге по х
next_vy = [0] * Ny          # Данные скоростей на i+1-ом шаге по х
next_py = [0] * Ny          # Данные давлений на i+1-ом шаге по y



# Задаём начальные условия
for i in range(Ny):
    data0_v_x.append([])
    curr_v_x.append([])
    data0_v_y.append([])
    curr_v_y.append([])
    data0_p.append([])
    curr_p.append([])
    for j in range(Nx):
        fvx, fvy, fp = 0, 0, 0              # Значение скоростей и давления в начальный момент времени
        if (0.4*Nx < j < 0.6*Nx and 0.4*Ny < i < 0.6*Ny):
            fvx = 0
            fvy = 0
            fp = 0
        data0_v_x[i].append(fvx)
        data0_v_y[i].append(fvy)

        curr_v_x[i].append(data0_v_x[i][j])
        curr_v_y[i].append(data0_v_y[i][j])

        data0_p[i].append(fp)
        curr_p[i].append(data0_p[i][j])



#Весь массив данных
V_All_data = []               # Весь массив данных
P_All_data = []

for k in range(steps+1):
    V_All_data.append([])
    P_All_data.append([])
    for i in range(Ny):
        V_All_data[k].append([])
        P_All_data[k].append([])
        for j in range(Nx):
            V_All_data[k][i].append(0)
            P_All_data[k][i].append(0)

for i in range(Ny):
    for j in range(Nx):
        V_All_data[0][i][j] = sqrt(data0_v_x[i][j]**2 + data0_v_y[i][j]**2)
P_All_data[0] = data0_p





# 4. Сам код


# w1_next - значение на n шаге в i+1 точке
# w2_prev - значение на n шаге в i-1 точке
# w1_new - значение на n+1 шаге в i-ой точке
# next => R, prev => L. (w1 справа, w2 слева)


for k in range(1, steps+1):

    # Считаем по оси х при определенном y
    for i in range(Ny):
        for j in range(Nx):
            w1_i, w2_i = Wx(curr_v_x[i][j], curr_p[i][j])

            if(j == 0):
                w1_next, _ = Wx(curr_v_x[i][j + 1], curr_p[i][j + 1])
                w1_new = w1_i - cx * dt / hx * (w1_i - w1_next)          # Честно посчитали т.к. справа знаем данные
                w2_new = Border_x_left(w1_new, 1, 0, k * dt, i)             # Из граничного условия aP + bV = f_left(t) на каждом шаге

            if(j > 0 and j < Nx-1):
                _, w2_prev = Wx(curr_v_x[i][j - 1], curr_p[i][j - 1])
                w1_next, _ = Wx(curr_v_x[i][j + 1], curr_p[i][j + 1])
                w1_new = w1_i - cx * dt / hx * (w1_i - w1_next)
                w2_new = w2_i - cx * dt / hx * (w2_i - w2_prev)

            if (j == Nx - 1):
                _, w2_prev = Wx(curr_v_x[i][j - 1], curr_p[i][j - 1])
                w2_new = w2_i - cx * dt / hx * (w2_i - w2_prev)
                w1_new = Border_x_right(w2_new, 1, 0, k * dt, i)  # Из граничного условия aP + bV = f_right(t) на каждом шаге

            vi_next, pi_next = W_1x(w1_new, w2_new)
            next_vx[j] = vi_next
            next_px[j] = pi_next
        for j in range(Nx):
            curr_v_x[i][j] = next_vx[j]
            curr_p[i][j] = next_px[j]

    # Считаем по оси y при определенном x
    for j in range(Nx):
        for i in range(Ny):
            w1_i, w2_i = Wy(curr_v_y[i][j], curr_p[i][j])

            if (i == 0):
                w1_next, _ = Wy(curr_v_y[i + 1][j], curr_p[i + 1][j])
                w1_new = w1_i - cy * dt / hy * (w1_i - w1_next)        # Честно посчитали т.к. справа знаем данные
                w2_new = Border_y_down(w1_new, 1, 0, k * dt, j)             # Из граничного условия aP + bV = f_left(t) на каждом шаге

            if (i > 0 and i < Ny - 1):
                _, w2_prev = Wy(curr_v_y[i - 1][j], curr_p[i - 1][j])
                w1_next, _ = Wy(curr_v_y[i + 1][j], curr_p[i + 1][j])
                w1_new = w1_i - cy * dt / hy * (w1_i - w1_next)
                w2_new = w2_i - cy * dt / hy * (w2_i - w2_prev)

            if (i == Ny - 1):
                _, w2_prev = Wy(curr_v_y[i - 1][j], curr_p[i - 1][j])
                w2_new = w2_i - cy * dt / hy * (w2_i - w2_prev)
                w1_new = Border_y_up(w2_new, 1, 0, k * dt, i)             # Из граничного условия aP + bV = f_left(t) на каждом шаге

            vi_next, pi_next = W_1y(w1_new, w2_new)
            next_vy[i] = vi_next
            next_py[i] = pi_next
        for i in range(Ny):
            curr_v_y[i][j] = next_vy[i]
            curr_p[i][j] = next_py[i]
            V_All_data[k][i][j] = sqrt(curr_v_y[i][j]**2 + curr_v_x[i][j]**2)
            P_All_data[k][i][j] = next_py[i]

print("Интерполяция завершена")





# 5. Отображение анимации

fig_V, ax_V = plt.subplots()
fig_P, ax_P = plt.subplots()

ax_V.scatter(x, y, c=V_All_data[0], cmap='autumn')
ax_P.scatter(x, y, c=P_All_data[0], cmap='winter')


def update_V(frame):

    ax_V.clear()
    ax_V.scatter(x, y,
               c=V_All_data[frame],
               cmap='autumn', label='Шаг %.0f'%(frame))

    ax_V.set_title("Скорость")
    ax_V.set_xlabel('x')
    ax_V.set_ylabel('y')
    ax_V.legend()

def update_P(frame):

    ax_P.clear()
    ax_P.scatter(x, y,
               c=P_All_data[frame],
               cmap='winter', label='Шаг %.0f'%(frame))

    ax_P.set_title("Давление")
    ax_P.set_xlabel('x')
    ax_P.set_ylabel('y')
    ax_P.legend()



ani_V = animation.FuncAnimation(fig=fig_V, func=update_V, frames=steps, interval=50)
ani_P = animation.FuncAnimation(fig=fig_P, func=update_P, frames=steps, interval=50)

plt.show()
#ani_V.save('gif/2D/Acoustic_2D_animation_V_sin_left_2.gif', writer='pillow')
ani_P.save('gif/2D/Acoustic_2D_animation_P_sin_left_4.gif', writer='pillow')
