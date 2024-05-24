from matplotlib import pyplot as plt
import matplotlib.animation as animation
from math import sin, cos, pi, floor, sqrt, log


# 1. Задание функций (матриц и граничных условий)
def Wx(Vx, Vy, Gxx, Gxy, Gyy):
    w1 = -Gxy*lamda/mu + Gyy
    w2 = Gxy*lamda/(2*mu) + Vy*lamda/(2*sqrt(mu/rho))
    w3 = Gxy*lamda/(2*mu) - Vy*lamda/(2*sqrt(mu/rho))
    w4 = Gxx/2 + Vx*(lamda + 2*mu)/(2*sqrt((lamda + 2*mu)/rho))
    w5 = Gxx/2 - Vx*(lamda + 2*mu)/(2*sqrt((lamda + 2*mu)/rho))

    return w1, w2, w3, w4, w5

def W_1x(w1, w2, w3, w4, w5):
    Vx = w4 * sqrt((lamda + 2 * mu) / rho) / (lamda + 2 * mu) - w5 * sqrt((lamda + 2 * mu) / rho) / (lamda + 2 * mu)
    Vy = w2 * sqrt(mu / rho) / lamda - w3 * sqrt(mu / rho) / lamda
    Gxx = w4 + w5
    Gxy = mu * w2 / lamda + mu * w3 / lamda
    Gyy = w1 + w2 + w3

    return Vx, Vy, Gxx, Gxy, Gyy


def Wy(Vx, Vy, Gxx, Gxy, Gyy):
    w1 = -Gxy * lamda / mu + Gxx
    w2 = Gxy * lamda / (2 * mu) + Vx * lamda / (2 * sqrt(mu / rho))
    w3 = Gxy * lamda / (2 * mu) - Vx * lamda / (2 * sqrt(mu / rho))
    w4 = Gyy / 2 + Vy * (lamda + 2 * mu) / (2 * sqrt((lamda + 2 * mu) / rho))
    w5 = Gyy / 2 - Vy * (lamda + 2 * mu) / (2 * sqrt((lamda + 2 * mu) / rho))

    return w1, w2, w3, w4, w5

def W_1y(w1, w2, w3, w4, w5):
    Vx = w2 * sqrt(mu / rho) / lamda - w3 * sqrt(mu / rho) / lamda
    Vy = w4 * sqrt((lamda + 2 * mu) / rho) / (lamda + 2 * mu) - w5 * sqrt((lamda + 2 * mu) / rho) / (lamda + 2 * mu)
    Gxx = w1 + w2 + w3
    Gxy = mu * w2 / lamda + mu * w3 / lamda
    Gyy = w4 + w5

    return Vx, Vy, Gxx, Gxy, Gyy


# Задаёт давление/скорость на левой стенке
def Border_x_left(w2, w4, tay, y):
    Gxx_left, Gxy_left = 1000000*sin(tay/dt * pi/20), 0     #tay/dt = k

    w3 = lamda/mu * Gxy_left - w2
    w5 = Gxx_left - w4
    return w3, w5



def Border_x_right(w3, w5, tay, y):
    Gxx_right, Gxy_right = 1000000*sin(tay/dt * pi/20), 0

    w2 = lamda/mu * Gxy_right - w3
    w4 = Gxx_right - w5
    return w2, w4



def Border_y_down(w2, w4, tay, x):
    Gyy_down, Gxy_down = 1000000*sin(tay/dt * pi/20), 0

    w3 =  lamda/mu * Gxy_down - w2
    w5 = Gyy_down - w4
    return w3, w5


def Border_y_up(w3, w5, tay, x):
    Gyy_up, Gxy_up = 1000000*sin(tay/dt * pi/20), 0

    w2 =  lamda/mu * Gxy_up - w3
    w4 = Gyy_up - w5
    return w2, w4


# 2. Задание параметров

Nx, Ny = 40, 40
rho = 7800
E = 200*10**9
nu = 0.3
lamda = nu*E/((1+nu)*(1-2*nu))
mu = E/(2*(1+nu))
cp = sqrt((lamda + 2*mu)/rho)   # Скорость продольной волны
cs = sqrt(mu/rho)               # Скорость поперечной волны
L = 1
Ku = 1
hx = L/Nx
hy = L/Ny
dt = Ku * min(hx, hy)/max(cp, cs)
steps = 300
print(dt)
print("Сталь:")
print("E =", E)
print("ν =", nu)
print("λ = ", int(lamda))
print("μ = ", int(mu))
print("ρ = ", rho)
print("cp = ", int(cp))
print("cs = ", int(cs))
print()



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


# Начальные данные
data0_Vx = []
data0_Vy = []
data0_Gxx = []
data0_Gxy = []
data0_Gyy = []


# Данные на n-ом шаге
curr_Vx = []
curr_Vy = []
curr_Gxx = []
curr_Gxy = []
curr_Gyy = []


# Данные на n+1-ом шаге для движения по оси х
next_x_Vx = [0] * Nx
next_x_Vy = [0] * Nx
next_x_Gxx = [0] * Nx
next_x_Gxy = [0] * Nx
next_x_Gyy = [0] * Nx

# Данные на n+1-ом шаге для движения по оси y
next_y_Vx = [0] * Nx
next_y_Vy = [0] * Nx
next_y_Gxx = [0] * Nx
next_y_Gxy = [0] * Nx
next_y_Gyy = [0] * Nx


# Задаём начальные условия
for i in range(Ny):
    data0_Vx.append([])
    data0_Vy.append([])
    data0_Gxx.append([])
    data0_Gxy.append([])
    data0_Gyy.append([])

    curr_Vx.append([])
    curr_Vy.append([])
    curr_Gxx.append([])
    curr_Gxy.append([])
    curr_Gyy.append([])


    for j in range(Nx):
        Vx0, Vy0, Gxx0, Gxy0, Gyy0 = 0, 0, 0, 0, 0              # Значение в начальный момент времени (нули)
        if (0.4*Nx < j < 0.6*Nx and 0.4*Ny < i < 0.6*Ny):       # Тут можно задать не нули
            Vx0, Vy0, Gxx0, Gxy0, Gyy0 = 0, 0, 0, 0, 0

        data0_Vx[i].append(Vx0)
        data0_Vy[i].append(Vy0)
        data0_Gxx[i].append(Gxx0)
        data0_Gxy[i].append(Gxy0)
        data0_Gyy[i].append(Gyy0)

        curr_Vx[i].append(data0_Vx[i][j])
        curr_Vy[i].append(data0_Vy[i][j])
        curr_Gxx[i].append(data0_Gxx[i][j])
        curr_Gxy[i].append(data0_Gxy[i][j])
        curr_Gyy[i].append(data0_Gyy[i][j])


#Весь массив данных
V_All_data = []
Gxx_All_data = []
Gxy_All_data = []
Gyy_All_data = []


for k in range(steps+1):
    V_All_data.append([])
    Gxx_All_data.append([])
    Gxy_All_data.append([])
    Gyy_All_data.append([])

    for i in range(Ny):
        V_All_data[k].append([])
        Gxx_All_data[k].append([])
        Gxy_All_data[k].append([])
        Gyy_All_data[k].append([])

        for j in range(Nx):
            V_All_data[k][i].append(0)
            Gxx_All_data[k][i].append(0)
            Gxy_All_data[k][i].append(0)
            Gyy_All_data[k][i].append(0)


for i in range(Ny):
    for j in range(Nx):
        V_All_data[0][i][j] = sqrt(data0_Vx[i][j]**2 + data0_Vy[i][j]**2)
Gxx_All_data[0] = data0_Gxx
Gxy_All_data[0] = data0_Gxy
Gyy_All_data[0] = data0_Gyy






# 4. Сам код

# w1_i - значение на n шаге в i точке
# w1_next - значение на n шаге в i+1 точке
# w1_prev - значение на n шаге в i-1 точке
# w1_new - значение на n+1 шаге в i-ой точке


for k in range(1, steps+1):

    # Считаем по оси х при определенном y
    for i in range(Ny):
        for j in range(Nx):
            w1_i, w2_i, w3_i, w4_i, w5_i = Wx(curr_Vx[i][j], curr_Vy[i][j], curr_Gxx[i][j], curr_Gxy[i][j], curr_Gyy[i][j])

            if(j == 0):
                _, w2_next, _, w4_next, _ = Wx(curr_Vx[i][j+1], curr_Vy[i][j+1], curr_Gxx[i][j+1], curr_Gxy[i][j+1], curr_Gyy[i][j+1])
                w2_new = w2_i - cs * dt / hx * (w2_i - w2_next)
                w4_new = w4_i - cp * dt / hx * (w4_i - w4_next)

                w3_new, w5_new = Border_x_left(w2_new, w4_new, k * dt, i)

            if(j > 0 and j < Nx-1):
                _, w2_next, _, w4_next, _ = Wx(curr_Vx[i][j+1], curr_Vy[i][j+1], curr_Gxx[i][j+1], curr_Gxy[i][j+1], curr_Gyy[i][j+1])
                _, _, w3_prev, _, w5_prev = Wx(curr_Vx[i][j-1], curr_Vy[i][j-1], curr_Gxx[i][j-1], curr_Gxy[i][j-1], curr_Gyy[i][j-1])
                w2_new = w2_i - cs * dt / hx * (w2_i - w2_next)
                w4_new = w4_i - cp * dt / hx * (w4_i - w4_next)
                w3_new = w3_i - cs * dt / hx * (w3_i - w3_prev)
                w5_new = w5_i - cp * dt / hx * (w5_i - w5_prev)


            if (j == Nx - 1):
                _, _, w3_prev, _, w5_prev = Wx(curr_Vx[i][j-1], curr_Vy[i][j-1], curr_Gxx[i][j-1], curr_Gxy[i][j-1], curr_Gyy[i][j-1])
                w3_new = w3_i - cs * dt / hx * (w3_i - w3_prev)
                w5_new = w5_i - cp * dt / hx * (w5_i - w5_prev)

                w2_new, w4_new = Border_x_right(w3_new, w5_new, k * dt, i)

            next_x_Vx[j], next_x_Vy[j], next_x_Gxx[j], next_x_Gxy[j], next_x_Gyy[j] = W_1x(w1_i, w2_new, w3_new, w4_new, w5_new)

        for j in range(Nx):
            curr_Vx[i][j] = next_x_Vx[j]
            curr_Vy[i][j] = next_x_Vy[j]
            curr_Gxx[i][j] = next_x_Gxx[j]
            curr_Gxy[i][j] = next_x_Gxy[j]
            curr_Gyy[i][j] = next_x_Gyy[j]




    # Считаем по оси y при определенном x
    for j in range(Nx):
        for i in range(Ny):
            w1_i, w2_i, w3_i, w4_i, w5_i = Wy(curr_Vx[i][j], curr_Vy[i][j], curr_Gxx[i][j], curr_Gxy[i][j], curr_Gyy[i][j])
            if (i == 0):
                _, w2_next, _, w4_next, _ = Wy(curr_Vx[i + 1][j], curr_Vy[i+1][j], curr_Gxx[i + 1][j], curr_Gxy[i + 1][j], curr_Gyy[i + 1][j])
                w2_new = w2_i - cs * dt / hy * (w2_i - w2_next)  # Честно посчитали т.к. справа знаем данные
                w4_new = w4_i - cp * dt / hy * (w4_i - w4_next)

                w3_new, w5_new = Border_y_down(w2_new, w4_new, k * dt, j)

            if (i > 0 and i < Ny - 1):
                _, w2_next, _, w4_next, _ = Wy(curr_Vx[i + 1][j], curr_Vy[i+1][j], curr_Gxx[i + 1][j], curr_Gxy[i + 1][j], curr_Gyy[i + 1][j])
                _, _, w3_prev, _, w5_prev = Wy(curr_Vx[i - 1][j], curr_Vy[i - 1][j], curr_Gxx[i - 1][j], curr_Gxy[i - 1][j], curr_Gyy[i - 1][j])
                w2_new = w2_i - cs * dt / hy * (w2_i - w2_next)
                w4_new = w4_i - cp * dt / hy * (w4_i - w4_next)
                w3_new = w3_i - cs * dt / hy * (w3_i - w3_prev)
                w5_new = w5_i - cp * dt / hy * (w5_i - w5_prev)

            if (i == Ny - 1):
                _, _, w3_prev, _, w5_prev = Wy(curr_Vx[i - 1][j], curr_Vy[i - 1][j], curr_Gxx[i - 1][j], curr_Gxy[i - 1][j], curr_Gyy[i - 1][j])
                w3_new = w3_i - cs * dt / hy * (w3_i - w3_prev)
                w5_new = w5_i - cp * dt / hy * (w5_i - w5_prev)

                w2_new, w4_new = Border_y_up(w3_new, w5_new, k * dt, j)

            next_y_Vx[i], next_y_Vy[i], next_y_Gxx[i], next_y_Gxy[i], next_y_Gyy[i] = W_1y(w1_i, w2_new, w3_new, w4_new, w5_new)

        for i in range(Ny):
            curr_Vx[i][j] = next_y_Vx[i]
            curr_Vy[i][j] = next_y_Vy[i]
            curr_Gxx[i][j] = next_y_Gxx[i]
            curr_Gxy[i][j] = next_y_Gxy[i]
            curr_Gyy[i][j] = next_y_Gyy[i]


            V_All_data[k][i][j] = sqrt(curr_Vx[i][j]**2 + curr_Vy[i][j]**2)
            Gxx_All_data[k][i][j] = curr_Gxx[i][j]
            Gxy_All_data[k][i][j] = curr_Gxy[i][j]
            Gyy_All_data[k][i][j] = curr_Gyy[i][j]



print("Интерполяция завершена")



# 5. Отображение анимации

fig_V, ax_V = plt.subplots()
fig_Gxx, ax_Gxx = plt.subplots()
fig_Gxy, ax_Gxy = plt.subplots()
fig_Gyy, ax_Gyy = plt.subplots()


ax_V.scatter(x, y, c=V_All_data[0], cmap='autumn')
ax_Gxx.scatter(x, y, c=Gxx_All_data[0], cmap='winter')
ax_Gxy.scatter(x, y, c=Gxy_All_data[0], cmap='winter')
ax_Gyy.scatter(x, y, c=Gyy_All_data[0], cmap='winter')



def update_V(frame):

    ax_V.clear()
    ax_V.scatter(x, y,
               c=V_All_data[frame],
               cmap='autumn', label='Шаг %.0f'%(frame))

    ax_V.set_title("Скорость")
    ax_V.set_xlabel('x')
    ax_V.set_ylabel('y')
    ax_V.legend()

def update_Gxx(frame):

    ax_Gxx.clear()
    ax_Gxx.scatter(x, y,
               c=Gxx_All_data[frame],
               cmap='winter', label='Шаг %.0f'%(frame))

    ax_Gxx.set_title("Gxx")
    ax_Gxx.set_xlabel('x')
    ax_Gxx.set_ylabel('y')
    ax_Gxx.legend()

def update_Gxy(frame):

    ax_Gxy.clear()
    ax_Gxy.scatter(x, y,
               c=Gxy_All_data[frame],
               cmap='winter', label='Шаг %.0f'%(frame))

    ax_Gxy.set_title("Gxy")
    ax_Gxy.set_xlabel('x')
    ax_Gxy.set_ylabel('y')
    ax_Gxy.legend()

def update_Gyy(frame):

    ax_Gyy.clear()
    ax_Gyy.scatter(x, y,
               c=Gyy_All_data[frame],
               cmap='winter', label='Шаг %.0f'%(frame))

    ax_Gyy.set_title("Gyy")
    ax_Gyy.set_xlabel('x')
    ax_Gyy.set_ylabel('y')
    ax_Gyy.legend()


interval = 100
ani_V = animation.FuncAnimation(fig=fig_V, func=update_V, frames=steps, interval=interval)
ani_Gxx = animation.FuncAnimation(fig=fig_Gxx, func=update_Gxx, frames=steps, interval=interval)
ani_Gxy = animation.FuncAnimation(fig=fig_Gxy, func=update_Gxy, frames=steps, interval=interval)
ani_Gyy = animation.FuncAnimation(fig=fig_Gyy, func=update_Gyy, frames=steps, interval=interval)

plt.show()
#ani_V.save('gif/2D/Acoustic_2D_animation_V_100.gif', writer='pillow')
#ani_P.save('gif/2D/Acoustic_2D_animation_P_100.gif', writer='pillow')
