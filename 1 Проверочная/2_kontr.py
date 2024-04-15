from matplotlib import pyplot as plt
from math import pi, sin

def W(v, p):
    w1 = -c*rho/2*v + 1/2*p
    w2 = c*rho/2*v + 1/2*p
    return w1, w2


def W_1(w1, w2):
    v = -1/c/rho*w1 + 1/c/rho*w2
    p = w1 + w2
    return v, p

def Border_left(tay):               # Задаёт давление/скорость на левой стенке
    if(tay< 0.1):
        return sin(20*pi*tay)
    else:
        return 0

def Border_right(tay):               # Задаёт давление/скорость на правой стенке
    if (tay < 0.1):
        return sin(20 * pi * tay)
    else:
        return 0


# w1_i_next - значение на n шаге в i+1 точке
# w1_i_prev - значение на n шаге в i-1 точке
# w1_next_i - значение на n+1 шаге в i-ой точке
def interp(P_or_V_left, P_or_V_right, t):
    temp = 1
    for j in range(int(t/dt)):
        for i in range(N):
            w1_i, w2_i = W(curr_v[i], curr_p[i])            # Считаем честно из матрицы
            if(P_or_V_left == "P"):                         # Смотрим какое граничное условие на левой стенке
                P_Border_left = Border_left(j*dt)
            if(P_or_V_right == "P"):                        # Смотрим какое граничное условие на правой стенке
                P_Border_right = Border_right(j*dt)

            if(i == 0):
                w1_i_next, w2_i_next = W(curr_v[i+1], curr_p[i+1])
                w1_next_i = w1_i - c * dt / dx * (w1_i - w1_i_next)   # Честно посчитали т.к. справа знаем данные
                w2_next_i = P_Border_left - w1_next_i                 # Из граничного условия P = w1 + w2 на каждом шаге

            if(i > 0 and i < N-1):
                w1_i_prev, w2_i_prev = W(curr_v[i-1], curr_p[i-1])
                w1_i_next, w2_i_next = W(curr_v[i+1], curr_p[i+1])
                w1_next_i = w1_i - c * dt / dx * (w1_i - w1_i_next)
                w2_next_i = w2_i - c * dt / dx * (w2_i - w2_i_prev)

            if (i == N - 1):
                w1_i_prev, w2_i_prev = W(curr_v[i-1], curr_p[i-1])
                w2_next_i = w2_i - c * dt / dx * (w2_i - w2_i_prev)
                w1_next_i = P_Border_right - w2_next_i

            vi_next, pi_next = W_1(w1_next_i, w2_next_i)
            next_v[i] = vi_next
            next_p[i] = pi_next
        for i in range(N):
            curr_v[i] = next_v[i]
            curr_p[i] = next_p[i]
        tay = 0.3                          # Отображаем график перед столкновением волн в этот момент времени
        if(j*dt >= tay and temp==1):       # Ориентировнчно здесь
            temp = 0 # Чтобы отобразить 1 раз
            for i in range(1, N):
                plt.scatter(i / N * L, curr_p[i], c='#f12813')
            plt.scatter(0, curr_p[0], c='#f12813', label='Давление прямо перед столкновением, t = %.3f' %(tay))


    for i in range(1, N):
        plt.scatter(i / N * L, curr_p[i], c='#343ff4')
    plt.scatter(0, curr_p[0], c='#343ff4', label='Давление в момент при столкновении, t = %.3f' %(t))
    print("Интерполяция завершена")


N = 400
rho = 1000
c = 1500
L = 1000
t = 0.407
dx = L/N
Ku = 0.9
dt = Ku*dx/c
curr_v = [0] * N
curr_p = [0] * N
next_v = [0] * N
next_p = [0] * N




# Задание начальных условий
for i in range(N):
    curr_v[i] = 0
    curr_p[i] = 0


interp("P", "P", t)    # Преобразование



plt.minorticks_on()
plt.grid(which='major',
         color='grey',
         linewidth=1)
plt.grid(which='minor',
         color='k',
         linewidth=0.3,
         linestyle=":")
plt.legend()
plt.xlabel("x, м")
plt.ylabel("P")
plt.title("Интерференция двух синусов")



plt.show()