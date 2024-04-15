from matplotlib import pyplot as plt
import time

# 1. Задание функций (матриц и граничных условий)
def W(v, p):
    w1 = -c*rho/2*v + 1/2*p
    w2 = c*rho/2*v + 1/2*p
    return w1, w2

def W_1(w1, w2):
    v = -1/c/rho*w1 + 1/c/rho*w2
    p = w1 + w2
    return v, p

def Border_left():
    return int(0)

def Border_right():
    return int(0)



# 2. Задание начальных параметров

N = 150
rho = 1
c = 3
L = 2
T = 0.6 * L / c     # 1 полный оборот начального возмущения
hx = L/N
Ku = 0.8
dt = Ku*hx/c




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

for i in range(50):
    data0_v[50+i] = 1
    data0_p[50+i] = 1   #*rho*c

    curr_v[50 + i] = data0_v[50 + i]
    curr_p[50 + i] = data0_p[50 + i]



# 4. Сам код


# w1_next - значение на n шаге в i+1 точке
# w2_prev - значение на n шаге в i-1 точке
# w1_new - значение на n+1 шаге в i-ой точке
# next => R, prev => L. (w1 справа, w2 слева)

plt.ion()
fig, ax = plt.subplots()
ax.set(xlim=[0, L], ylim=[-2.5*max(max(curr_p), max(curr_v)), 2.5*max(max(curr_p), max(curr_v))], xlabel='L', ylabel='P, V')
line1, = ax.plot(x, curr_p, c = 'blue', label = 'Давление')
line2, = ax.plot(x, curr_v, c = 'red', label = 'Скорость')

ax.legend()
ax.grid()

for j in range(int(T/dt)):
    for i in range(N):
        w1_i, w2_i = W(curr_v[i], curr_p[i])
        P_Border_left = Border_left()
        P_Border_right = Border_right()

        if(i == 0):
            w1_next, _ = W(curr_v[i+1], curr_p[i+1])
            w1_new = w1_i - c * dt / hx * (w1_i - w1_next)   # Честно посчитали т.к. справа знаем данные
            w2_new = P_Border_left - w1_new                 # Из граничного условия P = w1 + w2 на каждом шаге

        if(i > 0 and i < N-1):
            _, w2_prev = W(curr_v[i-1], curr_p[i-1])
            w1_next, _ = W(curr_v[i+1], curr_p[i+1])
            w1_new = w1_i - c * dt / hx * (w1_i - w1_next)
            w2_new = w2_i - c * dt / hx * (w2_i - w2_prev)

        if (i == N - 1):
            _, w2_prev = W(curr_v[i-1], curr_p[i-1])
            w2_new = w2_i - c * dt / hx * (w2_i - w2_prev)
            w1_new = P_Border_right - w2_new

        vi_next, pi_next = W_1(w1_new, w2_new)
        next_v[i] = vi_next
        next_p[i] = pi_next
    for i in range(N):
        curr_v[i] = next_v[i]
        curr_p[i] = next_p[i]

    line1.set_ydata(curr_p)
    line2.set_ydata(curr_v)
    plt.draw()
    plt.gcf().canvas.flush_events()

    time.sleep(0.0001)

print("Интерполяция завершена")
plt.ioff()
plt.show()



# 5. График начальных и конечных данных

plt.figure(1)
for i in range(1, N):
    plt.scatter(i, data0_v[i], c='#6b1006')
    plt.scatter(i, data0_p[i], c='#070e7e')
plt.scatter(0, data0_v[0], c='#6b1006', label='Начальная скорость')
plt.scatter(0, data0_p[0], c='#070e7e', label='Начальное давление')

for i in range(1, N):
    plt.scatter(i, curr_v[i], c='#f12813')
    plt.scatter(i, curr_p[i], c='#343ff4')
plt.scatter(0, curr_v[0], c='#f12813', label='Конечная скорость')
plt.scatter(0, curr_p[0], c='#343ff4', label='Конечное давление')


plt.minorticks_on()
plt.grid(which='major',
         color='grey',
         linewidth=1)
plt.grid(which='minor',
         color='k',
         linewidth=0.3,
         linestyle=":")
plt.legend()
plt.show()