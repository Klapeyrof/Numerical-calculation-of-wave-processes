from matplotlib import pyplot as plt
from math import log
def MNK(x, y, n):
    a1, a2, a3, mx, my = 0, 0, 0, 0, 0
    for i in range(n):
        a1 += (x[i] * y[i])
        a2 += ((x[i]) ** 2)
        a3 += ((y[i]) ** 2)
        mx += (x[i])
        my += (y[i])
        if (i == n - 1):
            a1 /= n
            a2 /= n
            a3 /= n
            mx /= n
            my /= n
    k = (a1 - mx * my) / (a2 - (mx ** 2))
    b = my - k * mx
    sig_k = 1 / ((n) ** 0.5) * (((a3 - (my ** 2)) / (a2 - (mx ** 2))) - k ** 2) ** 0.5
    print("k=", '{0:.5}'.format(k), " sig_k=", '{0:.5}'.format(sig_k), " b=", '{0:.5}'.format(b))
    return k, sig_k, b


curr = [0] * 50 + [1] * 100 + [0] * 50  # Массив начальных данных
M = len(curr)
dt = 0.001
L = 5                                   # Длина "трубы"
dx = [L/M*0.92, L/M*0.94, L/M*0.96, L/M*0.98, L/M, L/M*1.02, L/M*1.04, L/M*1.06,L/M*1.08]

K = 0.9 # Число Куранта
# K = 0.95 => коэф-ты 76 и 39   (при shift = 35)
# K = 0.9 => коэф-ты 40 и 20
# K = 0.6 => коэф-ты 5.48 и 1.78
# K = 0.41 => коэф-ты 2.06 и 1.82

shift = 35 # Сдвиг от начала исходных данных, которые в конце будем сравнивать
# Например, при shift = 35 будут сравниваться данные с 50+35 до 200-50-35 индексами
# От этого параметра результат вли


v = dx[0]/dt*K                          # В результате разное число Куранта
T = [L/v*0.92, L/v*0.94, L/v*0.96, L/v*0.98, L/v, L/v*1.02,L/v*1.04,L/v*1.06,L/v*1.08]
for i in range(len(dx)):
    print(v/dx[i]*dt)


E1 = []
E2 = []


# Строим график для полинома 1ой и 2ой степени
for k in range(len(dx)):
    curr_1, curr_2 = [], []
    for i in range(M):
        curr_1.append(curr[i])
        curr_2.append(curr[i])
    next_1 = [0] * M
    next_2 = [0] * M
    #print(T)
    for i in range(int(T[k] / dt)):
        next_1[0] = curr_1[len(curr_1) - 1]
        next_2[0] = curr_2[len(curr_2) - 1]

        for j in range(1, M-1):
            next_1[j] = curr_1[j] - v * dt / dx[k] * (curr_1[j] - curr_1[j - 1])
            next_2[j] = ((v*dt/dx[k])**2)/2*(curr_2[j - 1] + curr_2[j + 1] - 2*curr_2[j]) - (v*dt/dx[k])/2*(curr_2[j + 1] - curr_2[j - 1]) + curr_2[j]
        next_1[M-1] = curr_1[M-1] - v * dt / dx[k] * (curr_1[M-1] - curr_1[M - 1 - 1])
        next_2[M-1] = ((v * dt / dx[k]) ** 2) / 2 * (curr_2[M-1 - 1] + curr_2[0] - 2 * curr_2[M-1]) - (v * dt / dx[k]) / 2 * (
                    curr_2[0] - curr_2[M-1 - 1]) + curr_2[M-1]

        for j in range(M):
            curr_1[j] = next_1[j]
            curr_2[j] = next_2[j]

    temp1,temp2 = [],[]
    for i in range(50+shift,M-50-shift):
        temp1.append(abs(curr_1[i] - curr[i]))
        temp2.append(abs(curr_2[i] - curr[i]))
    E1.append(max(temp1))
    E2.append(max(temp2))

    # Отображение k-го графика

    plt.figure(k)
    for j in range(M):
        plt.scatter(j, curr[j], c='r')
        plt.scatter(j, curr_1[j], c='b')
        plt.scatter(j, curr_2[j], c='green')

    plt.scatter(0, curr[j], c='r', label='Исходные данные')
    plt.scatter(0, curr_1[j], c='b', label='Полином 1-ой степени')
    plt.scatter(0, curr_2[j], c='green', label='Полином 2-ой степени')

    plt.minorticks_on()
    plt.grid(which='major',
             color='grey',
             linewidth=1)
    plt.grid(which='minor',
             color='k',
             linewidth=0.3,
             linestyle=":")
    plt.legend()
    plt.title("dx/(L/M) %.2f"%(dx[k]*M/L))




print(E1)
print(E2)
log_dx=[]
log_E1=[]
log_E2=[]
for i in range(len(dx)):
    log_dx.append(log(1/dx[i]))
    log_E1.append(abs(log(E1[i])))
    log_E2.append(abs(log(E2[i])))



plt.figure(len(dx))
k1, sig1, b1 = MNK(log_dx, log_E1, len(log_E1))
k2, sig2, b2 = MNK(log_dx, log_E2, len(log_E2))

plt.plot([min(log_dx),max(log_dx)], [k1 * min(log_dx) + b1, k1 * max(log_dx) + b1],c='b')
plt.plot([min(log_dx),max(log_dx)], [k2 * min(log_dx) + b2, k2 * max(log_dx) + b2],c='g')


for i in range(len(dx)):
    plt.scatter(log_dx[i], log_E1[i], c='b')          # Первый полином синий
    plt.scatter(log_dx[i], log_E2[i], c='g')          # Второй полином синий
plt.xlabel("$ln(1/dx)$")
plt.ylabel("|ln(E)|")


plt.minorticks_on()
plt.grid(which='major',
         color='grey',
         linewidth=1)
plt.grid(which='minor',
         color='k',
         linewidth=0.3,
         linestyle=":")







# График в исходных координатах
plt.figure(len(dx)+1)
for i in range(len(dx)):
    plt.scatter(dx[i], E1[i], c='b')
    plt.scatter(dx[i], E2[i], c='g')
plt.xlabel("$dx$")
plt.ylabel("E")


plt.minorticks_on()
plt.grid(which='major',
         color='grey',
         linewidth=1)
plt.grid(which='minor',
         color='k',
         linewidth=0.3,
         linestyle=":")

plt.show()