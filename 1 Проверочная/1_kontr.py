from matplotlib import pyplot as plt
from math import sin, pi
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




L = 1               # Длина "Трубы"
h = [L/100, L/200, L/400, L/600, L/800, L/1000]
Ku = 100/101        # Число Куранта. Нужно брать 100/целое чтобы T/dt было целым числом (иначе нужное кол-во шагов не получится сделать). (100/101, 100/200, 100/158,...)
v = 1               # Скорость волны (вообще ни на что не влияет)
dt = []             # Массив шагов по времени чтобы число куранта было константой
T = 1 * L / v       # 1 полный оборот начального возмущения

for i in range(len(h)):
    dt.append(Ku*h[i]/v)
    print(T / dt[i], round(T/dt[i]))

print("Если отличия от двух чисел невелико, то волна пройдёт полный период", '\n')

N = []                      # Массив кол-ва точек на участке длины L
for i in range(len(h)):
    N.append(int(L/h[i]))
print(N, '\n')


data = [[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []]     # Начальные данные
for k in range(len(h)):
    for i in range(N[k]):
        data[k].append(sin(i/N[k]*pi)**4)


# Массивы наибольших погрешностей для каждого полинома
E1 = []
E2 = []
E3 = []



# Получаем данные 1ой, 2ой и 3ей степени
# curr_i - данные i-го полинома в момент времени t
# next_i - данные i-го полинома в момент времени t + dt
for k in range(len(h)):
    curr_1, curr_2, curr_3 = [], [], []
    for i in range(N[k]):
        curr_1.append(data[k][i])
        curr_2.append(data[k][i])
        curr_3.append(data[k][i])
    next_1 = [0] * N[k]
    next_2 = [0] * N[k]
    next_3 = [0] * N[k]

    for i in range(round(T / dt[k]) ):
        for j in range(N[k]):
            next_1[j] = curr_1[j] - v * dt[k] / h[k] * (curr_1[j] - curr_1[j - 1])
            next_2[j] = ((v * dt[k] / h[k]) ** 2) / 2 * (curr_2[j - 1] + curr_2[(j + 1) % N[k]] - 2 * curr_2[j]) - (
                        v * dt[k] / h[k]) / 2 * (curr_2[(j + 1) % N[k]] - curr_2[j - 1]) + curr_2[j]

            # Третья степень
            a_3 = v*dt[k]/h[k]/6
            b_3 = (v*dt[k]/h[k])**2/2
            c_3 = (v*dt[k]/h[k])**3/6
            next_3[j] = curr_3[j] - a_3 * (2 * curr_3[(j + 1) % N[k]] + 3 * curr_3[j] - 6 * curr_3[j - 1] + curr_3[j - 2])\
                        + b_3 * (curr_3[(j + 1) % N[k]] - 2 * curr_3[j] + curr_3[j - 1])\
                        - c_3 * (curr_3[(j + 1) % N[k]] + 3 * curr_3[j - 1] - 3 * curr_3[j] - curr_3[j - 2])


        for j in range(N[k]):
            curr_1[j] = next_1[j]
            curr_2[j] = next_2[j]
            curr_3[j] = next_3[j]


    temp1, temp2, temp3 = [], [], []
    for i in range(N[k]):
        temp1.append(abs(curr_1[i] - data[k][i]))
        temp2.append(abs(curr_2[i] - data[k][i]))
        temp3.append(abs(curr_3[i] - data[k][i]))

    E1.append(max(temp1))
    E2.append(max(temp2))
    E3.append(max(temp3))


    # Отображение k-го графика

    if(k == 3):
        plt.figure(k)
        for j in range(1,N[k]):
            plt.scatter(j, data[k][j], c='r')
            plt.scatter(j, curr_1[j], c='b')
            plt.scatter(j, curr_2[j], c='green')
            plt.scatter(j, curr_3[j], c='red')

        plt.scatter(0, data[k][0], c='black', label='Исходные данные')
        plt.scatter(0, curr_1[0], c='b', label='Полином 1-ой степени')
        plt.scatter(0, curr_2[0], c='green', label='Полином 2-ой степени')
        plt.scatter(0, curr_3[0], c='red', label='Полином 3-ой степени')

        plt.minorticks_on()
        plt.grid(which='major',
                 color='grey',
                 linewidth=1)
        plt.grid(which='minor',
                 color='k',
                 linewidth=0.3,
                 linestyle=":")
        plt.legend()

    #plt.title("dx/(L/M) %.2f"%(dx[k]*M/L))


log_dx=[]
log_E1=[]
log_E2=[]
log_E3=[]

for i in range(len(h)):
    log_dx.append(log(1/h[i]))
    log_E1.append(abs(log(E1[i])))
    log_E2.append(abs(log(E2[i])))
    log_E3.append(abs(log(E3[i])))


plt.figure(len(h))
k1, sig1, b1 = MNK(log_dx, log_E1, len(log_E1))
k2, sig2, b2 = MNK(log_dx, log_E2, len(log_E2))
k3, sig3, b3 = MNK(log_dx, log_E3, len(log_E3))


plt.plot([min(log_dx), max(log_dx)], [k1 * min(log_dx) + b1, k1 * max(log_dx) + b1], c='b')
plt.plot([min(log_dx), max(log_dx)], [k2 * min(log_dx) + b2, k2 * max(log_dx) + b2], c='g')
plt.plot([min(log_dx), max(log_dx)], [k3 * min(log_dx) + b3, k3 * max(log_dx) + b3], c='r')


for i in range(len(h)):
    plt.scatter(log_dx[i], log_E1[i], c='b')          # Первый полином синий
    plt.scatter(log_dx[i], log_E2[i], c='g')          # Второй полином зелёный
    plt.scatter(log_dx[i], log_E3[i], c='r')          # Третий полином красный

plt.xlabel("$ln(1/h)$")
plt.ylabel("|ln(E)|")


plt.minorticks_on()
plt.grid(which='major',
         color='grey',
         linewidth=1)
plt.grid(which='minor',
         color='k',
         linewidth=0.3,
         linestyle=":")



plt.show()