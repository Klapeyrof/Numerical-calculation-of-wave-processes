from matplotlib import pyplot as plt
import random

curr = [0] * 50 + [1] * 100 + [0] * 50  # Начальные данные

#curr=[]

# Рандомные значения
# for i in range(10):
#     rand_1 = random.randint(10,50)
#     rand_2 = random.randint(0, 1)
#     for j in range(1,rand_1):
#         curr.append(rand_2)

M = len(curr)
curr_1, curr_2 = [], []                 # Данные для полинома 1 и 2 степени
for i in range(M):
    curr_1.append(curr[i])
    curr_2.append(curr[i])
M = len(curr_1)
next_1 = [0] * M
next_2 = [0] * M

dt = 0.01
L = 5
dx = L / M
K = 0.9 # Число Куранта
v = dx/dt*K #
T = L/v*1  # (Полные период)



# Строим график для полинома 1ой и 2ой степени
for i in range(int(T / dt)):
    # for j in range(M):
    #     plt.scatter(j, curr[j], c='b')
    # plt.show()
    next_1[0] = curr_1[len(curr_1) - 1]
    next_2[0] = curr_2[len(curr_2) - 1]
    for j in range(1, M-1):
        next_1[j] = curr_1[j] - v * dt / dx * (curr_1[j] - curr_1[j - 1])
        next_2[j] = ((v*dt/dx)**2)/2*(curr_2[j - 1] + curr_2[j + 1] - 2*curr_2[j]) - (v*dt/dx)/2*(curr_2[j + 1] - curr_2[j - 1]) + curr_2[j]
    next_1[M-1] = curr_1[M-1] - v * dt / dx * (curr_1[M-1] - curr_1[M - 1 - 1])
    next_2[M-1] = ((v * dt / dx) ** 2) / 2 * (curr_2[M-1 - 1] + curr_2[0] - 2 * curr_2[M-1]) - (v * dt / dx) / 2 * (
                curr_2[0] - curr_2[M-1 - 1]) + curr_2[M-1]

    for j in range(M):
        curr_1[j] = next_1[j]
        curr_2[j] = next_2[j]


# Строим график для исходных данных
plt.figure(1)
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
plt.show()


