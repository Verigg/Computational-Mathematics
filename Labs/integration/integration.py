import math
import random

from tabulate import tabulate
from sympy import symbols, diff, log, solveset, Interval, Min, Max, minimum, maximum, nsolve, integrate


## Вариант 21: y = x^2-lg(0.5x), [a,b] = 0.5, 1.0
## x* = 0.92, x** = 0.53, x*** = 0.98, x**** = 0.77
## Лабораторная номер 4

# Часть 1, вариант a - левые прямоугольники

def f(x):
    return x ** 2 - math.log10(x / 2)


def f_integral(a, b):
    x_ = symbols('x')
    f_ = x_ ** 2 - log(x_ / 2) / log(10)
    integral = integrate(f_, x_)  # интеграл от f
    res = integral.subs(x_, b) - integral.subs(x_, a)
    return float(res)


def left_rectangle(n, f_array, h):
    res = 0
    for i in range(n):
        res += f_array[i] * h
    return res


def left_rectangle_err(a, b, h): # Подходит и для right_rectangle_err
    res = (M(a, b, 1) * (b - a)) / 2 * h
    return res


def M(a, b, n):  # максимум производной n-го порядка на интервале a,b
    x_ = symbols('x')
    f_ = x_ ** 2 - log(x_ / 2) / log(10)
    deriv_n = diff(f_, x_, n)  # Производная n-го порядка

    interv = Interval(a, b)  # Максимум производной на интервале
    max = (maximum(deriv_n, x_, interv))
    return float(max)


def right_rectangle(n, f_array, h):
    res = 0
    for i in range(n):
        res += f_array[i+1] * h
    return res


def center_rectangle(n, x_array, h):
    res = 0
    for i in range(n):
        res += f((x_array[i+1] + x_array[i])/2)*h
    return res


def center_rectangle_err(a, b, h):
    res = (M(a, b, 2) * (b - a)) / 24 * h**2
    return res


def trapezoid(n, f_array, h):
    res = 0
    for i in range(n):
        res += (f_array[i] + f_array[i+1]) / 2 * h
    return  res


def trapezoid_err(a, b, h):
    res = (M(a, b, 2) * (b - a)) / 12 * h**2
    return res


def simpson(n, h):
    res = 0
    x_array = []
    for i in range(n + 1):
        x_array.append(a + h * i)

    for i in range(n-1):
        res += h / 6 * (f(x_array[i]) + 4 * f(x_array[i+1]) + f(x_array[i+2]))
        # res += h / 6 * (f(i) + 4 * f(i+1/2) + f(i+1))
    return res


def simpson_err(a, b, h):
    res = (M(a, b, 4) * (b - a)) / 2880 * h**5
    return res


def monte_carlo(n, a, b, h):
    res = 0
    for i in range(n):
        res += f(random.uniform(a,b))
    res *= h
    return res


def monte_carlo_err(a, b, n):
    mu = 0
    for i in range(n):
        mu += (random.uniform(a,b))
    mu /= n

    D = 0
    for i in range(n):
        D += (f(random.uniform(a,b)) - mu)**2
    D /= n * n**0.5
    D *= (b-a)

    return D

a = 0.5
b = 1

# первая часть:

n_array = []

for i in range(10):
    n_array.append(2 ** (i + 1))

table = []
j = 1
I_delta_array = []
I_delta_array.append(0)
for n in n_array:
    h = (b - a) / n
    f_array = []
    x_array = []
    for i in range(n + 1):
        f_array.append(f(a + h * i))
        x_array.append(a + h * i)

    I = left_rectangle(n, f_array, h)
    I_delta = abs(f_integral(a, b) - I)  # абсолютная ошибка
    I_D = I_delta / abs(f_integral(a, b)) * 100

    ratio = I_delta_array[j - 1] / I_delta

    R_n = left_rectangle_err(a, b, h)

    row = [j, n, I, I_delta, I_D, R_n, ratio]
    table.append(row)

    j += 1
    I_delta_array.append(I_delta)

with open(f"inegration_table.txt", "w", encoding="utf-8") as fout:  ## сохраняем в txt
    fout.write(
        "Часть 1 \nЛевые прямоугольники: \n"
    )
    fout.write(
        tabulate(table, tablefmt="grid", headers=["j, №эксп", "n", "I", "ΔI_n", "δI_n", "R_n", "ratio"])
    )

# вторая часть:

n = 10000
h = (b - a) / n
f_array = []
x_array = []
table = []
for i in range(n + 1):
    f_array.append(f(a + h * i))
    x_array.append(a + h * i)

# left_rectangle
I = left_rectangle(n, f_array, h)
I_delta = abs(f_integral(a, b) - I)  # абсолютная ошибка
I_D = I_delta / abs(f_integral(a, b)) * 100
R_n = left_rectangle_err(a, b, h)

row = [left_rectangle.__name__, I, I_delta, I_D, R_n]
table.append(row)

# right_rectangle
I = right_rectangle(n, f_array, h)
I_delta = abs(f_integral(a, b) - I)  # абсолютная ошибка
I_D = I_delta / abs(f_integral(a, b)) * 100
R_n = left_rectangle_err(a, b, h) # == left_rectangle_err

row = [right_rectangle.__name__, I, I_delta, I_D, R_n]
table.append(row)

# center_rectangle
I = center_rectangle(n, x_array, h)
I_delta = abs(f_integral(a, b) - I)  # абсолютная ошибка
I_D = I_delta / abs(f_integral(a, b)) * 100
R_n = center_rectangle_err(a, b, h)

row = [center_rectangle.__name__, I, I_delta, I_D, R_n]
table.append(row)

# trapezoid
I = trapezoid(n, f_array, h)
I_delta = abs(f_integral(a, b) - I)  # абсолютная ошибка
I_D = I_delta / abs(f_integral(a, b)) * 100
R_n = trapezoid_err(a, b, h)

row = [trapezoid.__name__, I, I_delta, I_D, R_n]
table.append(row)

# simson
I = simpson(n, h)
I_delta = abs(f_integral(a, b) - I)  # абсолютная ошибка
I_D = I_delta / abs(f_integral(a, b)) * 100
R_n = simpson_err(a, b, h)

row = [simpson.__name__, I, I_delta, I_D, R_n]
table.append(row)

# monte carlo
I = monte_carlo(n, a, b, h)
I_delta = abs(f_integral(a, b) - I)  # абсолютная ошибка
I_D = I_delta / abs(f_integral(a, b)) * 100
R_n = monte_carlo_err(a,b,n)

row = [monte_carlo.__name__, I, I_delta, I_D, R_n]
table.append(row)

with open(f"inegration_table2.txt", "w", encoding="utf-8") as fout:  ## сохраняем в txt
    fout.write(
        "Часть 2\n"
    )
    fout.write(
        tabulate(table, tablefmt="grid", headers=["Метод", "I", "ΔI_n", "δI_n", "R_n"])
    )