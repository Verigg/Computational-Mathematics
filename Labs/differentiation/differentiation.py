import math
from tabulate import tabulate
from sympy import symbols, diff, log, solveset, Interval, Min, Max, minimum, maximum, nsolve

## Вариант 21: y = x^2-lg(0.5x), [a,b] = 0.5, 1.0
## x* = 0.92, x** = 0.53, x*** = 0.98, x**** = 0.77

## K = 2, N = 4, M = 1

def f(x):
    return x ** 2 - math.log10(x / 2)

def f_dx(x, n):  # функция находит значение производной n-го порядка в точке x
    x_ = symbols('x')
    f_ = x_ ** 2 - log(x_ / 2) / log(10)
    deriv_n = diff(f_, x_, n)  # Производная n-го порядка
    res = deriv_n.subs(x_, x)
    return float(res)

def Ln_dx(x_array, f_array, m, H):
    x, x0, x1, x2, x3, x4, f0, f1, f2, f3, f4, h = symbols('x, x0, x1, x2, x3, x4, f0, f1, f2, f3, f4, h')
    f = f0/(24*h**4) *  (x - x1) * (x - x2) * (x - x3) * (x - x4)  + \
    f1/(-6 * h**4) * (x - x0) * (x - x2) * (x - x3) * (x - x4) + \
    f2/(4 * h**4) * (x - x0) * (x - x1) * (x - x3) * (x - x4) + \
    f3/(-6 * h**4) * (x - x0) * (x - x1) * (x - x2) * (x - x4) + \
    f4/(24 * h**4) * (x - x0) * (x - x1) * (x - x2) * (x - x3) # *_* wanna die

    deriv_n = diff(f, x, 2)

    res = deriv_n.subs({x0: x_array[0],x1: x_array[1], x2: x_array[2], x3: x_array[3], x4: x_array[4],
                  f0: f_array[0], f1: f_array[1], f2: f_array[2], f3: f_array[3], f4: f_array[4],
                  x: x_array[m], h:H
                  })
    return float(res)



def R_n(n, a, b, x, h):
    x_ = symbols('x')
    f_ = x_ ** 2 - log(x_ / 2) / log(10)
    deriv_n = diff(f_, x_, n + 1)  # Производная n-го порядка

    interv = Interval(a, b)  # Минимум и максимум производной на интервале
    min = (minimum(deriv_n, x_, interv))
    max = (maximum(deriv_n, x_, interv))

    res_min = float(min) / math.factorial(n + 1) * w(n, x, h, a)  # минимальный и максимальный остаточный член
    res_max = float(max) / math.factorial(n + 1) * w(n, x, h, a)

    res = [res_min, res_max]
    res.sort()
    return res


def w(n, x, h, a):
    res = 1
    for i in range(n - 1):
        if x - (a + i * h) != 0:
            res *= x - (a + i * h)
    return res

def print_result(f, l_n, r, x_array):
    print(f"m = {m}, x_m = {x_array[m]}")
    print(f"f^(2)(x_m) = {f:.10f} and L^(2)(x_m) = {l_n:.10f}")
    print(f"Δ = {l_n-f:.5e} between {r[0]:.5e} and {r[1]:.5e}?")
    print("Yes") if r[0] < l_n-f < r[1] else print("No")

k = 2
n = 4
m = 1

a = 0.5
b = 1

h = (b-a)/n

f_array = []
x_array = []

for i in range(n+1):
    f_array.append(f(a+h*i))
    x_array.append(a+h*i)

F = f_dx(x_array[m], k)
Ln = Ln_dx(x_array, f_array, m, h)
Rn = R_n(n, a, b, x_array[m], h)

print_result(F, Ln, Rn, x_array)