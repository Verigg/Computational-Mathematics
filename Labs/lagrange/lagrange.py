import math
import numpy
from sympy import symbols, diff, log
from tabulate import tabulate

import csv
import json
import math
import os

## Вариант 21: y = x^2-lg(0.5x), [a,b] = 0.5, 1.0
## x* = 0.92, x** = 0.53, x*** = 0.98, x**** = 0.77


def f(x):
     return x**2 - math.log10(x/2)


def L(n, x, a, b):
    linspace = numpy.linspace(a, b, n)
    sum = 0
    for i in range(n):
        sum += f(linspace[i]) * Phi(n, linspace, i, x)
    return sum


def Phi(n, linspace, i, x):
    res = 1
    for j in range(n):
        if j != i:
            res *= (x - linspace[j]) / (linspace[i] - linspace[j])
    return res


def abs_err(f, l_n): #Абсолютная ошибка
    return abs(f - l_n)


def rel_err(f, l_n): #Относительная ошибка
    return abs_err(f, l_n) / abs(f) * 100.0


def R(n, x, a, b):
    #Вычисление производной
    x_ = symbols('x')
    f_ = x_ ** 2 - log(x_ / 2) / log(10)
    deriv_n = diff(f_, x_, n+1)

    #Вычисление остаточного члена
    r = abs(float(deriv_n.subs(x_, x))) / math.factorial(n+1) * (b-a)**(n+1)

    return r


a = 0.5
b = 1
x_values = [0.92, 0.53, 0.98, 0.77]
n_values = [3, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

dirout_path = f"data/output"
os.makedirs(dirout_path, exist_ok=True)

for x in x_values:
    print("x = ",str(x), ":")
    table = []

    fout_path = f"{dirout_path}/x_{x}.csv"
    with open(fout_path, 'w') as fout:
        writer = csv.writer(fout)
        writer.writerow(["n", "Abs_Err", "Rel_Err", "R_n"])

        for n in n_values:
            f_res = f(x)
            l_n = L(n, x, a, b)
            abs_err_result = abs_err(f_res, l_n)
            rel_err_result = rel_err(f_res, l_n)
            r_n_result = R(n, x, a, b)

            line = [n, abs_err_result, rel_err_result, r_n_result]
            writer.writerow(line)
            table.append(line)

    print(tabulate(table, headers=["N", "Abs_Err", "Rel_Err", "R_n"], tablefmt="outline",), "\n")



