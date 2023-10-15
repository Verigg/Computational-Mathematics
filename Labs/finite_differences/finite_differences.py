import math
from tabulate import tabulate
from sympy import symbols, diff, log, solveset, Interval, Min, Max, minimum, maximum, nsolve


## Вариант 21: y = x^2-lg(0.5x), [a,b] = 0.5, 1.0
## x* = 0.92, x** = 0.53, x*** = 0.98, x**** = 0.77


def f(x):
    return x ** 2 - math.log10(x / 2)


def finite_diffs(n, a, h):  # создаём таблицу конечных разностей
    table = []

    for i in range(n + 1):
        row = []
        for j in range(n + 1):
            row.append(0.0)
        table.append(row)

    for i in range(n + 1):
        table[i][0] = f(a + i * h)

    for j in range(1, n + 1):
        for i in range(n - j + 1):
            table[i][j] = table[i + 1][j - 1] - table[i][j - 1]

    with open(f"table.txt", "w", encoding="utf-8") as fout:  ## сохраняем в txt
        fout.write(
            tabulate(table, tablefmt="grid", )
        )

    return table


def newton_first(table, t):
    n = len(table)
    l_n_x = table[0][0]
    coeff = t
    for i in range(n - 1):
        l_n_x = l_n_x + (coeff / math.factorial(i + 1)) * table[0][i + 1]
        coeff = coeff * (t - (i + 1))
    return l_n_x


def newton_second(table, t):
    n = len(table)
    l_n_x = table[n - 1][0]
    coeff = t
    for i in range(n - 1):
        l_n_x = l_n_x + (coeff / math.factorial(i + 1)) * table[n - 2 - i][i + 1]
        coeff = coeff * (t + (i + 1))
    return l_n_x


def gauss_first(table, t, idx):
    n = len(table) - 1
    l_n_x = table[idx][0] + t * table[idx][1]
    coeff = t
    for i in range(2, n + 1):
        coeff *= (t + (-1) ** (i + 1) * (i // 2)) / i
        l_n_x += coeff * table[idx - i // 2][i]
    return l_n_x


def R_n(n, a, b, table, x, h):
    x_ = symbols('x')
    f_ = x_ ** 2 - log(x_ / 2) / log(10)
    deriv_n = diff(f_, x_, n + 1)  # Производная n-го порядка

    interv = Interval(a, b)  # Минимум и максимум производной на интервале
    min = (minimum(deriv_n, x_, interv))
    max = (maximum(deriv_n, x_, interv))

    res_min = abs(float(max)) / math.factorial(n + 1) * w(table, x, h, a)  # минимальный и максимальный остаточный член
    res_max = abs(float(min)) / math.factorial(n + 1) * w(table, x, h, a)

    res = [res_min, res_max]
    res.sort()
    return res


def w(table, x, h, a):
    n = len(table)
    res = 1
    for i in range(n - 1):
        res *= x - (a + i * h)
    return res


def print_result(x, f_x, l_n, err, r):
    print(f"x = {x}")
    print(f"f(x) = {f_x:.10f} and L(x) = {l_n:.10f}")
    print(f"Δ = {err:.5e} between {r[0]:.5e} and {r[1]:.5e}?")

    print("Yes") if r[0] < err < r[1] else print("No")


a = 0.5
b = 1
n = 15
x_values = [0.53, 0.98, 0.77]
h = (b - a) / n
table = finite_diffs(n, a, h)

# x**
x = x_values[0]
t = (x - a) / h
f_x = f(x)
l_n = newton_first(table, t)
err = abs(f_x - l_n)

r = R_n(n, a, b, table, x, h)
print_result(x, f_x, l_n, err, r)

print("")

# x***
x = x_values[1]
t = (x - b) / h
f_x = f(x)
l_n = newton_second(table, t)
err = abs(f_x - l_n)

r = R_n(n, a, b, table, x, h)
print_result(x, f_x, l_n, err, r)

print("")

# x****
x = x_values[2]
t = (x - (a + 6.0 * h)) / h
f_x = f(x)
l_n = gauss_first(table, t, 6)
err = abs(f_x - l_n)

r = R_n(n, a, b, table, x, h)
print_result(x, f_x, l_n, err, r)
