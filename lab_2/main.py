import math
import functools
import numpy as np


def f(x):
    return math.log(x**2) + 1 - math.sin(x)

def k(x):
    return x**4 - 10*x**3 + 31*x**2 - 30*x + 9
def g(x):
    return np.sin(x) + np.cos(x)
def dichotomy_method(f, a, b, eps, delta=None, max_iter=1000):
    if delta is None:
        delta = eps / 2

    iterations = 0
    while abs(b - a) > eps and iterations < max_iter:
        iterations += 1
        x1 = (a + b - delta) / 2
        x2 = (a + b + delta) / 2
        if f(x1) <= f(x2):
            b = x2
        else:
            a = x1

    xmin = (a + b) / 2
    return xmin, iterations
def golden_section_method(f, a, b, eps, max_iter=1000):
    phi = (1 + np.sqrt(5)) / 2
    iterations = 0
    x1 = b - (b - a) / phi
    x2 = a + (b - a) / phi
    f_x1 = f(x1)
    f_x2 = f(x2)

    while abs(b - a) > eps and iterations < max_iter:
        iterations += 1
        if f_x1 <= f_x2:
            b = x2
            x2 = x1
            f_x2 = f_x1
            x1 = b - (b - a) / phi
            f_x1 = f(x1)
        else:
            a = x1
            x1 = x2
            f_x1 = f_x2
            x2 = a + (b - a) / phi
            f_x2 = f(x2)

    xmin = (a + b) / 2
    return xmin, iterations

def fibonacci_method(f, a, b, n):
    # Вычисляем последовательность чисел Фибоначчи
    iterations = 0

    fib = [1, 1]
    while fib[-1] < n:
        fib.append(fib[-1] + fib[-2])
    # Вычисляем начальный интервал
    L = (b - a) / fib[-2]
    x1 = a + L * fib[-3]
    x2 = a + L * fib[-2]
    f1, f2 = f(x1), f(x2)
    # Итеративно сужаем интервал
    for i in range(2, len(fib)):
        iterations += 1
        if f1 > f2:
            a, x1 = x1, x2
            x2 = a + L * fib[-i]
            f1, f2 = f2, f(x2)
        else:
            b, x2 = x2, x1
            x1 = a + L * fib[-i-1]
            f1, f2 = f(x1), f1
    return (a + b) / 2, iterations
def parabola_search(f, x1, x2, x3, eps):
    iterations = 0
    f1, f2, f3 = f(x1), f(x2), f(x3)
    while x3 - x1 > eps:
        iterations += 1
        u = x2 - ((x2 - x1)**2*(f2 - f3) - (x2 - x3)**2*(f2 - f1))/(2*((x2 - x1)*(f2 - f3) - (x2 - x3)*(f2 - f1)))
        fu = f(u)

        if x2 <= u:
            if f2 <= fu:
                x1, x2, x3 = x1, x2, u
                f1, f2, f3 = f1, f2, fu
            else:
                x1, x2, x3 = x2, u, x3
                f1, f2, f3 = f2, fu, f3
        else:
            if fu <= f2:
                x1, x2, x3 = x1, u, x2
                f1, f2, f3 = f1, fu, f2
            else:
                x1, x2, x3 = u, x2, x3
                f1, f2, f3 = fu, f2, f3
    return (x1 + x3) / 2, iterations

def brent(f, a, b, tol = 1e-6, max_iter = 500):
    iterations = 0
    golden_ratio = (3 - math.sqrt(5)) / 2
    x = w = v = a + golden_ratio * (b - a)
    fx = fw = fv = f(x)

    d = e = b - a

    for _ in range(max_iter):
        iterations += 1
        g, e = e, d
        tol_1 = tol * abs(x) + tol
        middle = (a + b) / 2

        if abs(x - middle) + (b - a) / 2 <= 2 * tol_1:
            break

        u = None
        if x != w and x != v and w != v and fw != fx and fw != fv and fv != fx:
            # Параболическая аппроксимация
            u = x * (fx - fw) * (fx - fv) / ((fx - fw) * (fx - fv) * (fw - fv)) + \
                w * (fw - fx) * (fw - fv) / ((fw - fx) * (fw - fv) * (fx - fv)) + \
                v * (fv - fx) * (fv - fw) / ((fv - fx) * (fv - fw) * (fx - fw))

            if a <= u <= b and abs(u - x) < g / 2:
                pass
            else:
                u = None

        if u is None:
            if x < middle:
                u = x + golden_ratio * (b - x)
                e = b - x
            else:
                u = x - golden_ratio * (x - a)
                e = x - a

            if abs(u - x) < tol_1:
                u = x + math.copysign(tol_1, u - x)

        d = abs(u - x)
        fu = f(u)

        if fu <= fx:
            if u >= x:
                a = x
            else:
                b = x

            v, w, x = w, x, u
            fv, fw, fx = fw, fx, fu
        else:
            if u >= x:
                b = u
            else:
                a = u

            if fu <= fw or w == x:
                v, w = w, u
                fv, fw = fw, fu
            elif fu <= fv or v == x or v == w:
                v = u
                fv = fu

    return x, iterations, fx


def count_calls(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        wrapper.calls += 1
        return func(*args, **kwargs)
    wrapper.calls = 0
    return wrapper

f_counted = count_calls(f)

def test_method(method, name, *args, **kwargs):
    f_counted.calls = 0
    iterations = method(f_counted, *args, **kwargs)
    print(f"{name}: iterations = {iterations[1]}, function calls = {f_counted.calls}")


# Задаем начальный интервал и точность вычислений
a = 12
b = 16
eps = 1e-6
n = 1000
brent_x_min, iterations, brent_f_min = brent(f, a, b) # метод брента

# Выводим результат
print("Минимум функции по методу дихотомии y=ln(x^2)+1-sin(x) достигается в точке x =", dichotomy_method(f, a, b, eps)[0])
print("Минимум функции по методу золотого сечения y=ln(x^2)+1-sin(x) достигается в точке x =", golden_section_method(f, a, b, eps)[0])
print("Минимум функции по методу Фибоначчи y=ln(x^2)+1-sin(x) достигается в точке x =", fibonacci_method(f, a, b, n)[0])
print("Минимум функции по методу парабол y=ln(x^2)+1-sin(x) достигается в точке x =", parabola_search(f, a, (a+b)/2, b, eps)[0])
print(f"Минимум функции по методу Брента y=ln(x^2)+1-sin(x) достигается в точке x = {brent_x_min}, значение функции в этой точке: {brent_f_min}")
print()
print()


# Выводим результат
print("Минимум функции y=x^4 - 10x^3 + 31x^2 - 30x + 9 при методе дихотомии", dichotomy_method(k, a, b, eps)[0])
print("Минимум функции y=x^4 - 10x^3 + 31x^2 - 30x + 9 при методе золотого сечения", golden_section_method(k, a, b, eps)[0])
print("Минимум функции y=x^4 - 10x^3 + 31x^2 - 30x + 9 при методе Фибоначчи", fibonacci_method(k, a, b, n)[0])
print("Минимум функции y=x^4 - 10x^3 + 31x^2 - 30x + 9 при методе парабол", parabola_search(k, a, (a+b)/2, b, eps)[0])
print("Минимум функции y=x^4 - 10x^3 + 31x^2 - 30x + 9 при методе Брента", brent(k, a, b, eps, n)[0])
print()
print("Минимум функции y=sin(x)+cos(x) при методе дихотомии", dichotomy_method(g, a, b, eps)[0])
print("Минимум функции y=sin(x)+cos(x) при методе золотого сечения", golden_section_method(g, a, b, eps)[0])
print("Минимум функции y=sin(x)+cos(x) при методе Фибоначчи", fibonacci_method(g, a, b, n)[0])
print("Минимум функции y=sin(x)+cos(x) при методе парабол", parabola_search(g, a, (a+b)/2, b, eps)[0])
print("Минимум функции y=sin(x)+cos(x) при методе Брента", brent(g, a, b, eps, n)[0])

print()

eps_list = [1e-2, 1e-4, 1e-6]
methods = [
    (dichotomy_method, "Dichotomy method"),
    (golden_section_method, "Golden section method"),
    (fibonacci_method, "Fibonacci method"),
    (parabola_search, "Parabola search"),
    (brent, "Brent's method")
]

a, b = 13, 16

for eps in eps_list:
    print(f"eps = {eps}:")
    for method, name in methods:
        if method == fibonacci_method:
            n = 30
            test_method(method, name, a, b, n)
        elif method == parabola_search:
            x1, x2, x3 = a, (a + b) / 2, b
            test_method(method, name, x1, x2, x3, eps)
        else:
            test_method(method, name, a, b, eps)
    print()




