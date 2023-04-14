import math
import functools



def f(x):
    return math.log(x**2) + 1 - math.sin(x)
def dichotomy_method(f, a, b, eps):
    """
    Метод дихотомии для минимизации функции.
    f: функция, которую нужно минимизировать
    a: начальная левая граница интервала
    b: начальная правая граница интервала
    eps: точность вычислений
    return: точка минимума функции
    """
    while abs(b - a) > eps:
        x = (a + b) / 2
        y1 = f(x - eps)
        y2 = f(x + eps)
        if y1 < y2:
            b = x
        else:
            a = x
    return (a + b) / 2

def golden_section_method(f, a, b, eps):
    """
    Метод золотого сечения для минимизации функции.
    :param f: функция, которую нужно минимизировать
    :param a: начальная левая граница интервала
    :param b: начальная правая граница интервала
    :param eps: точность вычислений
    :return: точка минимума функции
    """
    phi = (1 + math.sqrt(5)) / 2  # коэффициент золотого сечения
    x1 = b - (b - a) / phi
    x2 = a + (b - a) / phi
    while abs(b - a) > eps:
        if f(x1) < f(x2):
            b = x2
        else:
            a = x1
        x1 = b - (b - a) / phi
        x2 = a + (b - a) / phi
    return (a + b) / 2


def fibonacci_method(f, a, b, n):
    """
    Метод Фибоначчи для минимизации функции.
    :param f: функция, которую нужно минимизировать
    :param a: начальная левая граница интервала
    :param b: начальная правая граница интервала
    :param n: количество точек разбиения интервала
    :return: точка минимума функции
    """
    # Вычисляем последовательность чисел Фибоначчи
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
        if f1 > f2:
            a, x1 = x1, x2
            x2 = a + L * fib[-i]
            f1, f2 = f2, f(x2)
        else:
            b, x2 = x2, x1
            x1 = a + L * fib[-i-1]
            f1, f2 = f(x1), f1
    return (a + b) / 2
def parabola_search(f, x1, x2, x3, eps):
    f1, f2, f3 = f(x1), f(x2), f(x3)
    while x3 - x1 > eps:
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
    return (x1 + x3) / 2

def brent(f, a, b, tol, max_iter):
    golden_ratio = (3 - math.sqrt(5)) / 2
    x = w = v = a + golden_ratio * (b - a)
    fx = fw = fv = f(x)

    d = e = b - a

    for _ in range(max_iter):
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

    return x, fx

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
    result = method(f_counted, *args, **kwargs)
    print(f"{name}: iterations = {method.__name__.count}, function calls = {f_counted.calls}")

eps_list = [1e-2, 1e-4, 1e-6]
methods = [
    (dichotomy_method, "Dichotomy method"),
    (golden_section_method, "Golden section method"),
    (fibonacci_method, "Fibonacci method"),
    (parabola_search, "Parabola search"),
    (brent, "Brent's method")
]

a, b = 0, 2

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
