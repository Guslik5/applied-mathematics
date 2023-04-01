import numpy as np
import matplotlib.pyplot as plt

# Определим функции
def f(x):
    return x**2

def g(x):
    return np.sin(x)

# Вычислим аналитические производные
def df(x):
    return 2*x

def dg(x):
    return np.cos(x)

# определим функции для интегралов
def f2(x):
    return x

def g2(x):
    return x**2

# определяем функции для вычисления численных производных
def left_diff(f, x, h):
    return (f(x) - f(x-h))/h

def right_diff(f, x, h):
    return (f(x+h) - f(x))/h

def central_diff(f, x, h):
    return (f(x+h) - f(x-h))/(2*h)

# Функции для определения определенного интеграла на отрезке
def left_rectangle_method(f, a, b, n):
    """Вычисление интеграла функции f(x) на отрезке [a,b] с использованием метода левых прямоугольников"""

    h = (b - a) / n  # длина интервала
    x = a  # начальное значение x
    integral = 0  # значение интеграла
    for i in range(n):
        integral += f(x) * h  # добавляем площадь прямоугольника
        x += h  # переходим к следующей точке
    return integral

def right_rectangle_method(f, a, b, n):
    """Вычисление интеграла функции f(x) на отрезке [a,b] с использованием метода правых прямоугольников"""

    dx = (b - a) / n  # длина интервала
    x = a + dx  # начальное значение x
    integral = 0  # значение интеграла
    for i in range(n):
        integral += f(x) * dx  # добавляем площадь прямоугольника
        x += dx  # переходим к следующей точке
    return integral

def midpoint_rectangle_method(f, a, b, n):
    """Вычисление интеграла функции f(x) на отрезке [a,b] с использованием метода средних прямоугольников"""

    dx = (b - a) / n  # длина интервала
    x = a + dx / 2  # начальное значение x (смещаем на половину интервала)
    integral = 0  # значение интеграла
    for i in range(n):
        integral += f(x) * dx  # добавляем площадь прямоугольника
        x += dx  # переходим к следующей точке
    return integral

def trapezoidal_rule(f, a, b, n):
    """Вычисление интеграла функции f(x) на отрезке [a,b] с использованием формулы трапеций"""

    dx = (b - a) / n  # длина интервала
    x = a + dx  # начальное значение x
    integral = (f(a) + f(b)) / 2  # значение интеграла (суммируем значения на краях отрезка)
    for i in range(1, n):
        integral += f(x)  # добавляем значение функции в текущей точке
        x += dx  # переходим к следующей точке
    integral *= dx  # умножаем на длину интервала
    return integral

def simpson_method(f, a, b, n):
    """
    Вычисление интеграла функции f(x) на интервале [a, b] с использованием метода Симпсона.
    """
    h = (b - a) / n
    s = f(a) + f(b)

    for i in range(1, n):
        x = a + i * h
        if i % 2 == 0:
            s += 2 * f(x)
        else:
            s += 4 * f(x)

    return h / 3 * s

# Массив с уменьшением шага(в 2,4,8 и 16 раз)
h_values = [0.1, 0.05, 0.025, 0.0125, 0.00625]

# Массивы для хранения среднеквадратичных отклонений
df_rms_left_values = []
df_rms_right_values = []
df_rms_central_values = []
dg_rms_left_values = []
dg_rms_right_values = []
dg_rms_central_values = []

for h in h_values:
    x_grid = np.arange(-np.pi, np.pi+h, h)

    # вычисляем значения аналитических производных в узлах сетки
    f_values = f(x_grid)
    g_values = g(x_grid)

    # вычисляем значения аналитических производных в узлах сетки
    df_values = df(x_grid)
    dg_values = dg(x_grid)

    # вычисляем значения численных производных в узлах сетки
    df_left_values = left_diff(f, x_grid, h)
    dg_left_values = left_diff(g, x_grid, h)

    df_right_values = right_diff(f, x_grid, h)
    dg_right_values = right_diff(g, x_grid, h)

    df_central_values = central_diff(f, x_grid, h)
    dg_central_values = central_diff(g, x_grid, h)

    # вычисляем значения среднеквадратичного отклонениячисленных производных от аналитических
    df_rms_left = np.sqrt(np.mean((df_left_values - df_values)**2))
    dg_rms_left = np.sqrt(np.mean((dg_left_values - dg_values)**2))

    df_rms_right = np.sqrt(np.mean((df_right_values - df_values)**2))
    dg_rms_right = np.sqrt(np.mean((dg_right_values - dg_values)**2))

    df_rms_central = np.sqrt(np.mean((df_central_values - df_values)**2))
    dg_rms_central = np.sqrt(np.mean((dg_central_values - dg_values)**2))

    #print("Левые разностные производные функции f(x) в узлах сетки:\n", df_left_values)
    #print("Правые разностные производные функции f(x) в узлах сетки:\n", df_right_values)
    #print("Центральные разностные производные функции f(x) в узлах сетки:\n", df_central_values)

    #print("Левые разностные производные функции g(x) в узлах сетки:\n", dg_left_values)
    #print("Правые разностные производные функции g(x) в узлах сетки:\n", dg_right_values)
    #print("Центральные разностные производные функции g(x) в узлах сетки:\n", dg_central_values)


    # Добавление значений в массивы
    df_rms_left_values.append(df_rms_left)
    df_rms_right_values.append(df_rms_right)
    df_rms_central_values.append(df_rms_central)
    dg_rms_left_values.append(dg_rms_left)
    dg_rms_right_values.append(dg_rms_right)
    dg_rms_central_values.append(dg_rms_central)

    # вывод результатов
    print("Шаг сетки:", h)
    print("RMS ошибка для функции x^2:")
    print("Левая производная:", df_rms_left)
    print("Правая производная:", df_rms_right)
    print("Центральная производная:", df_rms_central)

    print("RMS ошибка для функции sin(x):")
    print("Левая производная:", dg_rms_left)
    print("Правая производная:", dg_rms_right)
    print("Центральная производная:", dg_rms_central)
    print()


# Вторая часть работы
a, b = 0, 1  # Интервал
n = 1000  # количество разбиений
n_values = [1000, 2000, 4000, 8000, 16000]  # Массив для уменьшения шага

# массивы для хранения определенного интеграла
if_rms_left_values = []
if_rms_right_values = []
if_rms_central_values = []
if_rms_trap_values = []
if_rms_simpson_values = []

ig_rms_left_values = []
ig_rms_right_values = []
ig_rms_central_values = []
ig_rms_trap_values = []
ig_rms_simpson_values = []

for n in n_values:
    # Вычисляем определенный интеграл функции x метод симпсона
    print("\nМетод Симпсона")
    integral_x = simpson_method(f2, a, b, n)
    print("Определенный интеграл для функции x на отрезке [0;1]:", integral_x)
    # Сравниваем с аналитическим ответом
    analytic_x = 1/2
    difference_x_simpson = np.sqrt(np.mean((integral_x - analytic_x) ** 2))
    print("Разность с аналитическим ответом:", difference_x_simpson)


    # Вычисляем определенный интеграл функции x^2
    integral_x2 = simpson_method(g2, a, b, n)
    print("Определенный интеграл для функции x^2 на отрезке [0;1]:", integral_x2)
    # Сравниваем с аналитическим ответом
    analytic_x2 = 1/3
    difference_x2_simpson = np.sqrt(np.mean((integral_x2 - analytic_x2) ** 2))
    print("Разность с аналитическим ответом:", difference_x2_simpson)

    print("\nФормула Трапеции")
    # Вычисляем определенный интеграл функции x
    integral_x = trapezoidal_rule(f2, a, b, n)
    print("Определенный интеграл для функции x на отрезке [0;1]:", integral_x)
    # Сравниваем с аналитическим ответом
    analytic_x = 1/2
    difference_x_trap = np.sqrt(np.mean((integral_x - analytic_x) ** 2))
    print("Разность с аналитическим ответом:", difference_x_trap)

    # Вычисляем определенный интеграл функции x^2
    integral_x2 = trapezoidal_rule(g2, a, b, n)
    print("Определенный интеграл для функции x^2 на отрезке [0;1]:", integral_x2)
    # Сравниваем с аналитическим ответом
    analytic_x2 = 1/3
    difference_x2_trap = np.sqrt(np.mean((integral_x2 - analytic_x2) ** 2))
    print("Разность с аналитическим ответом:", difference_x2_trap)


    print("\nСредних Прямоугольников")
    # Вычисляем определенный интеграл функции x
    integral_x = midpoint_rectangle_method(f2, a, b, n)
    print("Определенный интеграл для функции x на отрезке [0;1]:", integral_x)
    # Сравниваем с аналитическим ответом
    analytic_x = 1/2
    difference_x_mid = np.sqrt(np.mean((integral_x - analytic_x) ** 2))
    print("Разность с аналитическим ответом:", difference_x_mid)

    # Вычисляем определенный интеграл функции x^2
    integral_x2 = midpoint_rectangle_method(g2, a, b, n)
    print("Определенный интеграл для функции x^2 на отрезке [0;1]:", integral_x2)
    # Сравниваем с аналитическим ответом
    analytic_x2 = 1/3
    difference_x2_mid = np.sqrt(np.mean((integral_x2 - analytic_x2) ** 2))
    print("Разность с аналитическим ответом:", difference_x2_mid)

    print("\nПравых прямоугольников")
    # Вычисляем определенный интеграл функции x
    integral_x = right_rectangle_method(f2, a, b, n)
    print("Определенный интеграл для функции x на отрезке [0;1]:", integral_x)
    # Сравниваем с аналитическим ответом
    analytic_x = 1/2
    difference_x_right = np.sqrt(np.mean((integral_x - analytic_x) ** 2))
    print("Разность с аналитическим ответом:", difference_x_right)

    # Вычисляем определенный интеграл функции x^2
    integral_x2 = right_rectangle_method(g2, a, b, n)
    print("Определенный интеграл для функции x^2 на отрезке [0;1]:", integral_x2)
    # Сравниваем с аналитическим ответом
    analytic_x2 = 1/3
    difference_x2_right = np.sqrt(np.mean((integral_x2 - analytic_x2) ** 2))
    print("Разность с аналитическим ответом:", difference_x2_right)

    print("\nЛевых прямоугольников")
    # Вычисляем определенный интеграл функции x
    integral_x = left_rectangle_method(f2, a, b, n)
    print("Определенный интеграл для функции x на отрезке [0;1]:", integral_x)
    # Сравниваем с аналитическим ответом
    analytic_x = 1/2
    difference_x_left = np.sqrt(np.mean((integral_x - analytic_x) ** 2))
    print("Разность с аналитическим ответом:", difference_x_left)

    # Вычисляем определенный интеграл функции x^2
    integral_x2 = left_rectangle_method(g2, a, b, n)
    print("Определенный интеграл для функции x^2 на отрезке [0;1]:", integral_x2)
    # Сравниваем с аналитическим ответом
    analytic_x2 = 1/3
    difference_x2_left = np.sqrt(np.mean((integral_x2 - analytic_x2) ** 2))
    print("Разность с аналитическим ответом:", difference_x2_left)

    # добавляем значения в массивы
    if_rms_left_values.append(difference_x_left)
    if_rms_right_values.append(difference_x_right)
    if_rms_central_values.append(difference_x_mid)
    if_rms_trap_values.append(difference_x_trap)
    if_rms_simpson_values.append(difference_x_simpson)

    ig_rms_left_values.append(difference_x2_left)
    ig_rms_right_values.append(difference_x2_right)
    ig_rms_central_values.append(difference_x2_mid)
    ig_rms_trap_values.append(difference_x2_trap)
    ig_rms_simpson_values.append(difference_x2_simpson)


# Построим графики функций и производных
x = np.linspace(-np.pi, np.pi, 100)
fig, ax = plt.subplots(2, 2, figsize=(8, 6))

ax[0, 0].plot(x, f(x), label=r'$f(x) = x^2$')
ax[0, 0].set_title('Функция f(x)')
ax[0, 0].set_xlabel('x')
ax[0, 0].set_ylabel('f(x)')
ax[0, 0].legend()

ax[0, 1].plot(x, df(x), label=r"$f'(x) = 2x$")
ax[0, 1].set_title('Производная f(x)')
ax[0, 1].set_xlabel('x')
ax[0, 1].set_ylabel("f'(x)")
ax[0, 1].legend()

ax[1, 0].plot(x, g(x), label=r"$g(x) = \sin(x)$")
ax[1, 0].set_title('Функция g(x)')
ax[1, 0].set_xlabel('x')
ax[1, 0].set_ylabel('g(x)')
ax[1, 0].legend()

ax[1, 1].plot(x, dg(x), label=r"$g'(x) = \cos(x)$")
ax[1, 1].set_title('Производная g(x)')
ax[1, 1].set_xlabel('x')
ax[1, 1].set_ylabel("g'(x)")
ax[1, 1].legend()

plt.tight_layout()
plt.show()


# Построим графики зависимости среднеквадратичного отклонения от величины шага
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

ax1.plot(h_values, df_rms_left_values, label='Left difference')
ax1.plot(h_values, df_rms_right_values, label='Right difference')
ax1.plot(h_values, df_rms_central_values, label='Central difference')
ax1.set_xlabel('Step size, h')
ax1.set_ylabel('RMS error')
ax1.set_title('RMS error for x^2')
ax1.legend()

ax2.plot(h_values, dg_rms_left_values, label='Left difference')
ax2.plot(h_values, dg_rms_right_values, label='Right difference')
ax2.plot(h_values, dg_rms_central_values, label='Central difference')
ax2.set_xlabel('Step size, h')
ax2.set_ylabel('RMS error')
ax2.set_title('RMS error for sin(x)')
ax2.legend()

plt.tight_layout()
plt.show()


# Построим графики зависимости отклонения от величены шага
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

ax1.plot(n_values, if_rms_left_values, label='Left difference')
ax1.plot(n_values, if_rms_right_values, label='Right difference')
ax1.plot(n_values, if_rms_central_values, label='Central difference')
ax1.plot(n_values, if_rms_trap_values, label='Trap difference')
ax1.plot(n_values, if_rms_simpson_values, label='Simpson difference')
ax1.set_xlabel('Step size, n')
ax1.set_ylabel('RMS error x')
ax1.set_title('RMS error for x')
ax1.legend()

ax2.plot(n_values, ig_rms_left_values, label='Left difference')
ax2.plot(n_values, ig_rms_right_values, label='Right difference')
ax2.plot(n_values, ig_rms_central_values, label='Central difference')
ax2.plot(n_values, ig_rms_trap_values, label='Trap difference')
ax2.plot(n_values, ig_rms_simpson_values, label='Simpson difference')
ax2.set_xlabel('Step size, n')
ax2.set_ylabel('RMS error x')
ax2.set_title('RMS error for x^2')
ax2.legend()

plt.tight_layout()
plt.show()