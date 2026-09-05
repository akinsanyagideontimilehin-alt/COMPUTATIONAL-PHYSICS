def f(x):
    return 0.04 * (2 * x**3 - 5 * x**2 - 13 * x + 9)

def df(x):
    return 0.04 * (6 * x**2 - 10 * x - 13)

eps = 1e-5

# 1. Метод дихотомии (деления отрезка пополам)
a, b = 0.0, 2.0
steps_dichotomy = 0
while (b - a) / 2 > eps:
    steps_dichotomy += 1
    c = (a + b) / 2
    if f(a) * f(c) < 0:
        b = c
    else:
        a = c
root_dichotomy = (a + b) / 2

# 2. Метод Ньютона (касательных)
x_newton = 0.0  # Начальное приближение
steps_newton = 0
while True:
    steps_newton += 1
    x_next = x_newton - f(x_newton) / df(x_newton)
    if abs(x_next - x_newton) < eps:
        root_newton = x_next
        break
    x_newton = x_next

# 3. Метод хорд (секущих)
x0, x1 = 0.0, 2.0
steps_secant = 0
while abs(x1 - x0) > eps:
    steps_secant += 1
    f_x1, f_x0 = f(x1), f(x0)
    x_next = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0)
    x0, x1 = x1, x_next
root_secant = x1

# 4. Метод простых итераций
# Приведем к виду x = phi(x): 13x = 2x^3 - 5x^2 + 9  =>  x = (2x^3 - 5x^2 + 9) / 13
def phi(x):
    return (2 * x**3 - 5 * x**2 + 9) / 13

x_iter = 1.0  # Начальное приближение
steps_iter = 0
while True:
    steps_iter += 1
    x_next = phi(x_iter)
    if abs(x_next - x_iter) < eps:
        root_iter = x_next
        break
    x_iter = x_next

# 5. Комбинированный метод (хорд + касательных)
a, b = 0.0, 2.0
steps_comb = 0
while (b - a) > eps:
    steps_comb += 1
    a_next = a - f(a) * (b - a) / (f(b) - f(a))
    b_next = b - f(b) / df(b)
    a, b = min(a_next, b_next), max(a_next, b_next)
root_comb = (a + b) / 2

# Вывод результатов
print(f"1. Метод дихотомии: Корень = {root_dichotomy:.5f}, Шагов = {steps_dichotomy}")
print(f"2. Метод Ньютона: Корень = {root_newton:.5f}, Шагов = {steps_newton}")
print(f"3. Метод хорд: Корень = {root_secant:.5f}, Шагов = {steps_secant}")
print(f"4. Метод итераций: Корень = {root_iter:.5f}, Шагов = {steps_iter}")
print(f"5. Комбинированный метод: Корень = {root_comb:.5f}, Шагов = {steps_comb}")
