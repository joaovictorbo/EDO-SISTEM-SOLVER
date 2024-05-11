import numpy as np

import matplotlib.pyplot as plt

def runge_kutta(f, x0, y0, h, n):
    """
    Função para calcular uma EDO utilizando o método de Runge-Kutta de quarta ordem.

    Parâmetros:
    - f: função que define a EDO (dy/dx = f(x, y))
    - x0: valor inicial de x
    - y0: array com os valores iniciais de y
    - h: tamanho do passo
    - n: número de iterações

    Retorna:
    - x: array com os valores de x
    - y: matriz com os valores de y (cada linha representa um valor de y)
    """

    m = len(y0)  # número de valores de y
    x = np.zeros(n+1)
    y = np.zeros((n+1, m))
    x[0] = x0
    y[0] = y0

    for i in range(n):
        k1 = h * f(x[i], y[i])
        k2 = h * f(x[i] + h/2, y[i] + k1/2)
        k3 = h * f(x[i] + h/2, y[i] + k2/2)
        k4 = h * f(x[i] + h, y[i] + k3)

        x[i+1] = x[i] + h
        y[i+1] = y[i] + ((h/6)*(k1 + 2*k2 + 2*k3 + k4))

    return x, y

# Exemplo de uso
def f(x, y):
    return np.array([x**2 + y[0], x**2 + y[1]])  # exemplo com duas equações

x0 = 0
y0 = np.array([1, 2])  # exemplo com dois valores iniciais de y
h = 0.1
n = 10

x, y = runge_kutta(f, x0, y0, h, n)

plt.plot(x, y[:, 0], label='y1')
plt.plot(x, y[:, 1], label='y2')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Solução da EDO')
plt.legend()
plt.grid(True)
plt.show()