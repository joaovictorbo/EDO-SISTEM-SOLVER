import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Definindo os parâmetros da esfera
R = 1

# Definindo os parâmetros do plano
a, b, c, d = 1, 1, 1, 0

# Função da esfera
def esfera(x, y, z):
    return x**2 + y**2 + z**2 - R**2

# Função do plano
def plano(x, y, z):
    return a*x + b*y + c*z - d

# Função para encontrar a interseção
def intersecao(vars):
    x, y = vars
    z = fsolve(lambda z: esfera(x, y, z) + plano(x, y, z), 0)[0]
    return [plano(x, y, z), esfera(x, y, z)]

# Geração de pontos na curva de interseção
theta = np.linspace(0, 2 * np.pi, 100)
x_vals = R * np.cos(theta)
y_vals = R * np.sin(theta)
z_vals = [fsolve(lambda z: esfera(x, y, z) + plano(x, y, z), 0)[0] for x, y in zip(x_vals, y_vals)]

# Plotando as superfícies e a curva de interseção
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Superfície da esfera
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = R * np.outer(np.cos(u), np.sin(v))
y = R * np.outer(np.sin(u), np.sin(v))
z = R * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='b', alpha=0.3)

# Superfície do plano
xx, yy = np.meshgrid(np.linspace(-1.5, 1.5, 100), np.linspace(-1.5, 1.5, 100))
zz = (d - a*xx - b*yy) / c
ax.plot_surface(xx, yy, zz, color='r', alpha=0.3)

# Curva de interseção
ax.plot(x_vals, y_vals, z_vals, color='k', linewidth=2)

# Configurações adicionais do gráfico
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Curva de Interseção entre Esfera e Plano')

plt.show()
