import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Funções dadas
def muw(c):
    return 1.0 + c

def muo():
    return 4.0

def mug():
    return 0.25

def D(u, v, c):
    w = 1 - u - v
    return u**2 / muw(c) + v**2 / muo() + w**2 / mug()

def f(u, v, c):
    return (u**2 / muw(c)) / D(u, v, c)

def g(u, v, c):
    return (v**2 / muo()) / D(u, v, c)

# Derivadas parciais de f e g
def df_du(u, v, c):
    return (2 * u / muw(c) * D(u, v, c) - u**2 / muw(c) * D_du(u, v, c)) / (D(u, v, c)**2)

def df_dv(u, v, c):
    return (-u**2 / muw(c) * D_dv(u, v, c)) / (D(u, v, c)**2)

def df_dc(u, v, c):
    return u**2 * (-muwc(c) / muw(c)**2 * D(u, v, c) - D_dc(u, v, c) / muw(c)) / (D(u, v, c)**2)

def dg_du(u, v, c):
    return (-v**2 / muo() * D_du(u, v, c)) / (D(u, v, c)**2)

def dg_dv(u, v, c):
    return (2 * v / muo() * D(u, v, c) - v**2 / muo() * D_dv(u, v, c)) / (D(u, v, c)**2)

def dg_dc(u, v, c):
    return -v**2 / muo() * D_dc(u, v, c) / (D(u, v, c)**2)

def D_du(u, v, c):
    w = 1 - u - v
    return 2 * u / muw(c) - 2 * w / mug()

def D_dv(u, v, c):
    w = 1 - u - v
    return 2 * v / muo() - 2 * w / mug()

def D_dc(u, v, c):
    return -u**2 * muwc(c) / muw(c)**2

def muwc(c):
    return 1.0

# Sistema de EDOs
def system(s, y):
    u, v, z = y
    p = (df_dc(u, v, z) - df_du(u, v, z) * g(u, v, z) + df_dv(u, v, z) * f(u, v, z)) / (df_dc(u, v, z) * g(u, v, z) - df_dv(u, v, z) * f(u, v, z))
    q = (dg_dc(u, v, z) - dg_du(u, v, z) * g(u, v, z) + dg_dv(u, v, z) * f(u, v, z)) / (dg_dc(u, v, z) * g(u, v, z) - dg_dv(u, v, z) * f(u, v, z))
    return [p, q, 1]

# Condições iniciais
u0, v0, z0 = 0.4, 0.5, 0.3
y0 = [u0, v0, z0]

# Intervalo de integração
s_span = (0, 1)

# Resolução do sistema
sol = solve_ivp(system, s_span, y0, method='RK45', t_eval=np.linspace(0, 1, 100))

# Verificação dos resultados
if sol.y.shape[0] != 3:
    raise ValueError("A solução não retornou os valores esperados.")

# Plot dos resultados em 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(sol.y[0], sol.y[1], sol.y[2], label='Trajetória')
ax.set_xlabel('u(s)')
ax.set_ylabel('v(s)')
ax.set_zlabel('z(s)')
ax.legend()
plt.title('Solução do sistema de EDOs em 3D')
plt.show()
