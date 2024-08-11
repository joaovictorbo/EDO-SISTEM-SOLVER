import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

muw0 = 1.0  # Viscosidade inicial sem polimero
muo = 4.0
mug = 0.25

def Du(u, v, c): # dD/du
    w = 1 - u - v
    return 2 * u / muw(c) - 2 * w / mug

def Dv(u, v, c): # dD/dv
    w = 1 - u - v
    return 2 * v / muo - 2 * w / mug

def Dz(u, v, c): # dD/dc
    return -muwc(c) * u**2 / muw(c)**2

def muw(c): # Viscosidade da água
    return muw0 + c

def muwc(c): # dmuw/dc
    return 1.0

def D(u, v, c): # Denominador
    w = 1 - u - v
    return u**2 / muw(c) + v**2 / muo + w**2 / mug
def a(z):
    return np.sqrt(z)
# Definição das funções f e g
def f(u, v, z):
    return (u**2 / muw(z)) / D(u, v, z)

def g(u, v, z):
    return (v**2 / muo) / D(u, v, z)

# Derivadas parciais de f e g
def df_du(u, v, z):
    return (2 * u / muw(z) * D(u, v, z) - u**2 / muw(z) * Du(u, v, z)) / (D(u, v, z)**2)

def df_dv(u, v, z):
    return (-u**2 / muw(z) * Dv(u, v, z)) / (D(u, v, z)**2)

def df_dz(u, v, z):
    return u**2 * (-(muwc(z) / muw(z)**2) * D(u, v, z) - Dz(u, v, z) / muw(z)) / (D(u, v, z)**2)

def dg_du(u, v, z):
    return (-v**2 / muo) * Du(u, v, z) / (D(u, v, z)**2)

def dg_dv(u, v, z):
    return (2 * v / muo * D(u, v, z) - v**2 / muo * Dv(u, v, z)) / (D(u, v, z)**2)

def dg_dz(u, v, z):
    return -v**2 / muo * Dz(u, v, z) / (D(u, v, z)**2)

def F(u, v, z, f0, v0, g0, u0):
    f_value = f(u, v, z)
    g_value = g(u, v, z)
    return (f_value - f0) * (v - v0) - (g_value - g0) * (u - u0)

def G(u, v, z, f0, z0, u0, a, alpha):
    f_value = f(u, v, z)
    a_z = a(z)
    a_z0 = a(z0)
    return (f_value - f0) * (u0 * (z - z0) + alpha * (a_z - a_z0)) - f0 * (z - z0) * (u - u0)

# Derivadas parciais de F e G
def dF_du(u, v, z, f0, v0, g0, u0):
    return df_du(u, v, z) * (v - v0) - dg_du(u, v, z) * (u - u0) - (g(u, v, z) - g0)

def dF_dv(u, v, z, f0, v0, g0, u0):
    return df_dv(u, v, z) * (v - v0) - dg_dv(u, v, z) * (u - u0) + (f(u, v, z) - f0)

def dF_dz(u, v, z, f0, v0, g0, u0):
    return df_dz(u, v, z) * (v - v0) - dg_dz(u, v, z) * (u - u0)

def dG_du(u, v, z, f0, z0, u0, a, alpha):
    return df_du(u, v, z) * (u0 * (z - z0) + alpha * (a(z) - a(z0))) - f0 * (z - z0) - dg_du(u, v, z) * (u - u0)

def dG_dv(u, v, z, f0, z0, u0, a, alpha):
    return df_dv(u, v, z) * (u0 * (z - z0) + alpha * (a(z) - a(z0))) - dg_dv(u, v, z) * (u - u0)

def dG_dz(u, v, z, f0, z0, u0, a, alpha):
    return df_dz(u, v, z) * (u0 * (z - z0) + alpha * (a(z) - a(z0))) + (f(u, v, z) - f0) * (u0 + alpha * a(z)) - f0 * (u - u0)
def det_z(u, v, z, f0, v0, g0, u0):
    return np.linalg.det([
        [dF_du(u, v, z, f0, v0, g0, u0), dF_dv(u, v, z, f0, v0, g0, u0)],
        [dG_du(u, v, z, f0, z0, u0, a, alpha), dG_dv(u, v, z, f0, z0, u0, a, alpha)]
    ])

def det_v(u, v, z, f0, v0, g0, u0):
    return np.linalg.det([
        [dF_du(u, v, z, f0, v0, g0, u0), dF_dz(u, v, z, f0, v0, g0, u0)],
        [dG_du(u, v, z, f0, z0, u0, a, alpha), dG_dz(u, v, z, f0, z0, u0, a, alpha)]
    ])

def det_u(u, v, z, f0, v0, g0, u0):
    return np.linalg.det([
        [dF_dv(u, v, z, f0, v0, g0, u0), dF_dz(u, v, z, f0, v0, g0, u0)],
        [dG_dv(u, v, z, f0, z0, u0, a, alpha), dG_dz(u, v, z, f0, z0, u0, a, alpha)]
    ])

def system_with_determinants(s, y, f0, v0, g0, u0, z0, a, alpha):
    u, v, z = y

    det_u_value = det_u(u, v, z, f0, v0, g0, u0)
    det_v_value = det_v(u, v, z, f0, v0, g0, u0)
    det_z_value = det_z(u, v, z, f0, v0, g0, u0)

    du_ds = det_u_value
    dv_ds = -det_v_value
    dz_ds = det_z_value

    return [du_ds, dv_ds, dz_ds]

# Condições iniciais
u0 = 0.1  # Exemplo
v0 = 0.6  # Exemplo
z0 = 0.2  # Exemplo
f0 = f(u0, v0, z0)
g0 = g(u0, v0, z0)
alpha = 10**-3  # Exemplo
y0 = [u0, v0, z0]

y0 = [0.10010735741660978, 0.5999219526580091, 0.2215289621825406]

# Intervalo de integração
s_span = (0, 0.05)

# Resolução do sistema
sol = solve_ivp(system_with_transition, s_span, y0, method='RK45', t_eval=np.linspace(0, 0.05, 2))

# Verificação dos resultados
if sol.y.shape[0] != 3:
    raise ValueError("A solução não retornou os valores esperados.")

# Plot dos resultados em 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Adicionar o triângulo que vai de c=0 a c=1
vertices = np.array([
    [0, 0, 0], [1, 0, 0], [0, 1, 0],
    [0, 0, 1], [1, 0, 1], [0, 1, 1]
])
edges = [
    (vertices[0], vertices[1]),
    (vertices[0], vertices[2]),
    (vertices[0], vertices[3]),
    (vertices[1], vertices[2]),
    (vertices[1], vertices[4]),
    (vertices[2], vertices[5]),
    (vertices[3], vertices[4]),
    (vertices[3], vertices[5]),
    (vertices[4], vertices[5])
]

for edge in edges:
    ax.plot(*zip(*edge), color='black')

ax.plot(sol.y[0], sol.y[1], sol.y[2], label='Trajetória')
ax.set_xlabel('u(s)')
ax.set_ylabel('v(s)')
ax.set_zlabel('z(s)')
ax.legend()
plt.title('Solução do sistema de EDOs em 3D')
plt.show()
