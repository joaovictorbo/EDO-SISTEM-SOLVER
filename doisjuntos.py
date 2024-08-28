import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from system import system

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

def muw(c):  # Viscosidade da agua
    return muw0 * 2**c

def muwc(c):  # dmuw/dc
    return np.log(2) * muw0 * 2**c

def D(u, v, c): # Denominador
    w = 1 - u - v
    return u**2 / muw(c) + v**2 / muo + w**2 / mug

def a(z):
    return np.sqrt(z)
def da_dz(z):
    return 1 / (2 * np.sqrt(z))

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
    return df_du(u, v, z) * (u0 * (z - z0) + (alpha * (a(z) - a(z0))))

def dG_dv(u, v, z, f0, z0, u0, a, alpha):
    return df_dv(u, v, z) * ((u0 * (z - z0)) + (alpha * (a(z) - a(z0)))) - (dg_dv(u, v, z) * (u - u0))

def dG_dz(u, v, z, f0, z0, u0, a, alpha):
    return df_dz(u, v, z) * ((u0 * (z - z0)) + (alpha * (a(z) - a(z0)))) + (f(u, v, z) - f0) * (u0 + (alpha * da_dz(z))) - (f0 * (u - u0))

def det_z(u, v, z):
    return np.linalg.det([
        [dF_du(u, v, z, f0, v0, g0, u0), dF_dv(u, v, z, f0, v0, g0, u0)],
        [dG_du(u, v, z, f0, z0, u0, a, alpha), dG_dv(u, v, z, f0, z0, u0, a, alpha)]
    ])

def det_v(u, v, z):
    return np.linalg.det([
        [dF_du(u, v, z, f0, v0, g0, u0), dF_dz(u, v, z, f0, v0, g0, u0)],
        [dG_du(u, v, z, f0, z0, u0, a, alpha), dG_dz(u, v, z, f0, z0, u0, a, alpha)]
    ])

def det_u(u, v, z):
    return np.linalg.det([
        [dF_dv(u, v, z, f0, v0, g0, u0), dF_dz(u, v, z, f0, v0, g0, u0)],
        [dG_dv(u, v, z, f0, z0, u0, a, alpha), dG_dz(u, v, z, f0, z0, u0, a, alpha)]
    ])

def system_with_determinants(s, y):
    u, v, z = y
    h = (det_u(u, v, z)**2 + det_v(u, v, z)**2 + det_z(u, v, z)**2)**0.5
    du_ds = det_u(u, v, z)/h
    dv_ds = -det_v(u, v, z)/h
    dz_ds = det_z(u, v, z)/h
    return [du_ds, dv_ds, dz_ds]

# Condições iniciais
u0 = 0.5
v0 = 0.4
z0 = 0.8
f0 = f(u0, v0, z0)
g0 = g(u0, v0, z0)
alpha = 10**-3

# Integração para obter os valores iniciais positivos e negativos
getinicialvalues = solve_ivp(system, (0, 0.2), [u0, v0, z0], method='LSODA', t_eval=np.linspace(0, 0.1, 4))
u1, v1, z1 = getinicialvalues.y[:, -1]
y1 = [u1, v1, z1]

getinicialvaluesneg = solve_ivp(system, (0, -0.2), [u0, v0, z0], method='LSODA', t_eval=np.linspace(0, -0.1, 4))
u1neg, v1neg, z1neg = getinicialvaluesneg.y[:, -1]
y1neg = [u1neg, v1neg, z1neg]

# Intervalo de integração positivo e negativo
s_span = (0, 10)
s_span2 = (0, -10)

# Integração do sistema com determinantes
sol = solve_ivp(system_with_determinants, s_span, y1, method='LSODA', t_eval=np.linspace(0, 10, 2000))
getinicialvalues.y = np.hstack((getinicialvalues.y, sol.y))

sol2 = solve_ivp(system_with_determinants, s_span2, y1neg, method='LSODA', t_eval=np.linspace(0, -10, 2000))
getinicialvaluesneg.y = np.hstack((getinicialvaluesneg.y, sol2.y))

t_span = (0,10)
sol3 = solve_ivp(system, t_span, [u0, v0, z0], method='LSODA', t_eval=np.linspace(0,10,2000))
t_span2 = (0,-10)
sol4 = solve_ivp(system, t_span2, [u0, v0, z0], method='LSODA', t_eval=np.linspace(0,-10,2000))

# Função para verificar se um ponto está dentro do triângulo
def dentro_do_triangulo(u, v, c):
    return u >= 0 and v >= 0 and u + v <= 1 and 0 <= c <= 1

# Função para dividir as trajetórias ao entrar e sair do triângulo
def dividir_trajetorias(sol):
    trajetorias = []
    traj_atual = []
    dentro = False
    
    for i in range(len(sol[0])):
        u, v, c = sol[0][i], sol[1][i], sol[2][i]
        if dentro_do_triangulo(u, v, c):
            if not dentro:
                if traj_atual:
                    trajetorias.append(traj_atual)
                traj_atual = []
            dentro = True
            traj_atual.append([u, v, c])
        else:
            if dentro:
                trajetorias.append(traj_atual)
                traj_atual = []
            dentro = False
    
    if traj_atual:
        trajetorias.append(traj_atual)
    
    return trajetorias

# Dividir as trajetórias
trajetorias1 = dividir_trajetorias(getinicialvalues.y)
trajetorias2 = dividir_trajetorias(getinicialvaluesneg.y)
trajetorias3 = dividir_trajetorias(sol3.y)
trajetorias4 = dividir_trajetorias(sol4.y)

# Plotar as trajetórias divididas
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for traj in trajetorias1:
    traj = np.array(traj)
    ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], color="red")

for traj in trajetorias2:
    traj = np.array(traj)
    ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], color="red")
for traj in trajetorias3:
    traj = np.array(traj)
    ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], color="blue")
for traj in trajetorias4:
    traj = np.array(traj)
    ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], color="blue")


ax.set_xlabel('u(s)')
ax.set_ylabel('v(s)')
ax.set_zlabel('c(s)')
ax.legend()

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
ax.scatter(u0, v0, z0, color='red', s=10, label='Ponto Inicial', edgecolor='black')

plt.title('Solução do sistema de EDOs em 3D')
plt.show()
