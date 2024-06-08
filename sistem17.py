import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

muw0 = 1.0  # Viscosidade inicial sem polimero
muo = 4.0  # 9.5
mug = 0.25  # 0.45

def Du(u, v, c): # dD/du
    w = 1 - u - v
    if mug == 0 or muw(c) == 0:
        return np.nan
    return 2 * u / muw(c) - 2 * w / mug

def Dv(u, v, c): # dD/dv
    w = 1 - u - v
    if mug == 0 or muo == 0:
        return np.nan
    return 2 * v / muo - 2 * w / mug

def Dz(u, v, c): # dD/dc
    if muw(c) == 0:
        return np.nan
    return -muwc(c) * u**2 / muw(c)**2

def muw(c): # Viscosidade da água
    return muw0 + c

def muwc(c): # dmuw/dc
    return 1.0

def D(u, v, c): # Denominador
    w = 1 - u - v
    return u**2 / muw(c) + v**2 / muo + w**2 / mug

# Definição das funções F e G
def F(u, v, z):
    return (u**2 / muw(z)) / D(u, v, z)

def G(u, v, z):
    return (v**2 / muo) / D(u, v, z)

# Derivadas parciais de F e G
def dF_du(u, v, z):
    denom = D(u, v, z)
    if denom == 0 or muw(z) == 0:
        return np.nan
    return (2 * u / muw(z) * denom - u**2 / muw(z) * Du(u, v, z)) / (denom**2)

def dF_dv(u, v, z):
    denom = D(u, v, z)
    if denom == 0 or muw(z) == 0:
        return np.nan
    return (-u**2 / muw(z) * Dv(u, v, z)) / (denom**2)

def dF_dz(u, v, z):
    denom = D(u, v, z)
    if denom == 0 or muw(z) == 0:
        return np.nan
    return u**2 * (-(muwc(z) / muw(z)**2) * denom - Dz(u, v, z) / muw(z)) / (denom**2)

def dG_du(u, v, z):
    denom = D(u, v, z)
    if denom == 0 or muo == 0:
        return np.nan
    return (-v**2 / muo) * Du(u, v, z) / (denom**2)

def dG_dv(u, v, z):
    denom = D(u, v, z)
    if denom == 0 or muo == 0:
        return np.nan
    return (2 * v / muo * denom - v**2 / muo * Dv(u, v, z)) / (denom**2)

def dG_dz(u, v, z):
    denom = D(u, v, z)
    if denom == 0:
        return np.nan
    return -v**2 / muo * Dz(u, v, z) / (denom**2)

# Definição do sistema de EDOs
def system(s, y):
    u, v, z = y
    # Calculando os determinantes
    det1 = np.linalg.det([[dF_dv(u, v, z), dF_dz(u, v, z)],
                          [dG_dv(u, v, z), dG_dz(u, v, z)]])
    
    det2 = np.linalg.det([[dF_du(u, v, z), dF_dz(u, v, z)],
                          [dG_du(u, v, z), dG_dz(u, v, z)]])
    
    det3 = np.linalg.det([[dF_du(u, v, z), dF_dv(u, v, z)],
                          [dG_du(u, v, z), dG_dv(u, v, z)]])
    
    if det2 == 0 or det3 == 0:
        return [np.nan, np.nan, np.nan]
    
    du_ds = det1 / det2
    dv_ds = -det2 / det3
    dz_ds = 1
    
    return [du_ds, dv_ds, dz_ds]

# Condições iniciais
u0 = 0.5  # Exemplo
v0 = 0.5  # Exemplo
z0 = 0  # Exemplo
y0 = [u0, v0, z0]

# Intervalo de integração
s_span = (0, 1)

# Resolução do sistema
sol = solve_ivp(system, s_span, y0, method='RK45', t_eval=np.linspace(0, 1, 10))

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
