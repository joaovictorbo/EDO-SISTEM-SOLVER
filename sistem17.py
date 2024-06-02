import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Definição das funções F e G
def F(u, v, z):
    return u * v - z

def G(u, v, z):
    return u + v + z

# Derivadas parciais de F e G
def dF_du(u, v, z):
    return v

def dF_dv(u, v, z):
    return u

def dF_dz(u, v, z):
    return -1

def dG_du(u, v, z):
    return 1

def dG_dv(u, v, z):
    return 1

def dG_dz(u, v, z):
    return 1

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
    
    # Evitando divisão por zero
    epsilon = 1e-10
    if abs(det2) < epsilon:
        det2 = epsilon
    if abs(det3) < epsilon:
        det3 = epsilon
    
    du_ds = det1 / det2
    dv_ds = -det2 / det3
    dz_ds = 1
    
    return [du_ds, dv_ds, dz_ds]

# Condições iniciais
u0 = 1  # Exemplo
v0 = 1  # Exemplo
z0 = 0  # Exemplo
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
