import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

alpha = 10**-3

muw0 = 1.0  # Viscosidade inicial sem polimero

# Funções dadas
def muw(c):  # Viscosidade da agua
    return muw0 * 2**c

def muwc(c):  # dmuw/dc
    return np.log(2) * muw0 * 2**c

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

def a(c):
    return np.sqrt(c)

def da_dc(c):
    return 1.0 / (2*np.sqrt(c))    

def lambdac(u,v,c):
    return f(u, v, c)/ (u + alpha* da_dc(c))

def detalpha(u, v, c):
    return ((df_du(u, v, c)-lambdac(u,v,c)) * (dg_dv(u, v, c)-lambdac(u,v,c)) - df_dv(u,v,c) * dg_du(u, v, c))

def det1(u, v, c):
    return (df_dv(u, v, c) * dg_dc(u, v, c) - df_dc(u,v,c) * (dg_dv(u, v, c)-lambdac(u,v,c)))

def det2(u,v,c):
    return (df_dc(u,v,c) * dg_du(u,v,c) - (df_du(u,v,c)-lambdac(u,v,c)) * dg_dc(u,v,c))

def p(u, v, c):
    return det1(u, v, c)
#remover os if e deixar somente o if de ponto singular e a formula 16
def q(u, v, c):
    return det2(u, v, c)

def r(u, v, c):
    return detalpha(u, v, c)

def N(u, v, c):
    return np.sqrt(p(u, v, c)**2 + q(u, v, c)**2 + r(u, v, c)**2)
# Sistema de EDOs
def system(s, y):
    u, v, c = y
    eq1 = p(u, v, c) / N(u, v, c) if N(u, v, c) else None
    eq2 = q(u, v, c) / N(u, v, c) if N(u, v, c) else None
    eq3 = r(u,v,c) / N(u, v, c) if N(u, v, c) else None
    
    if eq1 is None or eq2 is None or eq3 is None:
        raise ValueError("Ponto singular do sistema de EDOs")
    # caso seja 0 retorna 0.01 A fazer
    return [eq1, eq2, eq3]

# Condições iniciais
u0, v0, c0 = 0.45, 0.47, 0.2
y0 = [u0, v0, c0]

# Intervalo de integração
t_span = (0,10)

# Resolução do sistema
sol = solve_ivp(system, t_span, y0, method='LSODA', t_eval=np.linspace(0,1,2000))
t_span2 = (0,-10)

# Resolução do sistema
sol2 = solve_ivp(system, t_span2, y0, method='LSODA', t_eval2=np.linspace(0,-1,2000))
# Verificação dos resultados
if sol.y.shape[0] != 3:
    raise ValueError("A solução não retornou os valores esperados.")

# Plot dos resultados em 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(sol.y[0], sol.y[1], sol.y[2], label='Trajetória')
ax.plot(sol2.y[0], sol2.y[1], sol2.y[2], label='Trajetória2')
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

plt.title('Solução do sistema de EDOs em 3D')
plt.show()
