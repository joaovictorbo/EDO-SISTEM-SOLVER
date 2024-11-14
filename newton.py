import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from sistemamodular.system import system

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

def df_du(u, v, c):
    D_val = D(u, v, c)
    D_du_val = D_du(u, v, c)
    numerator = (2 * u / muw(c)) * D_val - (u**2 / muw(c)) * D_du_val
    denominator = D_val**2
    return numerator / denominator

def df_dv(u, v, c):
    D_val = D(u, v, c)
    D_dv_val = D_dv(u, v, c)
    numerator = -(u**2 / muw(c)) * D_dv_val
    denominator = D_val**2
    return numerator / denominator

def df_dc(u, v, c):
    D_val = D(u, v, c)
    D_dc_val = D_dc(u, v, c)
    numerator = u**2 * (-muwc(c) / muw(c)**2 * D_val - D_dc_val / muw(c))
    denominator = D_val**2
    return numerator / denominator

def dg_du(u, v, c):
    D_val = D(u, v, c)
    D_du_val = D_du(u, v, c)
    numerator = -(v**2 / muo()) * D_du_val
    denominator = D_val**2
    return numerator / denominator

def dg_dv(u, v, c):
    D_val = D(u, v, c)
    D_dv_val = D_dv(u, v, c)
    numerator = (2 * v / muo()) * D_val - (v**2 / muo()) * D_dv_val
    denominator = D_val**2
    return numerator / denominator

def dg_dc(u, v, c):
    D_val = D(u, v, c)
    D_dc_val = D_dc(u, v, c)
    numerator = -(v**2 / muo()) * D_dc_val
    denominator = D_val**2
    return numerator / denominator

def D_du(u, v, c):
    w = 1 - u - v
    return 2 * u / muw(c) - 2 * w / mug()

def D_dv(u, v, c):
    w = 1 - u - v
    return 2 * v / muo() - 2 * w / mug()

def D_dc(u, v, c):
    return -u**2 * muwc(c) / muw(c)**2

def a(c):
    return np.sin(c)

def da_dc(c):
    return np.cos(c)  

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

def q(u, v, c):
    return det2(u, v, c)

def r(u, v, c):
    return detalpha(u, v, c)

def N(u, v, c):
    return np.sqrt(p(u, v, c)**2 + q(u, v, c)**2 + r(u, v, c)**2)
# Condições iniciais
u0 = 0.3
v0 = 0.6
z0 = 0.4
f0 = f(u0, v0, z0)
g0 = g(u0, v0, z0)
alpha = 10**-3

def F(u, v, z):
    f_value = f(u, v, z)
    g_value = g(u, v, z)
    return (f_value - f0) * (v - v0) - (g_value - g0) * (u - u0)

def G(u, v, z):
    f_value = f(u, v, z)
    a_z = a(z)
    a_z0 = a(z0)
    return (f_value - f0) * (u0 * (z - z0) + alpha * (a_z - a_z0)) - f0 * (z - z0) * (u - u0)

# Cálculo numérico da Jacobiana ∂(F,G)/∂(u,v)
def compute_jacobian(u, v, c, h=1e-6):
    # Derivadas parciais de F e G calcular pelo sistema 17 (sistema17.py)
    F_u = (F(u + h, v, c) - F(u - h, v, c)) / (2*h)
    F_v = (F(u, v + h, c) - F(u, v - h, c)) / (2*h)
    G_u = (G(u + h, v, c) - G(u - h, v, c)) / (2*h)
    G_v = (G(u, v + h, c) - G(u, v - h, c)) / (2*h)
    J = np.array([[F_u, F_v],
                  [G_u, G_v]])
    return J

# Método de Newton para correção de (u, v) com c congelado
def newton_correction(u0, v0, c_prev, iterations=3):
    u = u0
    v = v0
    for i in range(iterations):
        B = np.array([F(u, v, c_prev), G(u, v, c_prev)])
        J = compute_jacobian(u, v, c_prev)
        try:
            S = np.linalg.solve(J, -B)
        except np.linalg.LinAlgError:
            print("Jacobiana singular na iteração", i)
            break
        print(S)
        u += S[0]
        v += S[1]
    return u, v

# Função para verificar se um ponto está dentro do triângulo
def dentro_do_triangulo(u, v, c):
    return u >= 0 and v >= 0 and u + v <= 1 and 0 <= c <= 1

# Função para dividir as trajetórias ao entrar e sair do triângulo
def dividir_trajetorias(trajetoria):
    trajetorias = []
    traj_atual = []
    dentro = False

    for point in trajetoria:
        u, v, c = point
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

# Parâmetros iniciais
y0 = [u0, v0, z0]

# Passo pequeno
step_size = 0.01
# Parametro da curva de Hugoniot (s) como u(s) v(s) c(s)
s_values = np.arange(0, step_size*1000, step_size)

# Listas para armazenar as trajetórias
trajectory = []
trajectory.append([u0, v0, z0])

# Loop para gerar a curva de Hugoniot usando o método de Newton
for s in s_values[1:]:
    # Últimos valores de u, v, c
    y_0 = trajectory[-1]

    t_span = (s, s + step_size)
    sol = solve_ivp(system, t_span, y_0, method='LSODA', t_eval=np.linspace(s,s + step_size,2))
    u_prev, v_prev, c_prev = sol.y[:, -1]

    # Aplica o método de Newton para corrigir u e v
    u_new, v_new = newton_correction(u_prev, v_prev, c_prev)

    # Armazena os novos valores
    # c_new = c_prev, c_prev é atualizado no sistema 12
    trajectory.append([u_new, v_new, c_prev])

# Converte a trajetória em array numpy
trajectory = np.array(trajectory)

# Dividir as trajetórias
trajetorias_divididas = dividir_trajetorias(trajectory)

# Plotagem das trajetórias divididas
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for traj in trajetorias_divididas:
    traj = np.array(traj)
    ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], label='Trajetória')

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
ax.scatter(u0, v0, z0, color='red', s=50, label='Ponto Inicial', edgecolor='black')
# Passo pequeno para trajetória reversa
s_values_rev = np.arange(0, -step_size*1000, -step_size)

# Listas para armazenar as trajetórias reversas
trajectory_rev = []
trajectory_rev.append([u0, v0, z0])

# Loop para gerar a curva de Hugoniot reversa
for s in s_values_rev[1:]:
    # Últimos valores de u, v, c
    y_0 = trajectory_rev[-1]
    t_span = (s, s - step_size)
    sol = solve_ivp(system, t_span, y_0, method='LSODA', t_eval=np.linspace(s,s - step_size,2))
    u_prev, v_prev, c_prev = sol.y[:, -1]

    # Aplica o método de Newton para corrigir u e v
    u_new, v_new = newton_correction(u_prev, v_prev, c_prev)


    # Armazena os novos valores
    # c_new = c_prev, c_prev é atualizado no sistema 12
    trajectory_rev.append([u_new, v_new, c_prev])

# Converte a trajetória reversa em array numpy
trajectory_rev = np.array(trajectory_rev)

# Dividir as trajetórias reversas
trajetorias_divididas_rev = dividir_trajetorias(trajectory_rev)

# Plotar as trajetórias reversas
for traj in trajetorias_divididas_rev:
    traj = np.array(traj)
    ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], label='Trajetória Reversa')

plt.show()
