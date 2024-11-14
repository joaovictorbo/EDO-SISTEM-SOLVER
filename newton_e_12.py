import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import sistemamodular.system as system

def muo():
    return 4.0

def mug():
    return 0.25


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

def F(u, v, z):
    f_value = f(u, v, z)
    g_value = g(u, v, z)
    return (f_value - f0) * (v - v0) - (g_value - g0) * (u - u0)

def G(u, v, z):
    f_value = f(u, v, z)
    a_z = a(z)
    a_z0 = a(z0)
    return (f_value - f0) * (u0 * (z - z0) + alpha * (a_z - a_z0)) - f0 * (z - z0) * (u - u0)





alpha = 10**-3
muw0 = 1.0  # Viscosidade inicial sem polimero
u0 = 0.1
v0 = 0.6
z0 = 0.2
alpha = 10**-3
def f(u, v, c):
    return (u**2 / muw(c)) / D(u, v, c)

def g(u, v, c):
    return (v**2 / muo()) / D(u, v, c)

f0 = f(u0, v0, z0)
g0 = g(u0, v0, z0)


# Cálculo numérico da Jacobiana ∂(F,G)/∂(u,v)
def compute_jacobian(u, v, c, h=1e-6):
    F_u = (F(u + h, v, c) - F(u - h, v, c)) / (2*h)
    F_v = (F(u, v + h, c) - F(u, v - h, c)) / (2*h)
    G_u = (G(u + h, v, c) - G(u - h, v, c)) / (2*h)
    G_v = (G(u, v + h, c) - G(u, v - h, c)) / (2*h)
    J = np.array([[F_u, F_v],
                  [G_u, G_v]])
    return J

def newton_correction(u0, v0, c_fixed, iterations=3):
    u = u0
    v = v0
    for i in range(iterations):
        B = np.array([F(u, v, c_fixed), G(u, v, c_fixed)])
        J = compute_jacobian(u, v, c_fixed)
        try:
            S = np.linalg.solve(J, -B)
        except np.linalg.LinAlgError:
            print("Jacobiana singular na iteração", i)
            break
        u += S[0]
        v += S[1]
    return u, v

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

def dentro_do_triangulo(u, v, c):
    return u >= 0 and v >= 0 and u + v <= 1 and 0 <= c <= 1

# Inicialização para o primeiro sistema
u0, v0, c0 = 0.1, 0.6, 0.2
step_size = 0.01
s_values = np.arange(0, 10, step_size)

trajectory = [[u0, v0, c0]]
for s in s_values[1:]:
    u_prev, v_prev, c_prev = trajectory[-1]
    c_fixed = c_prev
    u_new, v_new = newton_correction(u_prev, v_prev, c_fixed)

    N_val = N(u_new, v_new, c_fixed)
    if N_val == 0:
        break
    du_ds = p(u_new, v_new, c_fixed) / N_val
    dv_ds = q(u_new, v_new, c_fixed) / N_val
    dc_ds = r(u_new, v_new, c_fixed) / N_val

    c_new = c_fixed + dc_ds * step_size
    trajectory.append([u_new, v_new, c_new])

trajectory = np.array(trajectory)
trajetorias_divididas = dividir_trajetorias(trajectory)

# Inicialização para o segundo sistema
u0, v0, c0 = 0.1, 0.6, 0.3
y0 = [u0, v0, c0]
t_span = (0, 10)
sol = solve_ivp(system, t_span, y0, method='LSODA', t_eval=np.linspace(0, 10, 2000))
t_span2 = (0, -10)
sol2 = solve_ivp(system, t_span2, y0, method='LSODA', t_eval=np.linspace(0, -10, 2000))

trajetorias1 = dividir_trajetorias(np.column_stack((sol.y[0], sol.y[1], sol.y[2])))
trajetorias2 = dividir_trajetorias(np.column_stack((sol2.y[0], sol2.y[1], sol2.y[2])))

# Plotagem das trajetórias
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Primeiro sistema (trajetória com método de Newton)
for traj in trajetorias_divididas:
    traj = np.array(traj)
    ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], color='blue', label='Trajetória Sistema 1')

# Segundo sistema (ODE solver)
for traj in trajetorias1:
    traj = np.array(traj)
    ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], color='green', label='Trajetória Sistema 2')

for traj in trajetorias2:
    traj = np.array(traj)
    ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], color='red', label='Trajetória Reversa Sistema 2')

ax.set_xlabel('u(s)')
ax.set_ylabel('v(s)')
ax.set_zlabel('c(s)')
ax.legend()

# Adicionar triângulo
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
ax.scatter(u0, v0, c0, color='orange', s=50, label='Ponto Inicial', edgecolor='black')

plt.show()
