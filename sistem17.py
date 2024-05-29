import numpy as np
from scipy.integrate import solve_ivp

# Funções de fluxo
def f(u, v, z):
    lambda_w = (u**2) / (u**2 + v**2 + (1 - u - v)**2)
    return lambda_w / (lambda_w + lambda_o(v) + lambda_g(u, v))

def g(u, v, z):
    lambda_o_val = lambda_o(v)
    return lambda_o_val / (lambda_w(u, z) + lambda_o_val + lambda_g(u, v))

def lambda_o(v):
    return (v**2) / (v**2 + (1 - v)**2)

def lambda_w(u, z):
    return (u**2) / (u**2 + z**2)

def lambda_g(u, v):
    sg = 1 - u - v
    return (sg**2) / (sg**2 + (1 - sg)**2)

# Função de adsorção
def a(z):
    return np.sqrt(z)

# Campo de vetores X_alpha(U)
def X_alpha(U, alpha, epsilon_ratio, f_R_alpha, sigma_alpha):
    u, v, z = U
    du_ds = f(u, v, z) - f_R_alpha - sigma_alpha * (u - u_R_alpha)
    dv_ds = g(u, v, z) - g_R_alpha - sigma_alpha * (v - v_R_alpha)
    dz_ds = epsilon_ratio * (f_R_alpha * (z - z_R_alpha) - sigma_alpha * (u_R_alpha * (z - z_R_alpha) - alpha * (a(z) - a(z_R_alpha))))
    return [du_ds, dv_ds, dz_ds]

# Função para resolver o sistema de EDOs
def solve_system(U_L_alpha, U_R_alpha, alpha, epsilon_ratio, t_span, t_eval):
    # Valores iniciais
    u_L_alpha, v_L_alpha, z_L_alpha = U_L_alpha
    u_R_alpha, v_R_alpha, z_R_alpha = U_R_alpha
    sigma_alpha = f(U_L_alpha[0], U_L_alpha[1], U_L_alpha[2]) / (u_L_alpha + alpha * h(a, z_R_alpha, z_L_alpha))

    def system(s, U):
        return X_alpha(U, alpha, epsilon_ratio, f_R_alpha, sigma_alpha)

    # Resolver o sistema de EDOs
    sol = solve_ivp(system, t_span, U_L_alpha, t_eval=t_eval, method='RK45')
    return sol

# Função para h(z)
def h(a, z, z0):
    if z != z0:
        return (a(z) - a(z0)) / (z - z0)
    else:
        return a(z)

# Parâmetros
alpha = 1e-3
epsilon_ratio = 1.0
t_span = (0, 10)
t_eval = np.linspace(0, 10, 100)

# Condições iniciais
U_L_alpha = [0.2, 0.3, 0.1]  # Exemplo de valores iniciais
U_R_alpha = [0.1, 0.4, 0.2]  # Exemplo de valores finais

# Resolver o sistema
solution = solve_system(U_L_alpha, U_R_alpha, alpha, epsilon_ratio, t_span, t_eval)

# Plotar os resultados
import matplotlib.pyplot as plt

plt.plot(solution.t, solution.y[0], label='u(s)')
plt.plot(solution.t, solution.y[1], label='v(s)')
plt.plot(solution.t, solution.y[2], label='z(s)')
plt.xlabel('s')
plt.ylabel('Valores')
plt.legend()
plt.show()
