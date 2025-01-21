# newton_funcoes.py
import numpy as np
from scipy.integrate import solve_ivp
from system import system

alpha = 0.001  # Parâmetro de relaxamento
muw0 = 1.0  # Viscosidade inicial sem polimero

def muw(c):
    return muw0 * 2**c

def muwc(c):
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

#Pertencem ao sistema 30
def F(u, v, z, u0, v0, f0, g0):
    return (f(u, v, z) - f0) * (v - v0) - (g(u, v, z) - g0) * (u - u0)

def G(u, v, z, u0, f0, z0):
    return (f(u, v, z) - f0) * (u0 * (z - z0) + alpha * (np.sin(z) - np.sin(z0))) - f0 * (z - z0) * (u - u0)

def compute_jacobian(u, v, c, u0, v0, f0, g0, z0, h=1e-6):
#VER FORMULAS 31 a 36
    F_u = (F(u + h, v, c, u0, v0, f0, g0) - F(u - h, v, c, u0, v0, f0, g0)) / (2*h)
    F_v = (F(u, v + h, c, u0, v0, f0, g0) - F(u, v - h, c, u0, v0, f0, g0)) / (2*h)
    G_u = (G(u + h, v, c, u0, f0, z0) - G(u - h, v, c, u0, f0, z0)) / (2*h)
    G_v = (G(u, v + h, c, u0, f0, z0) - G(u, v - h, c, u0, f0, z0)) / (2*h)
    return np.array([[F_u, F_v], [G_u, G_v]])

def newton_correction(u0, v0, c_prev, f0, g0, z0, iterations=10):
    u, v = u0, v0
    for _ in range(iterations):
        B = np.array([F(u, v, c_prev, u0, v0, f0, g0), G(u, v, c_prev, u0, f0, z0)])
        J = compute_jacobian(u, v, c_prev, u0, v0, f0, g0, z0)
        try:
            S = np.linalg.solve(J, -B)
            u, v = u + S[0], v + S[1]
        except np.linalg.LinAlgError:
            print("Jacobiana singular. Interrompendo ajuste.")
            break
    return u, v

def resolver_trajetoria(u0, v0, c0, t_span=(0, 10), t_span2=(0, -10), num_steps=2000):
    """
    Resolve a trajetória a partir das condições iniciais fornecidas, aplicando a correção de Newton em cada passo.

    Parâmetros:
        u0, v0, c0: Condições iniciais para u, v, c.
        t_span, t_span2: Intervalos de tempo para integração positiva e negativa.
        num_steps: Número de passos para a integração.

    Retorna:
        sol, sol2: Soluções para a trajetória positiva e negativa.
    """
    # Integração no sentido positivo
    t_values = np.linspace(t_span[0], t_span[1], num_steps + 1)
    y_values = np.zeros((3, num_steps + 1))
    y_values[:, 0] = [u0, v0, c0]

    for i in range(num_steps):
        t_i = t_values[i]
        y_i = y_values[:, i]
        delta_t = t_values[i + 1] - t_i

        # Método de Runge-Kutta de quarta ordem (RK4)
        k1 = delta_t * system(t_i, y_i)
        k2 = delta_t * system(t_i + delta_t / 2, y_i + k1 / 2)
        k3 = delta_t * system(t_i + delta_t / 2, y_i + k2 / 2)
        k4 = delta_t * system(t_i + delta_t, y_i + k3)
        y_next = y_i + (k1 + 2 * k2 + 2 * k3 + k4) / 6

        # Preparar para a correção de Newton
        u0_prev, v0_prev, c0_prev = y_i
        f0 = f(u0_prev, v0_prev, c0_prev)
        g0 = g(u0_prev, v0_prev, c0_prev)
        z0 = c0_prev
        c_prev = y_next[2]

        # Aplicar a correção de Newton
        u_corr, v_corr = newton_correction(u0_prev, v0_prev, c_prev, f0, g0, z0)
        y_next[0], y_next[1] = u_corr, v_corr

        # Atualizar os valores para o próximo passo
        y_values[:, i + 1] = y_next

    # Integração no sentido negativo
    t_values_neg = np.linspace(t_span2[0], t_span2[1], num_steps + 1)
    y_values_neg = np.zeros((3, num_steps + 1))
    y_values_neg[:, 0] = [u0, v0, c0]

    for i in range(num_steps):
        t_i = t_values_neg[i]
        y_i = y_values_neg[:, i]
        delta_t = t_values_neg[i + 1] - t_i

        # Método de Runge-Kutta de quarta ordem (RK4)
        k1 = delta_t * system(t_i, y_i)
        k2 = delta_t * system(t_i + delta_t / 2, y_i + k1 / 2)
        k3 = delta_t * system(t_i + delta_t / 2, y_i + k2 / 2)
        k4 = delta_t * system(t_i + delta_t, y_i + k3)
        y_next = y_i + (k1 + 2 * k2 + 2 * k3 + k4) / 6

        # Preparar para a correção de Newton
        u0_prev, v0_prev, c0_prev = y_i
        f0 = f(u0_prev, v0_prev, c0_prev)
        g0 = g(u0_prev, v0_prev, c0_prev)
        z0 = c0_prev
        c_prev = y_next[2]

        # Aplicar a correção de Newton
        u_corr, v_corr = newton_correction(u0_prev, v0_prev, c_prev, f0, g0, z0)
        y_next[0], y_next[1] = u_corr, v_corr

        # Atualizar os valores para o próximo passo
        y_values_neg[:, i + 1] = y_next

    # Criar objetos de solução semelhantes ao retornado pelo solve_ivp
    class Solution:
        def __init__(self, t, y):
            self.t = t
            self.y = y

    sol = Solution(t_values, y_values)
    sol2 = Solution(t_values_neg, y_values_neg)

    return sol, sol2


# Função para verificar se um ponto está dentro do triângulo
def dentro_do_triangulo(u, v, c):
    return u >= 0 and v >= 0 and u + v <= 1 and 0 <= c <= 1

# Função para dividir a trajetória com base na condição de estar dentro do triângulo
def dividir_trajetorias(sol):
    """
    Divide a trajetória em partes, de acordo com a condição de estar dentro do triângulo.
    
    Parâmetros:
        sol: Solução obtida da integração.
    
    Retorna:
        Lista de trajetórias divididas.
    """
    trajetorias = []
    traj_atual = []
    dentro = False
    
    for i in range(len(sol.y[0])):
        u, v, c = sol.y[0][i], sol.y[1][i], sol.y[2][i]
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
