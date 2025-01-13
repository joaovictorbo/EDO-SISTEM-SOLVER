import numpy as np
from scipy.integrate import solve_ivp
import newton_funcoes as funcoes
from jacobiana import calcular_jacobiana

# Parâmetros globais
alpha = 0.001
epsilon_1 = 0.1
epsilon_2 = 0.01
muw0 = 1.0  # Viscosidade inicial sem polímero

# Funções dadas
def muw(z):  # Viscosidade da água
    return muw0 * 2**z

def muwz(z):  # dmuw/dz
    return np.log(2) * muw0 * 2**z

def muo():
    return 4.0

def mug():
    return 0.25

def D(u, v, z):
    w = 1 - u - v
    return u**2 / muw(z) + v**2 / muo() + w**2 / mug()

def f(u, v, z):
    return (u**2 / muw(z)) / D(u, v, z)

def g(u, v, z):
    return (v**2 / muo()) / D(u, v, z)

def df_du(u, v, z):
    return (2 * u / muw(z) * D(u, v, z) - u**2 / muw(z) * D_du(u, v, z)) / (D(u, v, z)**2)

def df_dv(u, v, z):
    return (-u**2 / muw(z) * D_dv(u, v, z)) / (D(u, v, z)**2)

def df_dz(u, v, z):
    return u**2 * (-muwz(z) / muw(z)**2 * D(u, v, z) - D_dz(u, v, z) / muw(z)) / (D(u, v, z)**2)

def dg_du(u, v, z):
    return (-v**2 / muo() * D_du(u, v, z)) / (D(u, v, z)**2)

def dg_dv(u, v, z):
    return (2 * v / muo() * D(u, v, z) - v**2 / muo() * D_dv(u, v, z)) / (D(u, v, z)**2)

def dg_dz(u, v, z):
    return -v**2 / muo() * D_dz(u, v, z) / (D(u, v, z)**2)

def D_du(u, v, z):
    w = 1 - u - v
    return 2 * u / muw(z) - 2 * w / mug()

def D_dv(u, v, z):
    w = 1 - u - v
    return 2 * v / muo() - 2 * w / mug()

def D_dz(u, v, z):
    return -u**2 * muwz(z) / muw(z)**2

def a(z):
    return np.sin(z)

def da_dz(z):
    return np.cos(z)

def lambdaz(u, v, z):
    return f(u, v, z) / (u + alpha * da_dz(z))

def detalpha(u, v, z):
    return ((df_du(u, v, z) - lambdaz(u, v, z)) * (dg_dv(u, v, z) - lambdaz(u, v, z)) - df_dv(u, v, z) * dg_du(u, v, z))

def det1(u, v, z):
    return (df_dv(u, v, z) * dg_dz(u, v, z) - df_dz(u, v, z) * (dg_dv(u, v, z) - lambdaz(u, v, z)))

def det2(u, v, z):
    return (df_dz(u, v, z) * dg_du(u, v, z) - (df_du(u, v, z) - lambdaz(u, v, z)) * dg_dz(u, v, z))

def p(u, v, z):
    return det1(u, v, z)

def q(u, v, z):
    return det2(u, v, z)

def r(u, v, z):
    return detalpha(u, v, z)

def h_L(z_R, z_L):
    if z_R != z_L:
        return (a(z_R) - a(z_L)) / (z_R - z_L)
    else:
        return da_dz(z_L)

# Função para calcular sigma_alpha
def sigma_alpha(u_L, f_R, z_R, z_L):
    h_L_value = h_L(z_R, z_L)
    denominator = u_L + alpha * h_L_value
    if denominator == 0:
        raise ValueError("Divisão por zero ao calcular sigma_alpha.")
    return f_R / denominator

# Sistema (50)
def system(s, y, u_R, v_R, z_R, f_R, g_R, sigma):
    u, v, z = y

    # Componentes do sistema (50)
    du_ds = f(u, v, z) - f_R - sigma * (u - u_R)
    dv_ds = g(u, v, z) - g_R - sigma * (v - v_R)
    dz_ds = (epsilon_1 / epsilon_2) * (
        f_R * (z - z_R) - sigma * (u_R * (z - z_R) - alpha * (a(z) - a(z_R)))
    )
    return [du_ds, dv_ds, dz_ds]

# Resolver as trajetórias do sistema
def resolver_trajetoria(u0, v0, z0, u_R, v_R, z_R, f_R, g_R, sigma, t_span=(0, 10), num_steps=500):
    def rhs(s, y):
        return system(s, y, u_R, v_R, z_R, f_R, g_R, sigma)

    t_values = np.linspace(t_span[0], t_span[1], num_steps)
    sol = solve_ivp(rhs, t_span, [u0, v0, z0], t_eval=t_values, method="RK45")

    return sol.t, sol.y

def calcular_autovalores_e_autovetores(jacobiana):
    """
    Calcula os autovalores e autovetores de uma matriz Jacobiana.
    Retorna os autovalores e autovetores positivos e negativos separadamente.
    """
    autovalores, autovetores = np.linalg.eig(jacobiana)
    indices_positivos = np.where(autovalores > 0)[0]
    indices_negativos = np.where(autovalores < 0)[0]
    return (autovalores[indices_positivos], autovetores[:, indices_positivos]), \
           (autovalores[indices_negativos], autovetores[:, indices_negativos])

if __name__ == "__main__":
    # Valores iniciais
    u_L, v_L, z_L = 0.1, 0.1, 0.1

    # Print initial values
    print(f"Valores iniciais: u_L = {u_L}, v_L = {v_L}, z_L = {z_L}")

    # Resolver trajetórias
    t_values, y_values = funcoes.resolver_trajetoria(u_L, v_L, z_L)
    u_R, v_R, z_R = t_values.y[0, 5], t_values.y[1, 5], t_values.y[2, 5]  # Últimos valores das trajetórias

    # Print final values
    print(f"Valores finais: u_R = {u_R}, v_R = {v_R}, z_R = {z_R}")
    f_L = f(u_L, v_L, z_L)
    g_L = g(u_L, v_L, z_L)
    f_R = f(u_R, v_R, z_R)
    g_R = g(u_R, v_R, z_R)
    sigmaLR = sigma_alpha(u_L, f_R, z_R, z_L)
    print(f"sigmaLR = {sigmaLR}")
    print(f"F/u{f_R/u_R}")
    # Calcular jacobianas
    jacoL = calcular_jacobiana(u_L, v_L, z_L, f_L, g_L, sigmaLR, epsilon_1, epsilon_2, alpha)
    jacoR = calcular_jacobiana(u_R, v_R, z_R, f_R, g_R, sigmaLR, epsilon_1, epsilon_2, alpha)

    # Print Jacobians in a visually appealing way
    np.set_printoptions(precision=3, suppress=True)
    print("Jacobiana em U^L:")
    print(jacoL)
    print("\nJacobiana em U^R:")
    print(jacoR)
    # Determinar autovalores e autovetores
    (autovalores_pos_L, autovetores_pos_L), (autovalores_neg_L, autovetores_neg_L) = calcular_autovalores_e_autovetores(jacoL)
    (autovalores_pos_R, autovetores_pos_R), (autovalores_neg_R, autovetores_neg_R) = calcular_autovalores_e_autovetores(jacoR)

    # Exibir autovalores e autovetores (já normalizados)
    print("Autovalores positivos de J(U^L):", autovalores_pos_L)
    print("Autovetores positivos de J(U^L):\n", autovetores_pos_L)

    print("\nAutovalores negativos de J(U^R):", autovalores_neg_R)
    print("Autovetores negativos de J(U^R):\n", autovetores_neg_R)
    
    import matplotlib.pyplot as plt

    # Simular as órbitas
    # tem que fazer sair todas as variacoes de U0 na mesma foto tanto somando quanto subtraindo os autovetores
    for i in range(autovetores_pos_L.shape[1]):
        U0_pos = u_L + 0.0001 * autovetores_pos_L[:, i]
        U0_neg = u_L - 0.0001 * autovetores_neg_L[:, i]

        print(f"\nSimulando órbita para U_k^0 positivo[{i}]:", U0_pos)
        t_values_pos, y_values_pos = resolver_trajetoria(U0_pos[0], U0_pos[1], U0_pos[2], u_R, v_R, z_R, f_R, g_R, sigmaLR,t_span=(0, 1000), num_steps=5000)
        print("Trajetória simulada positiva:", y_values_pos)
        min_diff_pos = np.min(np.linalg.norm(y_values_pos.T - [u_R, v_R, z_R], axis=1))
        print("Minima diferença em relação a U^R (positivo):", min_diff_pos)

        print(f"\nSimulando órbita para U_k^0 negativo[{i}]:", U0_pos)
        t_values_neg, y_values_neg = resolver_trajetoria(U0_pos[0], U0_pos[1], U0_pos[2], u_R, v_R, z_R, f_R, g_R, sigmaLR,t_span=(0, -1000), num_steps=5000)
        print("Trajetória simulada negativa:", y_values_neg)
        min_diff_neg = np.min(np.linalg.norm(y_values_neg.T - [u_R, v_R, z_R], axis=1))
        print("Minima diferença em relação a U^R (negativo):", min_diff_neg)

        # Plotar u por z e v por z
        plt.figure(figsize=(12, 6))

        plt.subplot(1, 3, 1)
        plt.plot(y_values_pos[0], y_values_pos[2], label=f'Órbita Positiva {i}')
        plt.plot(y_values_neg[0], y_values_neg[2], label=f'Órbita Negativa {i}')
        plt.plot(t_values.y[0], t_values.y[2], '--', label='Órbita U^L')
        plt.scatter([u_L], [z_L], color='red', label='Ponto Inicial U_L')
        plt.scatter([u_R], [z_R], color='blue', label='Ponto Inicial U_R')
        plt.scatter([0], [0], color='black')
        plt.scatter([1], [1], color='black')
        plt.xlabel('u')
        plt.ylabel('z')
        plt.title('u por z')
        plt.legend()

        plt.subplot(1, 3, 2)
        plt.plot(y_values_pos[1], y_values_pos[2], label=f'Órbita Positiva {i}')
        plt.plot(y_values_neg[1], y_values_neg[2], label=f'Órbita Negativa {i}')
        plt.plot(t_values.y[1],t_values.y[2], '--',label='Órbita U^L')
        plt.scatter([v_L], [z_L], color='red', label='Ponto Inicial V_L')
        plt.scatter([v_R], [z_R], color='blue', label='Ponto Inicial V_R')
        plt.scatter([0], [0], color='black')
        plt.scatter([1], [1], color='black')
        plt.xlabel('v')
        plt.ylabel('z')
        plt.title('v por z')
        plt.legend()

        plt.subplot(1, 3, 3)
        plt.plot(y_values_pos[0], y_values_pos[1], label=f'Órbita Positiva {i}')
        plt.plot(y_values_neg[0], y_values_neg[1], label=f'Órbita Negativa {i}')
        plt.plot(t_values.y[0],t_values.y[1], '--',label='Órbita U^L')
        plt.scatter([u_L], [v_L], color='red', label='Ponto Inicial V_L')
        plt.scatter([u_R], [v_R], color='blue', label='Ponto Inicial V_R')
        plt.scatter([0], [0], color='black')
        plt.scatter([1], [1], color='black')
        plt.xlabel('u')
        plt.ylabel('v')
        plt.title('u por v')
        plt.legend()

        plt.tight_layout()
        plt.show()
