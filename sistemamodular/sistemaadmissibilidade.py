import numpy as np
from scipy.integrate import solve_ivp
import newton_funcoes as funcoes
from jacobiana import calcular_jacobiana

# Parâmetros globais
alpha = 0.1
epsilon_1 = 1.0
epsilon_2 = 1.0
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
    epsilon_1, epsilon_2 = 1.0, 1.0

    # Resolver trajetórias
    t_values, y_values = funcoes.resolver_trajetoria(u_L, v_L, z_L)
    u_R, v_R, z_R = y_values.y[0, -1], y_values.y[1, -1], y_values.y[2, -1]  # Últimos valores das trajetórias

    f_R = f(u_R, v_R, z_R)
    g_R = g(u_R, v_R, z_R)
    sigmaLR = sigma_alpha(u_L, f_R, z_R, z_L)

    # Calcular jacobianas
    jacoL = calcular_jacobiana(u_L, v_L, z_L, f_R, g_R, epsilon_1, epsilon_2, alpha=0)
    jacoR = calcular_jacobiana(u_R, v_R, z_R, f_R, g_R, epsilon_1, epsilon_2, alpha=1)

    # Determinar autovalores e autovetores
    (autovalores_pos_L, autovetores_pos_L), (autovalores_neg_L, autovetores_neg_L) = calcular_autovalores_e_autovetores(jacoL)
    (autovalores_pos_R, autovetores_pos_R), (autovalores_neg_R, autovetores_neg_R) = calcular_autovalores_e_autovetores(jacoR)

    # Exibir autovalores e autovetores
    print("Autovalores positivos de J(U^L):", autovalores_pos_L)
    print("Autovetores positivos de J(U^L):\n", autovetores_pos_L)

    print("\nAutovalores negativos de J(U^R):", autovalores_neg_R)
    print("Autovetores negativos de J(U^R):\n", autovetores_neg_R)

    # Determinar os pontos iniciais U_k^0
    U_k0 = [u_L + autovetores_pos_L[:, i] for i in range(autovetores_pos_L.shape[1])]

    # Simular as órbitas
    sigma = 0.5
    for i, U0 in enumerate(U_k0):
        print(f"\nSimulando órbita para U_k^0[{i}]:", U0)
        t_values, y_values = resolver_trajetoria(U0[0], U0[1], U0[2], u_R, v_R, z_R, f_R, g_R, sigma)

        max_diff = np.max(np.abs(y_values.T - [u_R, v_R, z_R]), axis=0)
        print("Máxima diferença em relação a U^R:", max_diff)

        if np.all(max_diff < 1e-10):
            print(f"A órbita para U_k^0[{i}] aproxima U^R satisfatoriamente.")
        else:
            print(f"A órbita para U_k^0[{i}] NÃO aproxima U^R.")
