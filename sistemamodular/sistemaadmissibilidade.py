import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Import necessário para plot 3D

import newton_funcoes as funcoes 
from jacobiana import calcular_jacobiana

# ----- PARÂMETROS GLOBAIS -----
alpha = 0.001
epsilon_1 = 0.1
epsilon_2 = 0.01
muw0 = 1.0  # Viscosidade inicial sem polímero

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

def D_du(u, v, z):
    w = 1 - u - v
    return 2 * u / muw(z) - 2 * w / mug()

def D_dv(u, v, z):
    w = 1 - u - v
    return 2 * v / muo() - 2 * w / mug()

def D_dz(u, v, z):
    return -u**2 * muwz(z) / muw(z)**2

def df_du(u, v, z):
    return (2 * u / muw(z) * D(u, v, z) - (u**2 / muw(z))*D_du(u, v, z)) / (D(u, v, z)**2)

def df_dv(u, v, z):
    return (-(u**2 / muw(z)) * D_dv(u, v, z)) / (D(u, v, z)**2)

def df_dz(u, v, z):
    return u**2 * (
        -muwz(z) / muw(z)**2 * D(u, v, z) - D_dz(u, v, z) / muw(z)
    ) / (D(u, v, z)**2)

def dg_du(u, v, z):
    return (-(v**2 / muo()) * D_du(u, v, z)) / (D(u, v, z)**2)

def dg_dv(u, v, z):
    return (2 * v / muo() * D(u, v, z) - (v**2 / muo())*D_dv(u, v, z)) / (D(u, v, z)**2)

def dg_dz(u, v, z):
    return -(v**2 / muo()) * D_dz(u, v, z) / (D(u, v, z)**2)

def a(z):
    return np.sin(z)

def da_dz(z):
    return np.cos(z)

def lambdaz(u, v, z):
    return f(u, v, z) / (u + alpha * da_dz(z))

def detalpha(u, v, z):
    return ((df_du(u, v, z) - lambdaz(u, v, z))
            * (dg_dv(u, v, z) - lambdaz(u, v, z))
            - df_dv(u, v, z)*dg_du(u, v, z))

def det1(u, v, z):
    return (df_dv(u, v, z)*dg_dz(u, v, z)
            - df_dz(u, v, z)*(dg_dv(u, v, z) - lambdaz(u, v, z)))

def det2(u, v, z):
    return (df_dz(u, v, z)*dg_du(u, v, z)
            - (df_du(u, v, z) - lambdaz(u, v, z))*dg_dz(u, v, z))

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

def sigma_alpha(u_L, f_R, z_R, z_L):
    h_L_value = h_L(z_R, z_L)
    denominator = u_L + alpha * h_L_value
    if denominator == 0:
        raise ValueError("Divisão por zero ao calcular sigma_alpha.")
    return f_R / denominator

# Sistema (50)
def system(s, y, u_R, v_R, z_R, f_R, g_R, sigma):
    u, v, z = y
    du_ds = f(u, v, z) - f_R - sigma * (u - u_R)
    dv_ds = g(u, v, z) - g_R - sigma * (v - v_R)
    dz_ds = (epsilon_1 / epsilon_2) * (
        f_R * (z - z_R) - sigma * (u_R * (z - z_R) - alpha * (a(z) - a(z_R)))
    )
    return [du_ds, dv_ds, dz_ds]

def resolver_trajetoria(u0, v0, z0, u_R, v_R, z_R, f_R, g_R, sigma,
                        t_span=(0,10), num_steps=500):
    """
    Resolve o IVP para os valores iniciais (u0, v0, z0).
    """
    def rhs(s, y):
        return system(s, y, u_R, v_R, z_R, f_R, g_R, sigma)

    t_eval = np.linspace(t_span[0], t_span[1], num_steps)
    sol = solve_ivp(rhs, t_span, [u0, v0, z0], t_eval=t_eval, method='RK45')
    return sol.t, sol.y

def calcular_autovalores_e_autovetores(jacobiana):
    """
    Calcula os autovalores e autovetores de uma matriz Jacobiana.
    Retorna os autovalores e autovetores positivos e negativos separadamente.
    """
    autovalores, autovetores = np.linalg.eig(jacobiana)
    indices_positivos = np.where(autovalores > 0)[0]
    indices_negativos = np.where(autovalores < 0)[0]
    return ((autovalores[indices_positivos], autovetores[:, indices_positivos]),
            (autovalores[indices_negativos], autovetores[:, indices_negativos]))

def gerar_pontos_iniciais(u_c, v_c, z_c, N, N2, raio):
    """
    Gera N pontos iniciais em um círculo de raio 'raio' no plano (u,v)
    em torno de (u_c, v_c). Mantém z = z_c.
    arrumar pra ser 3d
    """
    thetas = np.linspace(0, 2*np.pi, N, endpoint=False)
    phis = np.linspace(0, np.pi, N2, endpoint=False)
    pontos = []
    for theta in thetas:
        for phi in phis:
            u_i = u_c + raio * np.cos(theta) * np.sin(phi)
            v_i = v_c + raio * np.sin(theta) * np.sin(phi)
            z_i = z_c + raio * np.cos(phi)
            pontos.append([u_i, v_i, z_i])
    return pontos

if __name__ == "__main__":
    # Valores iniciais
    u_L, v_L, z_L = 0.1, 0.1, 0.1
    print(f"Valores iniciais: u_L = {u_L}, v_L = {v_L}, z_L = {z_L}")

    t_vals, sol_prim = funcoes.resolver_trajetoria(u_L, v_L, z_L)

    # Pegamos o final da trajetória como ponto R
    u_R = sol_prim.y[0, 1]
    v_R = sol_prim.y[1, 1]
    z_R = sol_prim.y[2, 1]
    distance = np.sqrt((u_R - u_L)**2 + (v_R - v_L)**2 + (z_R - z_L)**2)
    print(f"Valores finais: u_R = {u_R}, v_R = {v_R}, z_R = {z_R}")
    
    f_L = f(u_L, v_L, z_L)
    g_L = g(u_L, v_L, z_L)
    f_R = f(u_R, v_R, z_R)
    g_R = g(u_R, v_R, z_R)

    # Cálculo de sigma
    sigmaLR = sigma_alpha(u_L, f_R, z_R, z_L)
    print(f"sigmaLR = {sigmaLR}")
    if u_R != 0:
        print(f"F/u: {f_R/u_R}")
    else:
        print("Divisão por zero em u_R para cálculo de F/u.")

    # Jacobianas
    jacoL = calcular_jacobiana(u_L, v_L, z_L, f_L, g_L, sigmaLR, epsilon_1, epsilon_2, alpha)
    jacoR = calcular_jacobiana(u_R, v_R, z_R, f_R, g_R, sigmaLR, epsilon_1, epsilon_2, alpha)

    np.set_printoptions(precision=3, suppress=True)
    print("Jacobiana em U^L:")
    print(jacoL)
    print("\nJacobiana em U^R:")
    print(jacoR)

    # Autovalores e autovetores
    (autovalores_pos_L, autovetores_pos_L), (autovalores_neg_L, autovetores_neg_L) = calcular_autovalores_e_autovetores(jacoL)
    (autovalores_pos_R, autovetores_pos_R), (autovalores_neg_R, autovetores_neg_R) = calcular_autovalores_e_autovetores(jacoR)

    print("Autovalores positivos de J(U^L):", autovalores_pos_L)
    print("Autovetores positivos de J(U^L):\n", autovetores_pos_L)

    print("\nAutovalores negativos de J(U^R):", autovalores_neg_R)
    print("Autovetores negativos de J(U^R):\n", autovetores_neg_R)

    # -------------- Exemplo de gerar N pontos iniciais --------------
    N = 6
    N2 = 4         # número de pontos
    raio = distance / 2
    pontos_iniciais = gerar_pontos_iniciais(u_L, v_L, z_L, N,N2, raio)

    # ======= PLOTS 2D =======
    fig, axes = plt.subplots(1, 3, figsize=(12, 5))

    # Converte a solução base para plotar via quiver (a cada "step" para não poluir)
    step = 5  
    base_u = sol_prim.y[0, ::step]
    base_v = sol_prim.y[1, ::step]
    base_z = sol_prim.y[2, ::step]

    # -- Subplot (u x z)
    ax_uz = axes[0]
    ax_uz.quiver(base_u[:-1], base_z[:-1],
                 base_u[1:] - base_u[:-1],
                 base_z[1:] - base_z[:-1],
                 angles='xy', scale_units='xy', scale=1,
                 color='gray', label='Trajet. base (U^L)')

    ax_uz.scatter([u_L], [z_L], color='red', label='Ponto L')
    ax_uz.scatter([u_R], [z_R], color='blue', label='Ponto R')
    ax_uz.set_xlabel('u')
    ax_uz.set_ylabel('z')
    ax_uz.set_title('u x z')
    ax_uz.legend()

    # -- Subplot (v x z)
    ax_vz = axes[1]
    ax_vz.quiver(base_v[:-1], base_z[:-1],
                 base_v[1:] - base_v[:-1],
                 base_z[1:] - base_z[:-1],
                 angles='xy', scale_units='xy', scale=1,
                 color='gray', label='Trajet. base (U^L)')

    ax_vz.scatter([v_L], [z_L], color='red', label='Ponto L')
    ax_vz.scatter([v_R], [z_R], color='blue', label='Ponto R')
    ax_vz.set_xlabel('v')
    ax_vz.set_ylabel('z')
    ax_vz.set_title('v x z')
    ax_vz.legend()

    # -- Subplot (u x v)
    ax_uv = axes[2]
    ax_uv.quiver(base_u[:-1], base_v[:-1],
                 base_u[1:] - base_u[:-1],
                 base_v[1:] - base_v[:-1],
                 angles='xy', scale_units='xy', scale=1,
                 color='gray', label='Trajet. base (U^L)')

    ax_uv.scatter([u_L], [v_L], color='red', label='Ponto L')
    ax_uv.scatter([u_R], [v_R], color='blue', label='Ponto R')
    ax_uv.set_xlabel('u')
    ax_uv.set_ylabel('v')
    ax_uv.set_title('u x v')
    ax_uv.legend()

    # ================================================================
    # Para cada ponto inicial do círculo, resolvemos a trajetória
    # e plotamos nos subplots 2D.
    # ================================================================
    for idx, (u0, v0, z0) in enumerate(pontos_iniciais):
        t_loc, sol_loc = resolver_trajetoria(
            u0, v0, z0, u_R, v_R, z_R, f_R, g_R, sigmaLR,
            t_span=(0, 100), num_steps=1000
        )

        stepQ = 20
        loc_u = sol_loc[0, ::stepQ]
        loc_v = sol_loc[1, ::stepQ]
        loc_z = sol_loc[2, ::stepQ]

        # Quiver (u x z)
        ax_uz.quiver(loc_u[:-1], loc_z[:-1],
                     loc_u[1:] - loc_u[:-1],
                     loc_z[1:] - loc_z[:-1],
                     angles='xy', scale_units='xy', scale=1,
                     label=f'Órbita {idx}' if idx < 2 else None)

        # Quiver (v x z)
        ax_vz.quiver(loc_v[:-1], loc_z[:-1],
                     loc_v[1:] - loc_v[:-1],
                     loc_z[1:] - loc_z[:-1],
                     angles='xy', scale_units='xy', scale=1)

        # Quiver (u x v)
        ax_uv.quiver(loc_u[:-1], loc_v[:-1],
                     loc_u[1:] - loc_u[:-1],
                     loc_v[1:] - loc_v[:-1],
                     angles='xy', scale_units='xy', scale=1)

    plt.tight_layout()
    plt.show()

    # ====== PLOT 3D ======
    # Agora criamos um plot em 3D para visualizar a trajetória no espaço (u, v, z).
    fig_3d = plt.figure(figsize=(8, 6))
    ax_3d = fig_3d.add_subplot(111, projection='3d')

    # Trajetória base (U^L -> U^R)
    ax_3d.plot(base_u, base_v, base_z, color='gray', label='Trajet. base (U^L)')

    # Marca pontos L e R
    ax_3d.scatter(u_L, v_L, z_L, color='red', s=50, label='Ponto L')
    ax_3d.scatter(u_R, v_R, z_R, color='blue', s=50, label='Ponto R')

    # Para cada ponto inicial do círculo, plotamos também a trajetória resultante em 3D:
    for idx, (u0, v0, z0) in enumerate(pontos_iniciais):
        t_loc, sol_loc = resolver_trajetoria(
            u0, v0, z0, u_R, v_R, z_R, f_R, g_R, sigmaLR,
            t_span=(0, 100), num_steps=1000
        )
        ax_3d.plot(sol_loc[0], sol_loc[1], sol_loc[2], alpha=0.8)

    ax_3d.set_xlabel('u')
    ax_3d.set_ylabel('v')
    ax_3d.set_zlabel('z')
    ax_3d.set_title('Trajetórias em 3D (u, v, z)')
    ax_3d.legend()
    plt.show()