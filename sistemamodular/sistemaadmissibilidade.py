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
    return (2 * u / muw(z) * D(u, v, z) - (u**2 / muw(z)) * D_du(u, v, z)) / (D(u, v, z)**2)

def df_dv(u, v, z):
    return (-(u**2 / muw(z)) * D_dv(u, v, z)) / (D(u, v, z)**2)

def df_dz(u, v, z):
    return u**2 * (
        -muwz(z) / muw(z)**2 * D(u, v, z) - D_dz(u, v, z) / muw(z)
    ) / (D(u, v, z)**2)

def dg_du(u, v, z):
    return (-(v**2 / muo()) * D_du(u, v, z)) / (D(u, v, z)**2)

def dg_dv(u, v, z):
    return (2 * v / muo() * D(u, v, z) - (v**2 / muo()) * D_dv(u, v, z)) / (D(u, v, z)**2)

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
    Retorna:
      - sol_neg: solução para s de 0 até -t_span[1]
      - sol_pos: solução para s de 0 até t_span[1]
    """
    def rhs(s, y):
        return system(s, y, u_R, v_R, z_R, f_R, g_R, sigma)

    # Tempo "positivo"
    t_eval_pos = np.linspace(t_span[0], t_span[1], num_steps)
    sol_pos = solve_ivp(rhs, t_span, [u0, v0, z0],
                        t_eval=t_eval_pos, method='RK45')

    # Tempo "negativo" (ex: 0 até -10)
    t_eval_neg = np.linspace(t_span[0], -t_span[1], num_steps)
    sol_neg = solve_ivp(rhs, (0, -t_span[1]), [u0, v0, z0],
                        t_eval=t_eval_neg, method='RK45')

    return sol_neg, sol_pos

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
    Gera N*N2 pontos na superfície de uma esfera de raio 'raio'
    em torno de (u_c, v_c, z_c).
    """
    thetas = np.linspace(0, 2*np.pi, N, endpoint=False)   # Azimute
    phis   = np.linspace(0, np.pi,   N2, endpoint=True)   # Colatitude

    pontos = []
    for theta in thetas:
        for phi in phis:
            u_i = u_c + raio * np.sin(phi) * np.cos(theta)
            v_i = v_c + raio * np.sin(phi) * np.sin(theta)
            z_i = z_c + raio * np.cos(phi)
            pontos.append([u_i, v_i, z_i])
    return pontos

if __name__ == "__main__":
    # Valores iniciais
    u_L, v_L, z_L = 0.1, 0.1, 0.1
    print(f"Valores iniciais: u_L = {u_L}, v_L = {v_L}, z_L = {z_L}")

    # Resolvendo trajetória "base" (L -> R) usando funcoes.resolver_trajetoria (já existente)
    sol_prim2, sol_prim = funcoes.resolver_trajetoria(u_L, v_L, z_L)
    print("Tamanho da lista de ul a ur", len(sol_prim.y[0]))

    # Pegamos o final da trajetória como ponto R
    u_R = sol_prim.y[0, 1000]
    v_R = sol_prim.y[1, 1000]
    z_R = sol_prim.y[2, 1000]
    distance = np.sqrt((u_R - u_L)**2 + (v_R - v_L)**2 + (z_R - z_L)**2)
    print(f"Ponto R = {u_R}, v_R = {v_R}, z_R = {z_R}")

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

    # -------------- Exemplo de gerar N pontos iniciais na esfera --------------
    N = 6
    N2 = 4
    raio = distance / 2
    print("Raio: ", raio)
    pontos_iniciais = gerar_pontos_iniciais(u_L, v_L, z_L, N, N2, raio)

    # ======= PLOTS 2D da trajetória base (já resolvida em funcoes.resolver_trajetoria) =======
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
    # Para cada ponto inicial da esfera, resolvemos as trajetórias
    # (0 -> +10) e (0 -> -10) e plotamos nos subplots 2D.
    # ================================================================
    for idx, (u0, v0, z0) in enumerate(pontos_iniciais):

        sol_neg, sol_pos = resolver_trajetoria(
            u0, v0, z0, u_R, v_R, z_R, f_R, g_R, sigmaLR,
            t_span=(0, 10), num_steps=600
        )

        # -- Extração das trajetórias (negativa e positiva)
        u_neg = sol_neg.y[0]
        v_neg = sol_neg.y[1]
        z_neg = sol_neg.y[2]

        u_pos = sol_pos.y[0]
        v_pos = sol_pos.y[1]
        z_pos = sol_pos.y[2]

        # Plotamos também a bolinha vermelha no ponto inicial
        ax_uz.scatter(u0, z0, color='red', s=10)
        ax_vz.scatter(v0, z0, color='red', s=10)
        ax_uv.scatter(u0, v0, color='red', s=10)

        # Filtro para plotar com menos setas (para não poluir)
        stepQ = 50

        # -- Quiver (u x z) para a parte NEGATIVA
        ax_uz.quiver(u_neg[:-1:stepQ], z_neg[:-1:stepQ],
                     (u_neg[1::stepQ] - u_neg[:-1:stepQ]),
                     (z_neg[1::stepQ] - z_neg[:-1:stepQ]),
                     angles='xy', scale_units='xy', scale=2, color='orange')

        # -- Quiver (u x z) para a parte POSITIVA
        ax_uz.quiver(u_pos[:-1:stepQ], z_pos[:-1:stepQ],
                     (u_pos[1::stepQ] - u_pos[:-1:stepQ]),
                     (z_pos[1::stepQ] - z_pos[:-1:stepQ]),
                     angles='xy', scale_units='xy', scale=2, color='green')

        # -- Quiver (v x z) NEGATIVO
        ax_vz.quiver(v_neg[:-1:stepQ], z_neg[:-1:stepQ],
                     (v_neg[1::stepQ] - v_neg[:-1:stepQ]),
                     (z_neg[1::stepQ] - z_neg[:-1:stepQ]),
                     angles='xy', scale_units='xy', scale=2, color='orange')

        # -- Quiver (v x z) POSITIVO
        ax_vz.quiver(v_pos[:-1:stepQ], z_pos[:-1:stepQ],
                     (v_pos[1::stepQ] - v_pos[:-1:stepQ]),
                     (z_pos[1::stepQ] - z_pos[:-1:stepQ]),
                     angles='xy', scale_units='xy', scale=2, color='green')

        # -- Quiver (u x v) NEGATIVO
        ax_uv.quiver(u_neg[:-1:stepQ], v_neg[:-1:stepQ],
                     (u_neg[1::stepQ] - u_neg[:-1:stepQ]),
                     (v_neg[1::stepQ] - v_neg[:-1:stepQ]),
                     angles='xy', scale_units='xy', scale=2, color='orange')

        # -- Quiver (u x v) POSITIVO
        ax_uv.quiver(u_pos[:-1:stepQ], v_pos[:-1:stepQ],
                     (u_pos[1::stepQ] - u_pos[:-1:stepQ]),
                     (v_pos[1::stepQ] - v_pos[:-1:stepQ]),
                     angles='xy', scale_units='xy', scale=2, color='green')

    plt.tight_layout()
    plt.show()

    # ====== PLOT 3D ======
    # Agora criamos um plot em 3D para visualizar a trajetória no espaço (u, v, z).
    fig_3d = plt.figure(figsize=(8, 6))
    ax_3d = fig_3d.add_subplot(111, projection='3d')

    # Trajetória "base" (U^L -> U^R) (já resolvida)
    ax_3d.plot(base_u, base_v, base_z, color='gray', label='Trajet. base (U^L)')

    # Marca pontos L e R
    ax_3d.scatter(u_L, v_L, z_L, color='red', s=50, label='Ponto L')
    ax_3d.scatter(u_R, v_R, z_R, color='blue', s=50, label='Ponto R')

    # Para cada ponto inicial (u0, v0, z0), plotamos as duas trajetórias (neg e pos).
    for idx, (u0, v0, z0) in enumerate(pontos_iniciais):
        sol_neg, sol_pos = resolver_trajetoria(
            u0, v0, z0, u_R, v_R, z_R, f_R, g_R, sigmaLR,
            t_span=(0, 10), num_steps=600
        )
        # Adiciona ponto inicial em vermelho
        ax_3d.scatter(u0, v0, z0, color='red', s=20)

        # Trajetória negativa
        ax_3d.plot(sol_neg.y[0], sol_neg.y[1], sol_neg.y[2], alpha=0.8, color='orange')
        # Trajetória positiva
        ax_3d.plot(sol_pos.y[0], sol_pos.y[1], sol_pos.y[2], alpha=0.8, color='green')

    ax_3d.set_xlabel('u')
    ax_3d.set_ylabel('v')
    ax_3d.set_zlabel('z')
    ax_3d.set_title('Trajetórias em 3D (u, v, z)')
    ax_3d.legend()
    plt.show()