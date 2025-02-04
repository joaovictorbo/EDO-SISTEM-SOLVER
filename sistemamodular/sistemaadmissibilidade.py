import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Import necessário para plot 3D
import newton_funcoes as funcoes 
from jacobiana import calcular_jacobiana

# ------------------- PARÂMETROS GLOBAIS -------------------
alpha = 0.001
epsilon_1 = 0.1
epsilon_2 = 0.01
muw0 = 1.0  # Viscosidade inicial sem polímero

# ------------------- FUNÇÕES DOS MODELOS -------------------
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

# ------------------- SISTEMA E RESOLUÇÃO -------------------
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
    def rhs(s, y):
        return system(s, y, u_R, v_R, z_R, f_R, g_R, sigma)

    t_eval_pos = np.linspace(t_span[0], t_span[1], num_steps)
    sol_pos = solve_ivp(rhs, t_span, [u0, v0, z0],
                        t_eval=t_eval_pos, method='RK45')

    t_eval_neg = np.linspace(t_span[0], -t_span[1], num_steps)
    sol_neg = solve_ivp(rhs, (0, -t_span[1]), [u0, v0, z0],
                        t_eval=t_eval_neg, method='RK45')
    return sol_neg, sol_pos

def calcular_autovalores_e_autovetores(jacobiana):
    autovalores, autovetores = np.linalg.eig(jacobiana)
    indices_positivos = np.where(autovalores > 0)[0]
    indices_negativos = np.where(autovalores < 0)[0]
    return ((autovalores[indices_positivos], autovetores[:, indices_positivos]),
            (autovalores[indices_negativos], autovetores[:, indices_negativos]))

def gerar_pontos_iniciais(u_c, v_c, z_c, N, N2, raio):
    thetas = np.linspace(0, 2*np.pi, N, endpoint=False)
    phis   = np.linspace(0, np.pi,   N2, endpoint=True)
    pontos = []
    for theta in thetas:
        for phi in phis:
            u_i = u_c + raio * np.sin(phi) * np.cos(theta)
            v_i = v_c + raio * np.sin(phi) * np.sin(theta)
            z_i = z_c + raio * np.cos(phi)
            pontos.append([u_i, v_i, z_i])
    return pontos

# ------------------- FUNÇÃO PARA PLOTAR A "SETA" USANDO AS FÓRMULAS DO TEXTO -------------------
def plot_arrow_structure_text(P, Q, d, ax, r=None, color='purple', lw=0.5):
    """
    Dadas as coordenadas de P = (uP, vP, zP) e Q = (uQ, vQ, zQ), calcula o ponto
    M no segmento PQ, distante d do extremo Q (0 < d ≤ 1/3):
         M = P + (1-d)*(Q-P) = Q - d*(Q-P)
    e, usando um raio r (se não informado, r = d/3), constrói a circunferência
    no plano Π ortogonal a PQ (definido por:
         (uQ-uP)(u-uM) + (vQ-vP)(v-vM) + (zQ-zP)(z-zM) = 0)
    com as seguintes equações paramétricas:
      • Caso 1: se |zQ-zP| > tol:
             u = uM + r*cos(θ)
             v = vM + r*sin(θ)
             z = zM - (r/(zQ-zP))*[(uQ-uP)*cos(θ) + (vQ-vP)*sin(θ)]
      • Caso 2: se |zQ-zP| < tol e |vQ-vP| > tol:
             u = uM + r*cos(θ)
             v = vM - ((uQ-uP)/(vQ-vP))*r*cos(θ)
             z = zM + r*sin(θ)
      • Caso 3: se |zQ-zP| < tol e |vQ-vP| < tol:
             u = uM
             v = vM + r*cos(θ)
             z = zM + r*sin(θ)
    
    Em seguida, traça os segmentos entre os quatro pontos correspondentes a 
    θ = 0, π/2, π e 3π/2, formando o quadrilátero, e os segmentos ligando cada
    um desses pontos a Q.
    """
    tol = 1e-6
    P = np.array(P, dtype=float)
    Q = np.array(Q, dtype=float)
    D = Q - P
    M = Q - d * D
    if r is None:
        # Reduzindo o tamanho da seta: em vez de d/3, usamos d/6
        r = d / 6.0

    uP, vP, zP = P
    uQ, vQ, zQ = Q
    uM, vM, zM = M

    circle_pts = []
    angles = [0, np.pi/2, np.pi, 3*np.pi/2]

    if abs(zQ - zP) > tol:  # Caso 1
        for theta in angles:
            u_val = uM + r * np.cos(theta)
            v_val = vM + r * np.sin(theta)
            z_val = zM - (r/(zQ - zP)) * ((uQ - uP)*np.cos(theta) + (vQ - vP)*np.sin(theta))
            circle_pts.append([u_val, v_val, z_val])
    elif abs(zQ - zP) <= tol and abs(vQ - vP) > tol:  # Caso 2
        for theta in angles:
            u_val = uM + r * np.cos(theta)
            v_val = vM - ((uQ - uP)/(vQ - vP)) * r * np.cos(theta)
            z_val = zM + r * np.sin(theta)
            circle_pts.append([u_val, v_val, z_val])
    else:  # Caso 3
        for theta in angles:
            u_val = uM
            v_val = vM + r * np.cos(theta)
            z_val = zM + r * np.sin(theta)
            circle_pts.append([u_val, v_val, z_val])
    circle_pts = np.array(circle_pts)

    for i in range(len(circle_pts)):
        j = (i + 1) % len(circle_pts)
        ax.plot([circle_pts[i, 0], circle_pts[j, 0]],
                [circle_pts[i, 1], circle_pts[j, 1]],
                [circle_pts[i, 2], circle_pts[j, 2]],
                color=color, lw=lw)
    for pt in circle_pts:
        ax.plot([pt[0], Q[0]],
                [pt[1], Q[1]],
                [pt[2], Q[2]],
                color=color, lw=lw)
    # Removemos a marcação dos pontos M e Q (não se quer os pontos iniciais vermelhos)
    # ax.scatter(M[0], M[1], M[2], color=color, s=50)
    # ax.scatter(Q[0], Q[1], Q[2], color='red', s=50)

# ------------------- BLOCO PRINCIPAL -------------------
if __name__ == "__main__":
    # Valores iniciais
    u_L, v_L, z_L = 0.1, 0.1, 0.1
    print(f"Valores iniciais: u_L = {u_L}, v_L = {v_L}, z_L = {z_L}")

    sol_prim2, sol_prim = funcoes.resolver_trajetoria(u_L, v_L, z_L)
    print("Tamanho da lista de ul a ur", len(sol_prim.y[0]))

    u_R = sol_prim.y[0, 1000]
    v_R = sol_prim.y[1, 1000]
    z_R = sol_prim.y[2, 1000]
    distance = np.sqrt((u_R - u_L)**2 + (v_R - v_L)**2 + (z_R - z_L)**2)
    print(f"Ponto R = {u_R}, v_R = {v_R}, z_R = {z_R}")

    f_L = f(u_L, v_L, z_L)
    g_L = g(u_L, v_L, z_L)
    f_R = f(u_R, v_R, z_R)
    g_R = g(u_R, v_R, z_R)

    sigmaLR = sigma_alpha(u_L, f_R, z_R, z_L)
    print(f"sigmaLR = {sigmaLR}")
    if u_R != 0:
        print(f"F/u: {f_R/u_R}")
    else:
        print("Divisão por zero em u_R para cálculo de F/u.")

    jacoL = calcular_jacobiana(u_L, v_L, z_L, f_L, g_L, sigmaLR, epsilon_1, epsilon_2, alpha)
    jacoR = calcular_jacobiana(u_R, v_R, z_R, f_R, g_R, sigmaLR, epsilon_1, epsilon_2, alpha)

    np.set_printoptions(precision=3, suppress=True)
    print("Jacobiana em U^L:")
    print(jacoL)
    print("\nJacobiana em U^R:")
    print(jacoR)

    (autovalores_pos_L, autovetores_pos_L), (autovalores_neg_L, autovetores_neg_L) = calcular_autovalores_e_autovetores(jacoL)
    (autovalores_pos_R, autovetores_pos_R), (autovalores_neg_R, autovetores_neg_R) = calcular_autovalores_e_autovetores(jacoR)

    print("Autovalores positivos de J(U^L):", autovalores_pos_L)
    print("Autovetores positivos de J(U^L):\n", autovetores_pos_L)
    print("\nAutovalores negativos de J(U^R):", autovalores_neg_R)
    print("Autovetores negativos de J(U^R):\n", autovetores_neg_R)

    N = 6
    N2 = 4
    raio = distance / 2
    print("Raio: ", raio)
    pontos_iniciais = gerar_pontos_iniciais(u_L, v_L, z_L, N, N2, raio)

    # ===== PLOTS 2D (trajetória base) =====
    fig, axes = plt.subplots(1, 3, figsize=(12, 5))
    step = 5  
    base_u = sol_prim.y[0, ::step]
    base_v = sol_prim.y[1, ::step]
    base_z = sol_prim.y[2, ::step]

    ax_uz = axes[0]
    ax_uz.quiver(base_u[:-1], base_z[:-1],
                 base_u[1:] - base_u[:-1],
                 base_z[1:] - base_z[:-1],
                 angles='xy', scale_units='xy', scale=1,
                 color='gray', label='Trajet. base (U^L)')
    # Removidos os pontos vermelhos iniciais:
    # ax_uz.scatter([u_L], [z_L], color='red', label='Ponto L')
    ax_uz.scatter([u_R], [z_R], color='blue', label='Ponto R')
    ax_uz.set_xlabel('u')
    ax_uz.set_ylabel('z')
    ax_uz.set_title('u x z')
    ax_uz.legend()

    ax_vz = axes[1]
    ax_vz.quiver(base_v[:-1], base_z[:-1],
                 base_v[1:] - base_v[:-1],
                 base_z[1:] - base_z[:-1],
                 angles='xy', scale_units='xy', scale=1,
                 color='gray', label='Trajet. base (U^L)')
    # ax_vz.scatter([v_L], [z_L], color='red', label='Ponto L')
    ax_vz.scatter([v_R], [z_R], color='blue', label='Ponto R')
    ax_vz.set_xlabel('v')
    ax_vz.set_ylabel('z')
    ax_vz.set_title('v x z')
    ax_vz.legend()

    ax_uv = axes[2]
    ax_uv.quiver(base_u[:-1], base_v[:-1],
                 base_u[1:] - base_u[:-1],
                 base_v[1:] - base_v[:-1],
                 angles='xy', scale_units='xy', scale=1,
                 color='gray', label='Trajet. base (U^L)')
    # ax_uv.scatter([u_L], [v_L], color='red', label='Ponto L')
    ax_uv.scatter([u_R], [v_R], color='blue', label='Ponto R')
    ax_uv.set_xlabel('u')
    ax_uv.set_ylabel('v')
    ax_uv.set_title('u x v')
    ax_uv.legend()

    plt.tight_layout()
    plt.show()

    # ===== PLOT 3D com setas (a cada 10 pontos) =====
    fig_3d = plt.figure(figsize=(8, 6))
    ax_3d = fig_3d.add_subplot(111, projection='3d')

    # Trajetória base (U^L -> U^R)
    ax_3d.plot(base_u, base_v, base_z, color='gray', label='Trajet. base (U^L)')
    d_param = 0.2   # 0 < d ≤ 1/3
    arrow_step = 10
    # Desenha a seta na trajetória base a cada 10 pontos
    for i in range(0, len(base_u) - arrow_step, arrow_step):
        P_arrow = (base_u[i], base_v[i], base_z[i])
        Q_arrow = (base_u[i + arrow_step], base_v[i + arrow_step], base_z[i + arrow_step])
        # Passa r = d/6 e lw menor para setas menores
        plot_arrow_structure_text(P_arrow, Q_arrow, d_param, ax_3d, r=d_param/6.0, color='gray', lw=0.5)

    # Removidos os pontos vermelhos iniciais; mantemos apenas R (ponto azul)
    # ax_3d.scatter(u_L, v_L, z_L, color='red', s=50, label='Ponto L')
    ax_3d.scatter(u_R, v_R, z_R, color='blue', s=50, label='Ponto R')

    for idx, (u0, v0, z0) in enumerate(pontos_iniciais):
        sol_neg, sol_pos = resolver_trajetoria(
            u0, v0, z0, u_R, v_R, z_R, f_R, g_R, sigmaLR,
            t_span=(0, 10), num_steps=600
        )
        # Removida a marcação do ponto inicial vermelho:
        # ax_3d.scatter(u0, v0, z0, color='red', s=20)
        ax_3d.plot(sol_neg.y[0], sol_neg.y[1], sol_neg.y[2], alpha=0.8, color='orange')
        for i in range(0, sol_neg.y.shape[1]-arrow_step, arrow_step):
            P_arrow = (sol_neg.y[0, i], sol_neg.y[1, i], sol_neg.y[2, i])
            Q_arrow = (sol_neg.y[0, i+arrow_step], sol_neg.y[1, i+arrow_step], sol_neg.y[2, i+arrow_step])
            plot_arrow_structure_text(P_arrow, Q_arrow, d_param, ax_3d, r=d_param/6.0, color='orange', lw=0.5)
        ax_3d.plot(sol_pos.y[0], sol_pos.y[1], sol_pos.y[2], alpha=0.8, color='green')
        for i in range(0, sol_pos.y.shape[1]-arrow_step, arrow_step):
            P_arrow = (sol_pos.y[0, i], sol_pos.y[1, i], sol_pos.y[2, i])
            Q_arrow = (sol_pos.y[0, i+arrow_step], sol_pos.y[1, i+arrow_step], sol_pos.y[2, i+arrow_step])
            plot_arrow_structure_text(P_arrow, Q_arrow, d_param, ax_3d, r=d_param/6.0, color='green', lw=0.5)

    ax_3d.set_xlabel('u')
    ax_3d.set_ylabel('v')
    ax_3d.set_zlabel('z')
    ax_3d.set_title('Trajetórias em 3D (u, v, z) com setas a cada 10 pontos')
    ax_3d.legend()
    plt.show()