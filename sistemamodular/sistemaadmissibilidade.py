import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Import necessário para plot 3D
import newton_funcoes as funcoes 
from jacobiana import calcular_jacobiana

def dentro_do_triangulo(u, v, c):
    return u >= 0 and v >= 0 and u + v <= 1 and 0 <= c <= 1

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
    return u**2 * (-muwz(z) / muw(z)**2 * D(u, v, z) - D_dz(u, v, z) / muw(z)) / (D(u, v, z)**2)

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

def system(s, y, u_R, v_R, z_R, f_R, g_R, sigma):
    u, v, z = y
    du_ds = f(u, v, z) - f_R - sigma * (u - u_R)
    dv_ds = g(u, v, z) - g_R - sigma * (v - v_R)
    dz_ds = (epsilon_1 / epsilon_2) * (f_R * (z - z_R) - sigma * (u_R * (z - z_R) - alpha * (a(z) - a(z_R))))
    return [du_ds, dv_ds, dz_ds]

def rk4_step(s, y, dt, u_R, v_R, z_R, f_R, g_R, sigma):
    k1 = np.array(system(s, y, u_R, v_R, z_R, f_R, g_R, sigma))
    k2 = np.array(system(s + dt/2, y + dt*k1/2, u_R, v_R, z_R, f_R, g_R, sigma))
    k3 = np.array(system(s + dt/2, y + dt*k2/2, u_R, v_R, z_R, f_R, g_R, sigma))
    k4 = np.array(system(s + dt, y + dt*k3, u_R, v_R, z_R, f_R, g_R, sigma))
    return y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

def resolver_trajetoria(u0, v0, z0, u_R, v_R, z_R, f_R, g_R, sigma,
                        t_span=(0,10), num_steps=10000):
    """
    Resolve o IVP de forma iterativa usando RK4.
    Para cada iteração do laço (num_steps), são realizados 2 passos de integração.
    Em cada passo há uma verificação (via dentro_do_triangulo) para garantir que o ponto permanece
    na região definida (u>=0, v>=0, u+v<=1 e 0<=z<=1).
    
    Parâmetros:
      - u0, v0, z0: condição inicial.
      - u_R, v_R, z_R, f_R, g_R, sigma: parâmetros usados na função do sistema.
      - t_span: intervalo de integração (usado para definir o tamanho do passo).
      - num_steps: número de iterações do laço (cada iteração realiza 2 passos, totalizando 2*num_steps passos).
    
    Retorna:
      - sol_neg: array dos pontos integrados para s negativo.
      - sol_pos: array dos pontos integrados para s positivo.
    """
    dt_pos = (t_span[1] - t_span[0])/(2*num_steps)
    dt_neg = - (t_span[1] - t_span[0])/(2*num_steps)

    def integrate_direction(dt, n_steps):
        s = 0.0
        y = np.array([u0, v0, z0])
        sol = [y.copy()]
        for i in range(n_steps):
            for j in range(2):  # Dois passos por iteração
                y = rk4_step(s, y, dt, u_R, v_R, z_R, f_R, g_R, sigma)
                s += dt
                # Verifica se o novo ponto permanece dentro do triângulo
                if not dentro_do_triangulo(y[0], y[1], y[2]):
                    return np.array(sol)
                sol.append(y.copy())
        return np.array(sol)

    sol_pos = integrate_direction(dt_pos, num_steps)
    sol_neg = integrate_direction(dt_neg, num_steps)
    return sol_neg, sol_pos

def calcular_autovalores_e_autovetores(jacobiana):
    autovalores, autovetores = np.linalg.eig(jacobiana)
    indices_positivos = np.where(autovalores > 0)[0]
    indices_negativos = np.where(autovalores < 0)[0]
    return ((autovalores[indices_positivos], autovetores[:, indices_positivos]),
            (autovalores[indices_negativos], autovetores[:, indices_negativos]))

def gerar_pontos_iniciais(u_c, v_c, z_c, N, N2, raio):
    thetas = np.linspace(0, 2*np.pi, N, endpoint=False)   # Azimute
    phis   = np.linspace(0, np.pi,   N2, endpoint=True)    # Colatitude

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

    # Resolvendo trajetória "base" (L -> R) usando a função já existente em funcoes
    sol_prim2, sol_prim = funcoes.resolver_trajetoria(u_L, v_L, z_L)
    print("Tamanho da lista de u de L a R:", len(sol_prim.y[0]))

    # Pegamos o final da trajetória como ponto R
    u_R = sol_prim.y[0, 1000]
    v_R = sol_prim.y[1, 1000]
    z_R = sol_prim.y[2, 1000]
    distance = np.sqrt((u_R - u_L)**2 + (v_R - v_L)**2 + (z_R - z_L)**2)
    print(f"Ponto R = {u_R}, {v_R}, {z_R}")

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

    # Gerar pontos iniciais na esfera
    N = 6
    N2 = 4
    raio = distance / 2
    print("Raio:", raio)
    pontos_iniciais = gerar_pontos_iniciais(u_L, v_L, z_L, N, N2, raio)

    # ======= PLOTS 2D da trajetória base (já resolvida) =======
    fig, axes = plt.subplots(1, 3, figsize=(12, 5))
    step = 5  
    base_u = sol_prim.y[0, ::step]
    base_v = sol_prim.y[1, ::step]
    base_z = sol_prim.y[2, ::step]

    # Subplot (u x z)
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

    # Subplot (v x z)
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

    # Subplot (u x v)
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

    # Para cada ponto inicial da esfera, resolvemos as trajetórias (negativa e positiva)
    for idx, (u0, v0, z0) in enumerate(pontos_iniciais):
        sol_neg, sol_pos = resolver_trajetoria(
            u0, v0, z0, u_R, v_R, z_R, f_R, g_R, sigmaLR,
            t_span=(0, 10), num_steps=20000
        )
        u_neg = sol_neg[:, 0]
        v_neg = sol_neg[:, 1]
        z_neg = sol_neg[:, 2]

        u_pos = sol_pos[:, 0]
        v_pos = sol_pos[:, 1]
        z_pos = sol_pos[:, 2]

        ax_uz.scatter(u0, z0, color='red', s=10)
        ax_vz.scatter(v0, z0, color='red', s=10)
        ax_uv.scatter(u0, v0, color='red', s=10)

        stepQ = 50

        ax_uz.quiver(u_neg[:-1:stepQ], z_neg[:-1:stepQ],
                     (u_neg[1::stepQ] - u_neg[:-1:stepQ]),
                     (z_neg[1::stepQ] - z_neg[:-1:stepQ]),
                     angles='xy', scale_units='xy', scale=2, color='orange')
        ax_uz.quiver(u_pos[:-1:stepQ], z_pos[:-1:stepQ],
                     (u_pos[1::stepQ] - u_pos[:-1:stepQ]),
                     (z_pos[1::stepQ] - z_pos[:-1:stepQ]),
                     angles='xy', scale_units='xy', scale=2, color='green')

        ax_vz.quiver(v_neg[:-1:stepQ], z_neg[:-1:stepQ],
                     (v_neg[1::stepQ] - v_neg[:-1:stepQ]),
                     (z_neg[1::stepQ] - z_neg[:-1:stepQ]),
                     angles='xy', scale_units='xy', scale=2, color='orange')
        ax_vz.quiver(v_pos[:-1:stepQ], z_pos[:-1:stepQ],
                     (v_pos[1::stepQ] - v_pos[:-1:stepQ]),
                     (z_pos[1::stepQ] - z_pos[:-1:stepQ]),
                     angles='xy', scale_units='xy', scale=2, color='green')

        ax_uv.quiver(u_neg[:-1:stepQ], v_neg[:-1:stepQ],
                     (u_neg[1::stepQ] - u_neg[:-1:stepQ]),
                     (v_neg[1::stepQ] - v_neg[:-1:stepQ]),
                     angles='xy', scale_units='xy', scale=2, color='orange')
        ax_uv.quiver(u_pos[:-1:stepQ], v_pos[:-1:stepQ],
                     (u_pos[1::stepQ] - u_pos[:-1:stepQ]),
                     (v_pos[1::stepQ] - v_pos[:-1:stepQ]),
                     angles='xy', scale_units='xy', scale=2, color='green')

    plt.tight_layout()
    plt.show()

    # ======= PLOT 3D =======
    fig_3d = plt.figure(figsize=(8, 6))
    ax_3d = fig_3d.add_subplot(111, projection='3d')
    ax_3d.plot(base_u, base_v, base_z, color='gray', label='Trajet. base (U^L)')
    ax_3d.scatter(u_L, v_L, z_L, color='red', s=50, label='Ponto L')
    ax_3d.scatter(u_R, v_R, z_R, color='blue', s=50, label='Ponto R')

    for idx, (u0, v0, z0) in enumerate(pontos_iniciais):
        sol_neg, sol_pos = resolver_trajetoria(
            u0, v0, z0, u_R, v_R, z_R, f_R, g_R, sigmaLR,
            t_span=(0, 100), num_steps=5000
        )
        ax_3d.scatter(u0, v0, z0, color='red', s=20)
        ax_3d.plot(sol_neg[:,0], sol_neg[:,1], sol_neg[:,2], alpha=0.8, color='orange')
        ax_3d.plot(sol_pos[:,0], sol_pos[:,1], sol_pos[:,2], alpha=0.8, color='green')

    ax_3d.set_xlabel('u')
    ax_3d.set_ylabel('v')
    ax_3d.set_zlabel('z')
    ax_3d.set_title('Trajetórias em 3D (u, v, z)')
    ax_3d.legend()
    plt.show()