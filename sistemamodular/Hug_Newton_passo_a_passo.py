#!/usr/bin/env python3
"""
    Hug_Newton_passo_a_passo_rescalado.py

Este script implementa o ramo da curva de Hugoniot do sistema perturbado pela adsorção
que aproxima a curva de contato via o método de Newton para três equações, conforme descrito
na Seção 4.2.8 do relatório.

Alterações:
  - Escalonamento da equação de arco-comprimento (H) dividindo por h.
  - Chute inicial (U_guess) sempre a uma distância ~h de Ua usando direção normalizada.

Autor: [Seu Nome]
Data: [Data]
"""

from matplotlib.pylab import LinAlgError
import numpy as np
import matplotlib.pyplot as plt

# Importa as funções definidas no módulo system:
# As funções F, G e seus derivados, além de f, g, a, alpha e system.
from system import F, G, dF_du, dF_dv, dF_dz, dG_du, dG_dv, dG_dz, f, g, a, alpha, system

# =============================================================================
# Equação de parametrização H escalonada e suas derivadas (dividida por h)
# =============================================================================
def H_scaled(U, Ua, h):
    # H/h = ((U - Ua)^2 - h^2)/h
    return ((U[0] - Ua[0])**2 + (U[1] - Ua[1])**2 + (U[2] - Ua[2])**2 - h**2) / h

def dH_scaled_du(U, Ua, h):
    return 2 * (U[0] - Ua[0]) / h

def dH_scaled_dv(U, Ua, h):
    return 2 * (U[1] - Ua[1]) / h

def dH_scaled_dz(U, Ua, h):
    return 2 * (U[2] - Ua[2]) / h

# =============================================================================
# Função que verifica se um ponto está no prisma definido:
# 0 <= u <= 1, 0 <= v <= 1, u+v <= 1, 0 <= z <= 1.
# =============================================================================
def in_prisma(U):
    u, v, z = U
    return (0 <= u <= 1) and (0 <= v <= 1) and ((u + v) <= 1) and (0 <= z <= 1)

# =============================================================================
# Monta o vetor de equações para o método de Newton: [F, G, H_scaled] = 0.
# =============================================================================
def FGH(U, Ua, h, U0):
    u, v, z = U
    u0, v0, z0 = U0
    # Valores de referência no ponto base:
    f0 = f(u0, v0, z0)
    g0 = g(u0, v0, z0)
    eqF = F(u, v, z, f0, v0, g0, u0)
    eqG = G(u, v, z, f0, z0, u0, a, alpha)
    eqH = H_scaled(U, Ua, h)
    return np.array([eqF, eqG, eqH])

# =============================================================================
# Monta a matriz Jacobiana do sistema [F, G, H_scaled] em U.
# =============================================================================
def jacobian_FGH(U, Ua, h, U0):
    u, v, z = U
    u0, v0, z0 = U0
    f0 = f(u0, v0, z0)
    g0 = g(u0, v0, z0)
    
    J = np.zeros((3, 3))
    # Derivadas parciais de F:
    J[0, 0] = dF_du(u, v, z, f0, v0, g0, u0)
    J[0, 1] = dF_dv(u, v, z, f0, v0, g0, u0)
    J[0, 2] = dF_dz(u, v, z, f0, v0, g0, u0)
    # Derivadas parciais de G:
    J[1, 0] = dG_du(u, v, z, f0, z0, u0, alpha)
    J[1, 1] = dG_dv(u, v, z, f0, z0, u0, alpha)
    J[1, 2] = dG_dz(u, v, z, f0, z0, u0, alpha)
    # Derivadas parciais de H_scaled:
    J[2, 0] = dH_scaled_du(U, Ua, h)
    J[2, 1] = dH_scaled_dv(U, Ua, h)
    J[2, 2] = dH_scaled_dz(U, Ua, h)
    return J

# =============================================================================
# Iteração de Newton para corrigir o palpite inicial
# =============================================================================
def newton_iteration(U_guess, Ua, h, U0, tol=1e-6, max_iter=20):
    U = U_guess.copy()
    for i in range(max_iter):
        F_val = FGH(U, Ua, h, U0)
        norm_F_val = np.linalg.norm(F_val)
        J = jacobian_FGH(U, Ua, h, U0)
        detJ  = np.linalg.det(J)
        condJ = np.linalg.cond(J)
        print(f"  it={i:2d} | ||F||={norm_F_val:.3e} | det(J)={detJ:.3e} | cond(J)={condJ:.3e}")
        if norm_F_val < tol:
            print(f"  → convergiu em {i} iterações\n")
            return U, True, i+1
        try:
            delta = np.linalg.solve(J, -F_val)
        except LinAlgError:
            print(f"  → Jacobiana singular (cond={condJ:.3e}). Abortando.\n")
            return U, False, i+1
        U += delta
    print(f"  → não convergiu após {max_iter} iterações (último cond={condJ:.3e})\n")
    return U, False, max_iter

# =============================================================================
# Função principal que implementa o algoritmo (Seção 4.2.8)
# =============================================================================
def main():
    U0 = np.array([0.42463291, 0.52495654, 0.50199074])
    h_step = 0.01                   # Passo de comprimento de arco desejado
    tol_newton = 0.001              # Tolerância para Newton
    max_iter_newton = 5             # Máx iterações de Newton
    num_steps = int(1.0/h_step)
    Upos = [U0.copy()]
    print('main: U0 =', Upos[0])

    # Primeiro palpite por Euler
    Y = system(h_step, U0)
    print('main: sistema 19 (P, Q, R): Y =', Y)
    U_1 = U0 + h_step * Y
    print('main: U_1 =', U_1)
    Ua = U0.copy()

    # Correção de U_1
    U_corr, converged, iters = newton_iteration(U_1, Ua, h_step, U0,
                                                tol=tol_newton, max_iter=max_iter_newton)
    if not converged or not in_prisma(U_corr):
        print("Falha na correção inicial.")
        return
    Upos.append(U_corr.copy())

    U_prev_prev = U0.copy()
    Ua = U_corr.copy()

    for step in range(1, num_steps+1):
        print(f"\n main: passo {step}")
        # Extrapolação normalizada:
        direction = Ua - U_prev_prev
        norm_dir = np.linalg.norm(direction)
        if norm_dir == 0:
            break
        direction /= norm_dir
        U_guess = Ua + h_step * direction

        U_new, converged, iters = newton_iteration(U_guess, Ua, h_step, U0,
                                                   tol=tol_newton, max_iter=max_iter_newton)
        if not converged or not in_prisma(U_new):
            print(f"Falha ou fora do prisma no passo {step}.")
            break
        Upos.append(U_new.copy())
        U_prev_prev, Ua = Ua, U_new.copy()

    # Plot 3D
    Upos = np.array(Upos)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(Upos[0,0], Upos[0,1], Upos[0,2], 'ro', label='Inicial')
    ax.plot(Upos[-1,0], Upos[-1,1], Upos[-1,2], 'bo', label='Final')
    ax.plot(Upos[:,0], Upos[:,1], Upos[:,2], 'm-', label='Hugoniot')
    ax.set_xlabel('u'); ax.set_ylabel('v'); ax.set_zlabel('z')
    ax.set_title('Hugoniot via Newton (H escalonado)')
    ax.legend()
    plt.show()

if __name__ == "__main__":
    main()
