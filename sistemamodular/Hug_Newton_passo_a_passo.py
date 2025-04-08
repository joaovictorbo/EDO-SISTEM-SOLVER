#!/usr/bin/env python3
"""
Hug_Newton_passo_a_passo.py

Este script implementa o ramo da curva de Hugoniot do sistema perturbado pela adsorção
que aproxima a curva de contato via o método de Newton para três equações, conforme descrito
na Seção 4.2.8 do relatório.

O algoritmo adotado é:
  (1) Obter um palpite inicial para U1 integrando o sistema (19) a partir de U0 (função system).
  (2) Corrigir o palpite usando iterações de Newton (função newton_iteration) para resolver:
          F(u,v,z)=0, G(u,v,z)=0, H(u,v,z; Ua, h)=0.
  (3) Para o próximo passo, usar a extrapolação U_guess = U_prev + h*(U_prev - U_prev_prev),
      com U_prev sendo o ponto corrigido do passo anterior, e corrigir com Newton.
  (4) Repetir até que o ramo esteja completo (ou o ponto saia do domínio).

Autor: [Seu Nome]
Data: [Data]
"""

import numpy as np
import matplotlib.pyplot as plt

# Importa as funções definidas no módulo system:
# As funções F, G e seus derivados, além de f, g, a, alpha e system.
from system import F, G, dF_du, dF_dv, dF_dz, dG_du, dG_dv, dG_dz, f, g, a, alpha, system

# =============================================================================
# Equação de parametrização H e suas derivadas (conforme (45) e (46)-(48))
# =============================================================================
def H_equac(U, Ua, h):
    """
    Calcula a equação de parametrização pelo comprimento de arco:
      H(u,v,z; ua,va,za, h) = (u - ua)**2 + (v - va)**2 + (z - za)**2 - h**2
    """
    return (U[0] - Ua[0])**2 + (U[1] - Ua[1])**2 + (U[2] - Ua[2])**2 - h**2

def dH_du(U, Ua):
    return 2 * (U[0] - Ua[0])

def dH_dv(U, Ua):
    return 2 * (U[1] - Ua[1])

def dH_dz(U, Ua):
    return 2 * (U[2] - Ua[2])

# =============================================================================
# Função que verifica se um ponto está no prisma definido:
# 0 <= u <= 1, 0 <= v <= 1, u+v <= 1, 0 <= z <= 1.
# =============================================================================
def in_prisma(U):
    u, v, z = U
    return (0 <= u <= 1) and (0 <= v <= 1) and ((u + v) <= 1) and (0 <= z <= 1)

# =============================================================================
# Monta o vetor de equações para o método de Newton: [F, G, H] = 0.
# Aqui, os parâmetros de referência (f0, g0, u0, v0, z0) são extraídos do ponto base U0.
# =============================================================================
def FGH(U, Ua, h, U0):
    u, v, z = U
    u0, v0, z0 = U0
    # Valores de referência no ponto base:
    f0 = f(u0, v0, z0)
    g0 = g(u0, v0, z0)
    eqF = F(u, v, z, f0, v0, g0, u0)
    eqG = G(u, v, z, f0, z0, u0, a, alpha)
    eqH = H_equac(U, Ua, h)
    return np.array([eqF, eqG, eqH])

# =============================================================================
# Monta a matriz Jacobiana do sistema [F, G, H] em U.
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
    # Derivadas parciais de H:
    J[2, 0] = dH_du(U, Ua)
    J[2, 1] = dH_dv(U, Ua)
    J[2, 2] = dH_dz(U, Ua)
    return J

# =============================================================================
# Iteração de Newton para corrigir o palpite inicial
# =============================================================================
def newton_iteration(U_guess, Ua, h, U0, tol=1e-6, max_iter=20):
    """
    Aplica o método de Newton para resolver o sistema:
         F(U, U0) = 0,  G(U, U0) = 0,  H(U, Ua, h) = 0.
    
    Parâmetros:
      U_guess  : palpite inicial para U (vetor [u, v, z])
      Ua       : ponto base para a equação de comprimento de arco (usualmente o último ponto corrigido)
      h        : passo (comprimento de arco desejado)
      U0       : ponto inicial fixo do ramo (para referência em F e G)
      tol      : tolerância para convergência
      max_iter : número máximo de iterações
      
    Retorna:
      U_corr, converged, iter_used
    """
    U = U_guess.copy()
    for i in range(max_iter):
        F_val = FGH(U, Ua, h, U0)
        norm_F = np.linalg.norm(F_val)
        if norm_F < tol:
            return U, True, i+1
        J = jacobian_FGH(U, Ua, h, U0)
        try:
            delta = np.linalg.solve(J, -F_val)
        except np.linalg.LinAlgError:
            print("Matriz singular na iteração de Newton.")
            return U, False, i+1
        U += delta
    return U, False, max_iter

# =============================================================================
# Função principal que implementa o algoritmo (seguindo a seção 4.2.8)
# =============================================================================
def main():
    # Ponto inicial no prisma
    U0 = np.array([0.1, 0.6, 0.2])  # U0 = (u0, v0, z0)
    h_step = 0.01                   # Passo de comprimento de arco desejado
    tol_newton = 1e-6               # Tolerância para o método de Newton
    max_iter_newton = 10        # Número máximo de iterações de Newton por passo
    num_steps = 100                 # Número máximo de pontos do ramo
    
    # Lista para armazenar os pontos do ramo (para h > 0)
    Upos = [U0.copy()]
    
    # -------------------------------------------------------------------------
    # Primeiro passo:
    # Obtenha o palpite inicial U_corr para U1 integrando o sistema (19) a partir de U0.
    # -------------------------------------------------------------------------
    U_integrado = system(h_step, U0)  # system: gera um palpite inicial aproximado
    # Use U0 como base para a correção do primeiro ponto
    U_corr = U0 + h_step * U_integrado # Um passo de integracao por Euler para o primeiro ponto da curva
                                       # a ser corrigido por Newton
    Ua = U0
    U_corr, converged, iters = newton_iteration(U_corr, Ua, h_step, U0,
                                                tol=tol_newton, max_iter=max_iter_newton)
    if not converged:
        print("Newton não convergiu no primeiro passo.")
        return
    if not in_prisma(U_corr):
        print("O ponto corrigido não está no prisma.")
        return
    Upos.append(U_corr.copy())
    
    # Para os próximos passos:
    # U_prev_prev guarda o penúltimo ponto corrigido e U_prev o último.
    U_prev_prev = U0.copy() # Passo (5) do documento 
    U_a = U_corr.copy() # Novo Ua do passo (5)
    
    for step in range(1, num_steps):
        # Gera o palpite para o próximo ponto usando extrapolação linear:
        U_guess = U_a+ h_step * (U_a - U_prev_prev)
        # Use o último ponto corrigido (U_prev) como base na correção (congela o valor de z, por exemplo)
        U_new, converged, iters = newton_iteration(U_guess, Ua, h_step, U0,
                                                   tol=tol_newton, max_iter=max_iter_newton)
        if not converged:
            print(f"Newton não convergiu no passo {step}. Encerrando iteração.")
            break
        if not in_prisma(U_new):
            print(f"Ponto fora do prisma no passo {step}. Encerrando iteração.")
            break
        Upos.append(U_new.copy())
        # Atualiza os pontos para a próxima extrapolação
        U_prev_prev = U_a
        U_a = U_new.copy()
    
    # Converte a lista de pontos para um array para a plotagem
    Upos = np.array(Upos)
    
    # Plotagem do ramo da curva de Hugoniot obtido (h > 0)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(Upos[:, 0], Upos[:, 1], Upos[:, 2], 'o-', label='Ramo Hugoniot (Newton 3 eq)')
    ax.set_xlabel('u')
    ax.set_ylabel('v')
    ax.set_zlabel('z')
    ax.set_title('Curva de Hugoniot via Newton para 3 equações')
    ax.legend()
    plt.show()

if __name__ == "__main__":
    main()