import numpy as np
from scipy.integrate import solve_ivp
from system import system_with_determinants

def resolver_trajetoria(u0, v0, c0, s_span=(0, 10), s_span2=(0, -10)):
    # Get initial values for positive and negative spans
    initial_pos = solve_ivp(system_with_determinants, s_span, [u0, v0, c0], method='LSODA', t_eval=np.linspace(s_span[0], s_span[1], 20000))
    initial_neg = solve_ivp(system_with_determinants, s_span2, [u0, v0, c0], method='LSODA', t_eval=np.linspace(s_span2[0], s_span2[1], 20000))
    return initial_pos, initial_neg

def dentro_do_triangulo(u, v, c):
    return u >= 0 and v >= 0 and u + v <= 1 and 0 <= c <= 1

def dividir_trajetorias(sol):
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
