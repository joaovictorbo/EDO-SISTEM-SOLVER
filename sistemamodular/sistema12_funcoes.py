# sistema12_funcoes.py
import numpy as np
from scipy.integrate import solve_ivp
from system import system

def resolver_trajetoria(u0, v0, c0, t_span=(0, 10), t_span2=(0, -10)):
    """
    Resolve a trajetória a partir das condições iniciais fornecidas.
    
    Parâmetros:
        u0, v0, c0: Condições iniciais para u, v, c.
        t_span, t_span2: Intervalos de tempo para integração positiva e negativa.
    
    Retorna:
        sol, sol2: Soluções para a trajetória positiva e negativa.
    """
    y0 = [u0, v0, c0]
    sol = solve_ivp(system, t_span, y0, method='LSODA', t_eval=np.linspace(t_span[0], t_span[1], 2000))
    sol2 = solve_ivp(system, t_span2, y0, method='LSODA', t_eval=np.linspace(t_span2[0], t_span2[1], 2000))
    return sol, sol2

def dentro_do_triangulo(u, v, c):
    return u >= 0 and v >= 0 and u + v <= 1 and 0 <= c <= 1

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
