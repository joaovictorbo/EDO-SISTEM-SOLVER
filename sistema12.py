import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from sistemamodular.system import system

u0, v0, c0 = 0.3,0.6, 0.4
y0 = [u0, v0, c0]
t_span = (0,10)
sol = solve_ivp(system, t_span, y0, method='LSODA', t_eval=np.linspace(0,10,2000))
t_span2 = (0,-10)
sol2 = solve_ivp(system, t_span2, y0, method='LSODA', t_eval=np.linspace(0,-10,2000))

# Função para verificar se um ponto está dentro do triângulo
def dentro_do_triangulo(u, v, c):
    return u >= 0 and v >= 0 and u + v <= 1 and 0 <= c <= 1

# Função para dividir as trajetórias ao entrar e sair do triângulo
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

# Dividir as trajetórias
trajetorias1 = dividir_trajetorias(sol)
trajetorias2 = dividir_trajetorias(sol2)

# Plotar as trajetórias divididas
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for traj in trajetorias1:
    traj = np.array(traj)
    ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], label='Trajetória')

for traj in trajetorias2:
    traj = np.array(traj)
    ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], label='Trajetória2')

ax.set_xlabel('u(s)')
ax.set_ylabel('v(s)')
ax.set_zlabel('c(s)')
ax.legend()

# Adicionar o triângulo que vai de c=0 a c=1
vertices = np.array([
    [0, 0, 0], [1, 0, 0], [0, 1, 0],
    [0, 0, 1], [1, 0, 1], [0, 1, 1]
])
edges = [
    (vertices[0], vertices[1]),
    (vertices[0], vertices[2]),
    (vertices[0], vertices[3]),
    (vertices[1], vertices[2]),
    (vertices[1], vertices[4]),
    (vertices[2], vertices[5]),
    (vertices[3], vertices[4]),
    (vertices[3], vertices[5]),
    (vertices[4], vertices[5])
]

for edge in edges:
    ax.plot(*zip(*edge), color='black')
ax.scatter(u0, v0, c0, color='red', s=10, label='Ponto Inicial', edgecolor='black')

plt.show()

