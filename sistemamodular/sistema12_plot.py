# sistema12_plot.py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plotar_trajetorias(trajetorias1, trajetorias2, u0, v0, c0):
    """
    Plota as trajetórias divididas e o triângulo.
    
    Parâmetros:
        trajetorias1, trajetorias2: Trajetórias divididas para os intervalos positivo e negativo.
        u0, v0, c0: Condições iniciais para destacar o ponto inicial.
    """
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
