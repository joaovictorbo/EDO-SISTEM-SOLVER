# newton_plot.py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def plotar_trajetorias(trajetorias):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Plota as trajetórias que estão dentro do triângulo
    for traj in trajetorias:
        traj = np.array(traj)
        ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], label="Trajetória Dentro do Triângulo")
        
    ax.set_xlabel("u(s)")
    ax.set_ylabel("v(s)")
    ax.set_zlabel("c(s)")
    ax.legend()
    
    # Adicionar o triângulo que vai de c=0 a c=1
    vertices = np.array([
        [0, 0, 0], [1, 0, 0], [0, 1, 0],
        [0, 0, 1], [1, 0, 1], [0, 1, 1]
    ])
    edges = [
        (vertices[0], vertices[1]), (vertices[0], vertices[2]), (vertices[0], vertices[3]),
        (vertices[1], vertices[2]), (vertices[1], vertices[4]), (vertices[2], vertices[5]),
        (vertices[3], vertices[4]), (vertices[3], vertices[5]), (vertices[4], vertices[5])
    ]
    for edge in edges:
        ax.plot(*zip(*edge), color='black')
    
    plt.show()