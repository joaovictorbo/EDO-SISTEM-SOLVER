import numpy as np

# Funções para os cálculos do sistema
alpha = 0
muw0 = 1.0  # Viscosidade inicial sem polimero

def muw(c):  # Viscosidade da agua
    return muw0 * 2**c

def muwc(c):  # dmuw/dc
    return np.log(2) * muw0 * 2**c

def muo():
    return 4.0

def mug():
    return 0.25

def D(u, v, c):
    w = 1 - u - v
    return u**2 / muw(c) + v**2 / muo() + w**2 / mug()

def f(u, v, c):
    return (u**2 / muw(c)) / D(u, v, c)

def g(u, v, c):
    return (v**2 / muo()) / D(u, v, c)

def df_du(u, v, c):
    return (2 * u / muw(c) * D(u, v, c) - u**2 / muw(c) * D_du(u, v, c)) / (D(u, v, c)**2)

def df_dv(u, v, c):
    return (-u**2 / muw(c) * D_dv(u, v, c)) / (D(u, v, c)**2)

def df_dc(u, v, c):
    return u**2 * (-muwc(c) / muw(c)**2 * D(u, v, c) - D_dc(u, v, c) / muw(c)) / (D(u, v, c)**2)

def dg_du(u, v, c):
    return (-v**2 / muo() * D_du(u, v, c)) / (D(u, v, c)**2)

def dg_dv(u, v, c):
    return (2 * v / muo() * D(u, v, c) - v**2 / muo() * D_dv(u, v, c)) / (D(u, v, c)**2)

def dg_dc(u, v, c):
    return -v**2 / muo() * D_dc(u, v, c) / (D(u, v, c)**2)

def D_du(u, v, c):
    w = 1 - u - v
    return 2 * u / muw(c) - 2 * w / mug()

def D_dv(u, v, c):
    w = 1 - u - v
    return 2 * v / muo() - 2 * w / mug()

def D_dc(u, v, c):
    return -u**2 * muwc(c) / muw(c)**2

def a(c):
    return np.sin(c)

def da_dc(c):
    return np.cos(c)  

def lambdac(u, v, c):
    return f(u, v, c) / (u + alpha * da_dc(c))

# Função para calcular a matriz Jacobiana
def calcular_jacobiana(u, v, c, f_R, g_R,sigma_alpha, epsilon_1, epsilon_2, alpha):
    # Derivadas parciais
    dfu = df_du(u, v, c)
    dfv = df_dv(u, v, c)
    dfc = df_dc(u, v, c)
    dgu = dg_du(u, v, c)
    dgv = dg_dv(u, v, c)
    dgc = dg_dc(u, v, c)
    
    # Lambda e sigma
    lambda_c = lambdac(u, v, c)
    
    # Jacobiana
    jacobiana = np.array([
        [dfu - sigma_alpha, dfv, dfc],
        [dgu, dgv - sigma_alpha, dgc],
        [0, 0, (epsilon_1 / epsilon_2) * (f_R - sigma_alpha * (u - alpha * da_dc(c)))]
    ])
    
    return jacobiana

# # Exemplo de uso
# u, v, c = 0.5, 0.3, 0.2  # Ponto de teste
# f_R, g_R = f(u,v,c), g(u,v,c)  # Valores de referência
# epsilon_1, epsilon_2 = 1.0, 1.0  # Parâmetros do sistema

# # Calcula a Jacobiana
# #jacobiana = calcular_jacobiana(u, v, c, f_R, g_R, epsilon_1, epsilon_2)

# # Calcula os autovalores
# autovalores = np.linalg.eigvals(jacobiana)

# # Exibe os resultados
# print("Matriz Jacobiana:")
# print(jacobiana)
# print("\nAutovalores:")
# print(autovalores)