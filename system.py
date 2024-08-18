import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

alpha = 10**-3
muw0 = 1.0  # Viscosidade inicial sem polimero

# Funções dadas
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
    return np.sqrt(c)

def da_dc(c):
    return 1.0 / (2*np.sqrt(c))    

def lambdac(u,v,c):
    return f(u, v, c)/ (u + alpha* da_dc(c))

def detalpha(u, v, c):
    return ((df_du(u, v, c)-lambdac(u,v,c)) * (dg_dv(u, v, c)-lambdac(u,v,c)) - df_dv(u,v,c) * dg_du(u, v, c))

def det1(u, v, c):
    return (df_dv(u, v, c) * dg_dc(u, v, c) - df_dc(u,v,c) * (dg_dv(u, v, c)-lambdac(u,v,c)))

def det2(u,v,c):
    return (df_dc(u,v,c) * dg_du(u,v,c) - (df_du(u,v,c)-lambdac(u,v,c)) * dg_dc(u,v,c))

def p(u, v, c):
    return det1(u, v, c)

def q(u, v, c):
    return det2(u, v, c)

def r(u, v, c):
    return detalpha(u, v, c)

def N(u, v, c):
    return np.sqrt(p(u, v, c)**2 + q(u, v, c)**2 + r(u, v, c)**2)
def system(s, y):
    u, v, c = y
    eq1 = p(u, v, c) / N(u, v, c) if N(u, v, c) else None
    eq2 = q(u, v, c) / N(u, v, c) if N(u, v, c) else None
    eq3 = r(u,v,c) / N(u, v, c) if N(u, v, c) else None
    
    if eq1 is None or eq2 is None or eq3 is None:
        raise ValueError("Ponto singular do sistema de EDOs")
    return [eq1, eq2, eq3]