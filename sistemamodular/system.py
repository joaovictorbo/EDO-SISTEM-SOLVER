import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

alpha = 0.0#0.001
muw0 = 1.0  # Viscosidade inicial sem polímero

# Funções dadas
def muw(c):  # Viscosidade da água
    return muw0 * 2**c

def muo():
    return 4.0

mug = 0.25

# Função de adsorção
def a(z):
    return np.sin(z)

def da_dz(z):
    return np.cos(z)

# Derivadas de D
def Du(u, v, c):  # dD/du
    w = 1 - u - v
    return 2 * u / muw(c) - 2 * w / mug

def Dv(u, v, c):  # dD/dv
    w = 1 - u - v
    return 2 * v / muo() - 2 * w / mug

def Dz(u, v, c):  # dD/dc
    return -muwc(c) * u**2 / muw(c)**2

def muwc(c):  # dmuw/dc
    return np.log(2) * muw0 * 2**c

def D(u, v, c):
    w = 1 - u - v
    return u**2 / muw(c) + v**2 / muo() + w**2 / mug

# Definição das funções f e g
def f(u, v, z):
    return (u**2 / muw(z)) / D(u, v, z)

def g(u, v, z):
    return (v**2 / muo()) / D(u, v, z)

# Derivadas parciais de f
def df_du(u, v, z):
    return (2 * u / muw(z) * D(u, v, z) - u**2 / muw(z) * Du(u, v, z)) / (D(u, v, z)**2)

def df_dv(u, v, z):
    return (-u**2 / muw(z) * Dv(u, v, z)) / (D(u, v, z)**2)

def df_dz(u, v, z):
    return u**2 * (-(muwc(z) / muw(z)**2) * D(u, v, z) - Dz(u, v, z) / muw(z)) / (D(u, v, z)**2)

# Derivadas parciais de g
def dg_du(u, v, z):
    return (-v**2 / muo()) * Du(u, v, z) / (D(u, v, z)**2)

def dg_dv(u, v, z):
    return (2 * v / muo() * D(u, v, z) - v**2 / muo() * Dv(u, v, z)) / (D(u, v, z)**2)

def dg_dz(u, v, z):
    return -v**2 / muo() * Dz(u, v, z) / (D(u, v, z)**2)

# Funções F e G do sistema de Rankine-Hugoniot
def F(u, v, z, f0, v0, g0, u0): # Eq. 29(a) da versao 2025-04-02
    f_value = f(u, v, z)
    g_value = g(u, v, z)
    return (f_value - f0) * (v - v0) - (g_value - g0) * (u - u0)

def G(u, v, z, f0, z0, u0, a, alpha): # Eq. 29(b) da versao 2025-04-02
    f_value = f(u, v, z)
    return (f_value - f0) * (u0 * (z - z0) + alpha * (a(z) - a(z0))) - f0 * (z - z0) * (u - u0)

# Derivadas parciais de F e G
def dF_du(u, v, z, f0, v0, g0, u0): #Eq (30)  da versao 2025-04-02
    return df_du(u, v, z) * (v - v0) - dg_du(u, v, z) * (u - u0) - (g(u, v, z) - g0)

def dF_dv(u, v, z, f0, v0, g0, u0): #Eq (31)  da versao 2025-04-02
    return df_dv(u, v, z) * (v - v0) - dg_dv(u, v, z) * (u - u0) + (f(u, v, z) - f0)

def dF_dz(u, v, z, f0, v0, g0, u0): #Eq (32)  da versao 2025-04-02
    return df_dz(u, v, z) * (v - v0) - dg_dz(u, v, z) * (u - u0)

def dG_du(u, v, z, f0, z0, u0, alpha): #Eq (33)  da versao 2025-04-02
    return df_du(u, v, z) * (u0 * (z - z0) + alpha * (a(z) - a(z0))) - f0 * (z - z0)

def dG_dv(u, v, z, f0, z0, u0, alpha): #Eq (34)  da versao 2025-04-02
    return df_dv(u, v, z) * (u0 * (z - z0) + alpha * (a(z) - a(z0)))

def dG_dz(u, v, z, f0, z0, u0, alpha): #Eq (35)  da versao 2025-04-02
    return (df_dz(u, v, z) * (u0 * (z - z0) + alpha * (a(z) - a(z0))) +
            (f(u, v, z) - f0) * (u0 + alpha * da_dz(z)) - f0 * (u - u0))

# Outras funções auxiliares relacionadas a D
def D_du(u, v, c):
    w = 1 - u - v
    return 2 * u / muw(c) - 2 * w / mug

def D_dv(u, v, c):
    w = 1 - u - v
    return 2 * v / muo() - 2 * w / mug

def D_dc(u, v, c):
    return -u**2 * muwc(c) / muw(c)**2

# Versão duplicada da função a (para variável c)
def a_c(c):
    return np.sin(c)

def da_dc(c):
    return np.cos(c)

def lambdac(u, v, c):
    return f(u, v, c) / (u + alpha * da_dc(c))

def detalpha(u, v, c): # D_r da versao 2025-04-02
    return ((df_du(u, v, c) - lambdac(u, v, c)) *
            (dg_dv(u, v, c) - lambdac(u, v, c)) - df_dv(u, v, c) * dg_du(u, v, c))

def det1(u, v, c): #D_p da versao 2025-04-02
    return (df_dv(u, v, c) * dg_dz(u, v, c) -
            df_dz(u, v, c) * (dg_dv(u, v, c) - lambdac(u, v, c)))

def det2(u, v, c): #D_q da versao 2025-04-02
    return (df_dz(u, v, c) * dg_du(u, v, c) -
            (df_du(u, v, c) - lambdac(u, v, c)) * dg_dz(u, v, c))

def p(u, v, c):
    return det1(u, v, c)

def q(u, v, c):
    return det2(u, v, c)

def r(u, v, c):
    return detalpha(u, v, c)

def N(u, v, c):
    return np.sqrt(p(u, v, c)**2 + q(u, v, c)**2 + r(u, v, c)**2)

def system(s, y): # Sistema (19) da versao 2025-04-02
    u, v, c = y
    eq1 = p(u, v, c) / N(u, v, c) if N(u, v, c) else None
    eq2 = q(u, v, c) / N(u, v, c) if N(u, v, c) else None
    eq3 = r(u, v, c) / N(u, v, c) if N(u, v, c) else None
    if eq1 is None or eq2 is None or eq3 is None:
        raise ValueError("Ponto singular do sistema de EDOs")
    return np.array([eq1, eq2, eq3])