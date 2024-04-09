import numpy as np
from scipy.integrate import odeint
import plotly.graph_objects as go
import Inicia as ini

# Concentracoes iniciais
cmin, cmax = ini.concentrations()

# Saturacoes inciais a esquerda no plano cl
# ul, vl, cl = ini.Ul()

# Saturacoes iniciais a direita no plano cr
# ur, vr, cr = ini.Ur()
"""
# Viscosidades
muwl = ini.muw(cl)
muwr = ini.muw(cr)

"""
muw0, muo, mug = ini.viscosidades()

# Define the viscosity functions
def muw(c):
    # Viscosity of water
    return muw0 + c

def muwc(c):
    # dmuw/dc
    return 1.0

# Define the denominator function
def D(u, v, c):
    # Denominator
    w = 1 - u - v
    return u**2 / muw(c) + v**2 / muo + w**2 / mug

# Define the partial derivatives of D
def Du(u, v, c):
    # dD/du
    w = 1 - u - v
    return 2 * u / muw(c) - 2 * w / mug

def Dv(u, v, c):
    # dD/dv
    w = 1 - u - v
    return 2 * v / muo - 2 * w / mug

def Dc(u, v, c):
    # dD/dc
    return -muwc(c) * u**2 / muw(c) ** 2
# Define the functions f, g, and a
def f(u, v, c):
    return u + v + c


def g(u, v, c):
    return u - v - c


def a(c):
    return c + 0.1


def fw(u, v, c):
    # fw
    return (u ** 2 / muw(c)) / D(u, v, c)


def fwu(u, v, c):
    # dfw/du
    return (2 * u / muw(c) * D(u, v, c) - u ** 2 / muw(c) * Du(u, v, c)) / (
        D(u, v, c) ** 2
    )


def fwv(u, v, c):
    # dfw/dv
    return (-(u ** 2) / muw(c) * Dv(u, v, c)) / (D(u, v, c) ** 2)


def fwc(u, v, c):
    # dfw/dc
    return (
        u ** 2
        * (-(muwc(c) / muw(c) ** 2) * D(u, v, c) - Dc(u, v, c) / muw(c))
        / (D(u, v, c) ** 2)
    )


def fo(u, v, c):
    # fo
    return (v ** 2 / muo) / D(u, v, c)


def fou(u, v, c):
    # dfo/du
    return (-(v ** 2) / muo) * Du(u, v, c) / (D(u, v, c) ** 2)


def fov(u, v, c):
    # dfo/dv
    return (2 * v / muo * D(u, v, c) - v ** 2 / muo * Dv(u, v, c)) / (
        D(u, v, c) ** 2
    )


def foc(u, v, c):
    # dfo/dc
    return -(v ** 2) / muo * Dc(u, v, c) / (D(u, v, c) ** 2)


def fg(u, v, c):
    # fg
    w = 1 - u - v
    return (w ** 2 / mug) / D(u, v, c)


def fgu(u, v, c):
    # dfg/du
    w = 1 - u - v
    return (-2 * w / mug * D(u, v, c) - w ** 2 / mug * Du(u, v, c)) / (
        D(u, v, c) ** 2
    )


def fgv(u, v, c):
    # dfg/dv
    w = 1 - u - v
    return (-2 * w / mug * D(u, v, c) - w ** 2 / mug * Dv(u, v, c)) / (
        D(u, v, c) ** 2
    )


def fgc(u, v, c):
    # dfg/dv
    w = 1 - u - v
    return (-(w ** 2) / mug * Dc(u, v, c)) / (D(u, v, c) ** 2)

# Define the right-hand side of the system of ODEs
def model(y, t, f, g, a, fR, ou, gR, ov, cR, sigma):
    u, v, c = y
    du_dt = f(u, v, c) - ou - (fR - ou) * (u - uR)
    dv_dt = g(u, v, c) - sigma * v - (gR - ov) * (v - vR)
    dc_dt = sigma * (a(c) - a(cR)) + (fR - ou) * (c - cR)
    return [du_dt, dv_dt, dc_dt]

# Set the parameter values
fR = 1
ou = 1
gR = 1
ov = 1
cR = 1
sigma = 1
uR = 1
vR = 1

# Set the initial conditions
u0 = 1
v0 = 1
c0 = 1

y0 = [u0, v0, c0] # insert initial values for u, v, and c here

T = 10 # insert value for T here
N = 1000 # insert value for N here
# Set the time points at which to solve the system
t = np.linspace(0, T, N) # insert values for T and N here

# Solve the system of ODEs
sol = odeint(model, y0, t, args=(f, g, a, fR, ou, gR, ov, cR, sigma))

# Extract the solutions for u, v, and c
u = sol[:,0]
v = sol[:,1]
c = sol[:,2]
# Plot the solutions
fig = go.Figure()
fig.add_trace(go.Scatter(x=t, y=u, name='u'))
fig.add_trace(go.Scatter(x=t, y=v, name='v'))
fig.add_trace(go.Scatter(x=t, y=c, name='c'))
fig.show()
