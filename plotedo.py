import numpy as np
from scipy.integrate import odeint
import plotly.graph_objects as go

# Define the functions f, g, and a
def f(u, v, c):
    return u+v+c 

def g(u, v, c):
    return u-v-c

def a(c):
    return c+0.1

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
