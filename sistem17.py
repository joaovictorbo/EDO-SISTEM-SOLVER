import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

muw0 = 1.0  # Viscosidade inicial sem polimero
muo = 4.0
mug = 0.25

def Du(u, v, c): # dD/du
    w = 1 - u - v
    return 2 * u / muw(c) - 2 * w / mug

def Dv(u, v, c): # dD/dv
    w = 1 - u - v
    return 2 * v / muo - 2 * w / mug

def Dz(u, v, c): # dD/dc
    return -muwc(c) * u**2 / muw(c)**2

def muw(c): # Viscosidade da água
    return muw0 + c

def muwc(c): # dmuw/dc
    return 1.0

def D(u, v, c): # Denominador
    w = 1 - u - v
    return u**2 / muw(c) + v**2 / muo + w**2 / mug

# Definição das funções F e G
def F(u, v, z):
    return (u**2 / muw(z)) / D(u, v, z)

def G(u, v, z):
    return (v**2 / muo) / D(u, v, z)

# Derivadas parciais de F e G
def dF_du(u, v, z):
    return (2 * u / muw(z) * D(u, v, z) - u**2 / muw(z) * Du(u, v, z)) / (D(u, v, z)**2)

def dF_dv(u, v, z):
    return (-u**2 / muw(z) * Dv(u, v, z)) / (D(u, v, z)**2)

def dF_dz(u, v, z):
    return u**2 * (-(muwc(z) / muw(z)**2) * D(u, v, z) - Dz(u, v, z) / muw(z)) / (D(u, v, z)**2)

def dG_du(u, v, z):
    return (-v**2 / muo) * Du(u, v, z) / (D(u, v, z)**2)

def dG_dv(u, v, z):
    return (2 * v / muo * D(u, v, z) - v**2 / muo * Dv(u, v, z)) / (D(u, v, z)**2)

def dG_dz(u, v, z):
    return -v**2 / muo * Dz(u, v, z) / (D(u, v, z)**2)

# Definição do sistema de EDOs
def system(s, y):
    u, v, z = y
    try:
        # Calculando os determinantes
        mat1 = np.array([[dF_dv(u, v, z), dF_dz(u, v, z)],
                         [dG_dv(u, v, z), dG_dz(u, v, z)]])
        mat2 = np.array([[dF_du(u, v, z), dF_dz(u, v, z)],
                         [dG_du(u, v, z), dG_dz(u, v, z)]])
        mat3 = np.array([[dF_du(u, v, z), dF_dv(u, v, z)],
                         [dG_du(u, v, z), dG_dv(u, v, z)]])
        
        if np.isnan(mat1).any() or np.isnan(mat2).any() or np.isnan(mat3).any():
            raise ValueError("Matrix contains NaN values.")
        
        det1 = np.linalg.det(mat1)
        det2 = np.linalg.det(mat2)
        det3 = np.linalg.det(mat3)
        
        if np.isclose(det2, 0):
            raise ZeroDivisionError("det2 is zero.")
        if np.isclose(det3, 0):
            raise ZeroDivisionError("det3 is zero.")
        
        du_ds = det1 / det2
        dv_ds = -det2 / det3
        dz_ds = 1
    except ZeroDivisionError as e:
        print(f"Exception: {e}")
        du_ds, dv_ds, dz_ds = np.nan, np.nan, np.nan
    except ValueError as e:
        print(f"Exception: {e}")
        du_ds, dv_ds, dz_ds = np.nan, np.nan, np.nan
    
    return [du_ds, dv_ds, dz_ds]

# Adicionando o caso de transição
def system_with_transition(s, y):
    u, v, z = y
    try:
        # Calculando os determinantes
        mat_uv = np.array([[dF_dv(u, v, z), dF_dz(u, v, z)],
                           [dG_dv(u, v, z), dG_dz(u, v, z)]])
        mat_uz = np.array([[dF_du(u, v, z), dF_dz(u, v, z)],
                           [dG_du(u, v, z), dG_dz(u, v, z)]])
        mat_vz = np.array([[dF_dv(u, v, z), dF_dz(u, v, z)],
                           [dG_dv(u, v, z), dG_dz(u, v, z)]])
        
        if np.isnan(mat_uv).any() or np.isnan(mat_uz).any() or np.isnan(mat_vz).any():
            raise ValueError(y)
        
        det_uv = np.linalg.det(mat_uv)
        det_uz = np.linalg.det(mat_uz)
        det_vz = np.linalg.det(mat_vz)
        
        if np.isclose(det_uv, 0):
            # Transição para o sistema em (20)
            dv_ds = 1
            if np.isclose(det_vz, 0):
                du_dv, dz_dv = np.nan, np.nan
            else:
                du_dv = -det_uz / det_vz
                dz_dv = det_vz / det_uz
            du_ds = du_dv * dv_ds
            dz_ds = dz_dv * dv_ds
        elif np.isclose(det_uz, 0):
            # Transição para o sistema em (22)
            if np.isclose(det_vz, 0):
                dv_du, dz_du = np.nan, np.nan
            else:
                dv_du = -det_vz / det_uz
                dz_du = det_uz / det_vz
            dv_ds = dv_du * du_ds
            dz_ds = dz_du * du_ds
            du_ds = 1
        else:
            # Sistema original
            du_ds = det_uv / det_uz
            dv_ds = -det_uz / det_vz
            dz_ds = 1
    except ZeroDivisionError as e:
        print(f"Exception: {e}")
        du_ds, dv_ds, dz_ds = np.nan, np.nan, np.nan
    except ValueError as e:
        print(f"Exception: {e}")
        du_ds, dv_ds, dz_ds = np.nan, np.nan, np.nan
    
    return [du_ds, dv_ds, dz_ds]

# Condições iniciais
u0 = 0.5  # Exemplo
v0 = 0.5  # Exemplo
z0 = 0.4  # Exemplo
y0 = [u0, v0, z0]
print(y0)
# Intervalo de integração
s_span = (0, 1)

# Resolução do sistema
sol = solve_ivp(system_with_transition, s_span, y0, method='RK45', t_eval=np.linspace(0, 1, 10))

# Verificação dos resultados
if sol.y.shape[0] != 3:
    raise ValueError("A solução não retornou os valores esperados.")

# Plot dos resultados em 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(sol.y[0], sol.y[1], sol.y[2], label='Trajetória')
ax.set_xlabel('u(s)')
ax.set_ylabel('v(s)')
ax.set_zlabel('z(s)')
ax.legend()
plt.title('Solução do sistema de EDOs em 3D')
plt.show()
