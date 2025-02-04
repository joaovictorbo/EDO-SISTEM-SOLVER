import sympy
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# 1) Definir constantes e variáveis simbólicas
# ----------------------------------------------------------------------------
u, v, z = sympy.symbols('u v z', real=True)

# Definir parâmetros fixos como simbólicos
u0, v0, f0, g0, z0 = sympy.symbols('u0 v0 f0 g0 z0', real=True)

# Parâmetros numéricos fixos
alpha = 0.01
muw0 = 1.0

# ----------------------------------------------------------------------------
# 2) Definir as funções simbólicas f(u,v,z) e g(u,v,z)
# ----------------------------------------------------------------------------

# muw(c) = muw0 * 2^c
muw = muw0 * 2**(z)
muo = 4.0
mug = 0.25

# Definir w = 1 - u - v
w = 1 - u - v

# D(u,v,z) = u^2/muw + v^2/muo + w^2/mug
D_expr = (u**2)/muw + (v**2)/muo + (w**2)/mug

# f(u,v,z) = (u^2 / muw) / D
f_expr = ((u**2)/muw) / D_expr

# g(u,v,z) = (v^2 / muo) / D
g_expr = ((v**2)/muo) / D_expr

# ----------------------------------------------------------------------------
# 3) Definir F e G com f0 e g0 já calculados
# ----------------------------------------------------------------------------
F_expr = (f_expr - f0)*(v - v0) - (g_expr - g0)*(u - u0)
G_expr = (f_expr - f0)*(
            u0*(z - z0) + alpha*(sympy.sin(z) - sympy.sin(z0))
         ) - f0*(z - z0)*(u - u0)

# ----------------------------------------------------------------------------
# 4) Atribuir valores a u0, v0, z0 e calcular f0, g0
# ----------------------------------------------------------------------------
valores_fixos = {
    u0: 0.2,  # Defina o valor desejado de u0
    v0: 0.3,  # Defina o valor desejado de v0
    z0: 0.5   # Defina o valor desejado de z0
}

# Substituir em f e g para calcular f0 e g0
f0_val = f_expr.subs({u: valores_fixos[u0], v: valores_fixos[v0], z: valores_fixos[z0]})
g0_val = g_expr.subs({u: valores_fixos[u0], v: valores_fixos[v0], z: valores_fixos[z0]})

# Atualizar os valores fixos com f0 e g0
valores_fixos[f0] = f0_val
valores_fixos[g0] = g0_val

# ----------------------------------------------------------------------------
# 5) Resolver F=0 e G=0 com os valores atribuídos
# ----------------------------------------------------------------------------
F_num = F_expr.subs(valores_fixos)
G_num = G_expr.subs(valores_fixos)

sol = sympy.solve((F_num, G_num), (u, v, z), dict=True)

print("Soluções numéricas para F=0 e G=0 com valores definidos:")
print(sol)

# ----------------------------------------------------------------------------
# 6) Parametrizar e resolver numericamente
# ----------------------------------------------------------------------------
t = sympy.Symbol('t', real=True)

# Substituir u -> t
F_sub = F_num.subs(u, t)
G_sub = G_num.subs(u, t)

# Resolver para (v, z) em função de t
sol_param = sympy.solve((F_sub, G_sub), (v, z), dict=True)

print("\nSoluções paramétricas v(t), z(t) onde u=t:")
print(sol_param)

# ----------------------------------------------------------------------------
# 7) Plotar as soluções paramétricas
# ----------------------------------------------------------------------------
if sol_param:
    sol_family = sol_param[0]
    v_of_t = sol_family[v]
    z_of_t = sol_family[z]

    # Transformar em funções NumPy
    v_func = sympy.lambdify(t, v_of_t, 'numpy')
    z_func = sympy.lambdify(t, z_of_t, 'numpy')

    # Definir range para t
    t_vals = np.linspace(-1.0, 1.0, 200)

    # Gerar pontos
    U_data = []
    V_data = []
    Z_data = []

    for tv in t_vals:
        try:
            vv = v_func(tv)
            zz = z_func(tv)
            if np.isreal(vv) and np.isreal(zz):
                U_data.append(tv)
                V_data.append(float(vv))
                Z_data.append(float(zz))
        except:
            pass

    # Converter para arrays
    U_data = np.array(U_data)
    V_data = np.array(V_data)
    Z_data = np.array(Z_data)

    # Plot 3D
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(U_data, V_data, Z_data, 'r-', label='Curva F=0, G=0')
    ax.set_xlabel('u')
    ax.set_ylabel('v')
    ax.set_zlabel('z')
    ax.set_title('Interseção das superfícies F=0 e G=0')
    ax.legend()
    plt.show()
else:
    print("Não foram encontradas soluções paramétricas para (v,z) em função de t.")
