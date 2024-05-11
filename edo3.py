import numpy as np
"""
    a: Este parâmetro representa o limite esquerdo do domínio do problema de valor de contorno (BVP). É um valor numérico que especifica o ponto de partida do intervalo onde o BVP é definido.

    b: Este parâmetro representa o limite direito do domínio do BVP. É um valor numérico que especifica o ponto final do intervalo onde o BVP é definido.

    c: Este parâmetro é uma função que representa o coeficiente do termo de derivada de primeira ordem no BVP. É uma função da variável independente x. No código fornecido, é representado como uma função lambda.

    d: Este parâmetro é uma função que representa o coeficiente do termo de derivada de segunda ordem no BVP. Similar a c, é uma função da variável independente x e é representada como uma função lambda.

    alpha: Este parâmetro representa o valor da solução no limite esquerdo (x = a). É um valor numérico que especifica o valor da solução no ponto de partida do intervalo.

    beta: Este parâmetro representa o valor da solução no limite direito (x = b). Similar a alpha, é um valor numérico que especifica o valor da solução no ponto final do intervalo.

    gamma: Este parâmetro é uma função que representa a condição de contorno nos pontos interiores do domínio. É uma função da variável independente x e é usada para especificar uma condição de contorno diferente de Dirichlet ou Neumann. No código fornecido, é representado como uma função lambda.

    N: Este parâmetro representa o número de intervalos usados para discretizar o domínio [a, b]. Ele determina a granularidade da discretização, com um N maior resultando em uma discretização mais fina e potencialmente uma solução mais precisa."""
def solve_bvp(a, b, c, d, alpha, beta, gamma, N):

    h = (b - a) / N
    x = np.linspace(a, b, N+1)
    A = np.zeros((N+1, N+1))
    b = np.zeros(N+1)

    # Dirichlet boundary conditions
    A[0, 0] = 1
    b[0] = alpha

    A[N, N] = 1
    b[N] = beta

    # Interior points
    for i in range(1, N):
        A[i, i-1] = 1 / h**2 - c(x[i]) / (2 * h)
        A[i, i] = -2 / h**2 + d(x[i])
        A[i, i+1] = 1 / h**2 + c(x[i]) / (2 * h)
        b[i] = -gamma(x[i])

    # Neumann boundary conditions
    A[0, 1] = -1 / h
    b[0] -= alpha / h

    A[N, N-1] = 1 / h
    b[N] += beta / h
    
    # Robin boundary conditions
    A[0, 1] = -1 / h + alpha
    A[0, 0] = 1
    b[0] = 0
    
    A[N, N-1] = 1 / h
    A[N, N] = -1
    b[N] = beta
    

    u = solver(A, b)

    return x, u
def solver(A, b):
    N = len(b)
    # Forward elimination
    for i in range(1, N):
        m = A[i, i-1] / A[i-1, i-1]
        A[i, i] -= m * A[i-1, i]
        b[i] -= m * b[i-1]
    
    # Back substitution
    u = np.zeros(N)
    u[N-1] = b[N-1] / A[N-1, N-1]
    for i in range(N-2, -1, -1):
        u[i] = (b[i] - A[i, i+1] * u[i+1]) / A[i, i]
    
    return u

# Example call
a = 10
b = 20
c = lambda x: 1
d = lambda x: 1
alpha = 1
beta = 1
gamma = lambda x: 0
N = 10

x, u = solve_bvp(a, b, c, d, alpha, beta, gamma, N)

# Print the solution
print("x:", x)
print("u:", u)
import matplotlib.pyplot as plt

# Plot the solution
plt.plot(x, u)
plt.xlabel('x')
plt.ylabel('u')
plt.title('Solution')
plt.grid(True)
plt.show()