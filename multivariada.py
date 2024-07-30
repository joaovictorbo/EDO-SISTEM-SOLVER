import numpy as np

# Definindo a matriz X
X = np.array([
    [1, 2, -1, 3],
    [2, 4, 1, 1],
    [1, 3, -2, 0],
    [1, 4, -1, 7],
    [1, -1, 3, 3],
    [1, 3, 2, 1]
])

# Calculando a matriz de covariância
covariance_matrix = np.cov(X, rowvar=False)
print(covariance_matrix)
# Calculando os autovalores e autovetores da matriz de covariância
eigenvalues, eigenvectors = np.linalg.eig(covariance_matrix)
print("Autovalores:")
print(eigenvalues)
print("Autovetores:")
print(eigenvectors)

# x vezes a matriz de autovetores
print("X vezes a matriz de autovetores:")
xautovetoresY = np.dot(X, eigenvectors)
print(xautovetoresY)
