# Exemplo 2: Sistema de recomendação via SVD

import numpy as np
from scipy.sparse.linalg import svds

# 1. Exemplo de matriz de avaliações (0 = sem avaliação)
#    Linhas = usuários, Colunas = itens  
R = np.array([
    [5, 3, 0, 1],
    [4, 0, 0, 1],
    [1, 1, 0, 5],
    [1, 0, 0, 4],
    [0, 1, 5, 4]
], dtype=float)

# 2. Subtrai a média de cada usuário (para lidar com viés)
user_means = np.true_divide(R.sum(axis=1), (R != 0).sum(axis=1))
R_centered = R - user_means.reshape(-1, 1)
R_centered[R == 0] = 0   # mantém zeros onde não havia avaliação

# 3. Computa SVD truncada (k componentes)
k = 2
U, sigma, Vt = svds(R_centered, k=k)
Sigma = np.diag(sigma)

# 4. Reconstrói a previsão aproximada
R_pred = U @ Sigma @ Vt + user_means.reshape(-1, 1)

# 5. Gera recomendações para um usuário (ex.: usuário 0)
user_id = 0
rated = R[user_id] > 0
preds = R_pred[user_id]
recommendations = [
    (item_id, preds[item_id])
    for item_id in range(R.shape[1])
    if not rated[item_id]
]
recommendations.sort(key=lambda x: x[1], reverse=True)

print("Recomendações para o usuário 0 (item, nota prevista):")
for item, score in recommendations:
    print(f"  Item {item}: {score:.2f}")
