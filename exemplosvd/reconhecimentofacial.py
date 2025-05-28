# Exemplo 1: Reconhecimento facial via SVD (Eigenfaces)

import numpy as np
from sklearn.datasets import fetch_olivetti_faces
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier

# 1. Carrega o dataset Olivetti Faces
faces = fetch_olivetti_faces(shuffle=True, random_state=0)
X = faces.data           # (400 imagens × 4096 pixels)
y = faces.target         # rótulos de pessoas

# 2. Divide em treino e teste
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.25, random_state=42
)

# 3. Centraliza (subtrai a face média do conjunto de treino)
mean_face = X_train.mean(axis=0)
X_train_centered = X_train - mean_face
X_test_centered  = X_test  - mean_face

# 4. Computa a SVD e obtém as k primeiras "eigenfaces"
k = 50
U, S, Vt = np.linalg.svd(X_train_centered, full_matrices=False)
eigenfaces = Vt[:k]      # cada linha é uma eigenface (dimensão 4096)

# 5. Projeta as imagens nos k componentes principais
#    (coordenadas dos coeficientes para cada face)
proj_train = X_train_centered @ eigenfaces.T   # (300 × k)
proj_test  = X_test_centered  @ eigenfaces.T   # (100 × k)

# 6. Usa k-NN simples sobre o espaço reduzido
clf = KNeighborsClassifier(n_neighbors=3)
clf.fit(proj_train, y_train)
accuracy = clf.score(proj_test, y_test)
print(f"Acurácia no conjunto de teste: {accuracy:.2%}")
