# -*- coding: utf-8 -*-
"""
Spyder Editor

Desenha o triangulo de saturacoes independente do nivel c
"""
# 
 
#import numpy as np  #funcoes matemáticas
#import sympy as sym #Funcoes simbolicas
import matplotlib.pyplot as plt #Para plotagem de graficos
#from mpl_toolkits.mplot3d import axes3d

import Functions as fun
import Inicia as ini


#if using a Jupyter notebook, include:
#matplotlib inline

#Dados de entrada

#mapeia para baricentrica ou não
mp = ini.baricentrica()

#Vertices do triangulo
Gw, Go = 0, 0
Ww, Wo = 1, 0
Ow, Oo = 0, 1

G1, G2 = fun.map(Gw, Go, mp)
W1, W2 = fun.map(Ww, Wo, mp)
O1, O2 = fun.map(Ow, Oo, mp)

##Ponto umbilico no plano cl
## Precisa inicializar as viscosidades
#Uwl, Uol = muwl/(muwl + muo + mug), muwl/(muwl + muo + mug)
#Uwl, Uol = fun.map(Uwl, Uol, mp)


####################################################

# O triangulo
# Prepara uma janela grafica
plt.figure('uv projection on saturation triangle', figsize=(8,6), dpi=80)
plt.xlabel('$u$')
plt.ylabel('$v$')
#plt.axis('off')
#plt.grid()
plt.plot([G1, W1], [G2, W2], 'k-') 
plt.plot([G1, O1], [G2, O2], 'k-') 
plt.plot([W1, O1], [W2, O2], 'k-') 

# O quadrado
#plt.plot([0, 1], [0, 0], 'k-') 
#plt.plot([0, 0], [0, 1], 'k-') 
#plt.plot([1, 1], [0, 1], 'k-') 
#plt.plot([1, 0], [1, 1], 'k-')

# Plota o grid
#plt.grid() 
