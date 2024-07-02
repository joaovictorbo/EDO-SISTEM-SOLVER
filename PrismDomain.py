#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 10:35:22 2022

@author: cido
"""

import Functions as fun
import Inicia as ini

#mapeia para baricentrica ou não
mp = ini.baricentrica()

# Define as concentracoes minima e maxima
cmin, cmax = ini.concentrations() 

# Para a figura tridimensional
ax = ini.ambiente3d()
ax.view_init(elev=30., azim=-130.) #Initial Camera Position
#ax.axis('off')

#Vertices do triangulo
Gw, Go = 0, 0
Ww, Wo = 1, 0
Ow, Oo = 0, 1

# Mapeia para coordenadas baricentricas de mp = 1
G1, G2 = fun.map(Gw, Go, mp)
W1, W2 = fun.map(Ww, Wo, mp)
O1, O2 = fun.map(Ow, Oo, mp)

# Triangulo no plano cmin
ax.plot([G1,W1], [G2, W2], [cmin,cmin], 'k')
ax.plot([G1,O1], [G2, O2], [cmin,cmin], 'k')
ax.plot([W1,O1], [W2, O2], [cmin,cmin], 'k')

#Triangulo no plano c = cmax
ax.plot([G1,W1], [G2, W2], [cmax,cmax], 'k')
ax.plot([G1,O1], [G2, O2], [cmax,cmax], 'k')
ax.plot([W1,O1], [W2, O2], [cmax,cmax], 'k')

#Arestas verticais
ax.plot([G1,G1], [G2,G2], [cmin,cmax], 'k')
ax.plot([W1,W1], [W2,W2], [cmin,cmax], 'k')
ax.plot([O1,O1], [O2,O2], [cmin,cmax], 'k')

# Se quiser identificar os eixos
ax.set_xlabel('$s_w$')
ax.set_ylabel('$s_o$')
ax.set_zlabel('$c$')