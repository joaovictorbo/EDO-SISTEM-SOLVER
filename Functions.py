# -*- coding: utf-8 -*-
"""
Spyder Editor

Definicao das principais funcoes dependentes do modelo

fronteira(u,v): testa se o ponto (u, v) estah fora da regiao de interesse

map(x, y, mp): mapeamento para o  triangulo equilatero
umap(x, y, mp): des-mapeamento do triangulo equilatero

muw(c): Viscosidade da agua

D(u,v,c): Denominador
Du(u,v,c): dD/du
Dv(u,v,c): dD/dv
Dc(u,v,c): dD/dc

fw(u,v,c): funcao fluxo agua
fwu(u,v,c): dfw/du
fwv(u,v,c): dfw/dv
fwc(u,v,c): dfw/dc

fw(u,v,c): funcao fluxo oleo
fou(u,v,c): dfo/du
fov(u,v,c): dfo/dv
foc(u,v,c): dfo/dc

fw(u,v,c): funcao fluxo gas
fou(u,v,c): dfg/du
fov(u,v,c): dfg/dv
foc(u,v,c): dfg/dc

jac2(u,v,c): Jacobiana nas saturacoes
jac3(u,v,c): Jacobiana 3x3

lbdc(u,v,c): lambdac (autovalor contato)
lbdas(u,v,c): lambdas (autovalow slow)
lbdaf(u,v,c): lambdaf (autovalow fast)
eigvls(u,v,c): ambos autovalores (slow e fast)

H(u,v,u0,v0,c0): Hugoniot no plano c = c0 por (u0,v0)
dHdu(u,v,u0,v0,c0): dH/du
dHdv(u,v,u0,v0,c0): dH/dv

efe(x,y,z): componente 1 do campo de contato
ge(x,y,z): componente 2 do campo de contato
aga(x,y,z): componente 3 do campo de contato

singularcurve(cmin, cmax, N): curva de singularidades do contato em 3d
Psingular(c): ponto na curva de singularidaes no nivel c
"""

import numpy as np  #funcoes matemáticas
import Inicia as ini

# Concentracoes iniciais
cmin, cmax = ini.concentrations()

# Saturacoes inciais a esquerda no plano cl
#ul, vl, cl = ini.Ul()

# Saturacoes iniciais a direita no plano cr
#ur, vr, cr = ini.Ur()
'''
# Viscosidades
muwl = ini.muw(cl)
muwr = ini.muw(cr)

'''
muw0, muo, mug = ini.viscosidades()

# Testa se o ponto (u,v) estah fora de um quadrado interior ao triangulo
def fronteira(u,v):
    umin, umax, vmin, vmax = ini.dominio()
    w = 1-u-v
    if u < umin or u > umax:
        return True
    elif v < vmin or v > vmax:
        return True
    elif w < 0.0 or w > 1.0:
        return True
    else:
        return False
    
#Mapeamento
def map(x, y, mp): #Funcao que faz o mapeamento do triangulo equilatero
    if mp == 1:
        x = x + 0.5*y
        y = np.sqrt(3)/2*y
    return x,y

#Desmapeamento
def umap(x, y, mp): #Funcao que faz o des-mapeamento do triangulo equilatero
    if mp == 1:
        x = x - 1/np.sqrt(3)*y
        y = 2/np.sqrt(3)*y
    return x,y

#Viscosidade da agua em funcao da concetracao do polimero   
def muw(c): #Viscosidade da agua
    return muw0 + c

def muwc(c): #dmuw/dc
    return 1.0

#Denominador (mobilidade total)
def D(u,v,c): #Denominador
    w = 1 - u - v
    return u**2/muw(c) + v**2/muo + w**2/mug

def Du(u,v,c): #dD/du
    w = 1 - u - v
    return 2*u/muw(c) - 2*w/mug

def Dv(u,v,c): #dD/dv
    w = 1 - u - v
    return 2*v/muo - 2*w/mug

def Dc(u,v,c): ##dD/dc
    return -muwc(c)*u**2/muw(c)**2

#Funcoes de Fluxo

def fw(u,v,c): #fw
    return (u**2/muw(c))/D(u,v,c)

def fwu(u,v,c): #dfw/du
    return (2*u/muw(c)*D(u,v,c) - u**2/muw(c)*Du(u,v,c)) / (D(u,v,c)**2)

def fwv(u,v,c): #dfw/dv
    return (- u**2/muw(c)*Dv(u,v,c)) / (D(u,v,c)**2)

def fwc(u,v,c): #dfw/dc
    return u**2*(- (muwc(c)/muw(c)**2)*D(u,v,c) - Dc(u,v,c)/muw(c)) / (D(u,v,c)**2)
                        
def fo(u,v,c): #fo
    return (v**2/muo)/D(u,v,c)

def fou(u,v,c): #dfo/du
    return (-v**2/muo)*Du(u,v,c) / (D(u,v,c)**2)

def fov(u,v,c): #dfo/dv
    return (2*v/muo*D(u,v,c) - v**2/muo*Dv(u,v,c)) / (D(u,v,c)**2)

def foc(u,v,c): #dfo/dc
    return -v**2/muo*Dc(u,v,c) / (D(u,v,c)**2)

def fg(u,v,c): #fg
    w = 1 - u - v
    return (w**2/mug)/D(u,v,c)

def fgu(u,v,c): #dfg/du
    w = 1 - u - v
    return (-2*w/mug*D(u,v,c)-w**2/mug*Du(u,v,c))/(D(u,v,c)**2)

def fgv(u,v,c): #dfg/dv
    w = 1 - u - v
    return (-2*w/mug*D(u,v,c)-w**2/mug*Dv(u,v,c))/(D(u,v,c)**2)

def fgc(u,v,c): #dfg/dv
    w = 1 - u - v
    return (-w**2/mug*Dc(u,v,c))/(D(u,v,c)**2)


###############################################
#Matriz jacobiana no plano c
def jac2(u,v,c):
    return [[fwu(u,v,c), fwv(u,v,c)], [fou(u,v,c), fov(u,v,c)]]

#Matriz jacobiana no espaco
def jac3(u,v,c):
    fwoversw = (u/muw(c))/D(u,v,c)
    return [[fwu(u,v,c), fwv(u,v,c), fwc(u,v,c)], [fou(u,v,c), fov(u,v,c), foc(u,v,c)],
             [0.0, 0.0, fwoversw]]
 
   
################################################
# Autovalor contato lambdac
def lbdc(u,v,c): #lambdac = fw/sw
    return (u/muw(c))/D(u,v,c) # fw(u,v,c)/u

################################################
#Autovalores nos planos
    
 
# Os dois autovalores separados
def lbdas(u,v,c):
    #Jacobiana de A
    a11 = fwu(u,v,c)
    a12 = fwv(u,v,c)
    a21 = fou(u,v,c)
    a22 = fov(u,v,c)
    
    vneg = 0.5*(a11 + a22 - np.sqrt((a22 - a11)**2 + 4*a21*a12))
    vpos = 0.5*(a11 + a22 + np.sqrt((a22 - a11)**2 + 4*a21*a12))
       
    if vneg.any() <= vpos.any(): # Usar o any, porque vai ser chaamada num vetor
        lambdas = vneg
    else:
        lambdas = vpos
          
    return lambdas

def lbdaf(u,v,c):
    #Jacobiana de A
    a11 = fwu(u,v,c)
    a12 = fwv(u,v,c)
    a21 = fou(u,v,c)
    a22 = fov(u,v,c)
    
    vneg = 0.5*(a11 + a22 - np.sqrt((a11 - a22)**2 + 4*a21*a12))
    vpos = 0.5*(a11 + a22 + np.sqrt((a11 - a22)**2 + 4*a21*a12))
    
    if vneg.any() <= vpos.any():
        lambdaf = vpos
    else:
        lambdaf = vneg
    
    return lambdaf

#Os dois autovalores juntos
def eigvls(u,v,c):
    #Jacobiana de A
    a11 = fwu(u,v,c)
    a12 = fwv(u,v,c)
    a21 = fou(u,v,c)
    a22 = fov(u,v,c)
    
    vneg = 0.5*(a11 + a22 - np.sqrt((a22 - a11)**2 + 4*a21*a12))
    vpos = 0.5*(a11 + a22 + np.sqrt((a22 - a11)**2 + 4*a21*a12))
    
    lambdas = vneg
    lambdaf = vpos
    
    if vneg >= vpos:
        lambdas = vpos
        lambdaf = vneg
    
    return lambdas, lambdaf

#Hugoniot por (u0,v0,c0) no plano c = c0
def H(u,v,u0,v0,c0):
    return (fw(u,v,c0) - fw(u0,v0,c0))*(v-v0) - (fo(u,v,c0) - fo(u0,v0,c0))*(u-u0)

#dH/du por (u0,v0,c0) no plano c = c0
def dHdu(u,v,u0,v0,c0):
    return fwu(u,v,c0)*(v-v0) - fou(u,v,c0)*(u-u0) - (fo(u,v,c0) - fo(u0,v0,c0))

#dH/dv por (u0,v0,c0) no plano c = c0
def dHdv(u,v,u0,v0,c0):
    return fwv(u,v,c0)*(v-v0) - fov(u,v,c0)*(u-u0) + (fw(u,v,c0) - fw(u0,v0,c0))

#Campo de contato

def efe(x,y,z):
    return (lbdc(x,y,z) - fov(x,y,z))*fwc(x,y,z) + fwv(x,y,z)*foc(x,y,z)

def ge(x,y,z):
    return (lbdc(x,y,z) - fwu(x,y,z))*foc(x,y,z) + fou(x,y,z)*fwc(x,y,z)

def aga(x,y,z):
    return (fwu(x,y,z) - lbdc(x,y,z))*(fov(x,y,z)- lbdc(x,y,z)) - fwv(x,y,z)*fou(x,y,z)

# Curva de singularidades do campo de contato   
def singularcurve(cmin, cmax, N):
# Vetores para a curva de pontos singulares do campo de contato
    singw = []
    singo = []
    c = np.linspace(cmin, cmax, N)

# Curva de pontos singulares
    for n in range(0,N):
        muwsing = muw(c[n])
        singw.append(2*muwsing/(muo + 2*muwsing + mug))
        singo.append(muo/(muo + 2*muwsing + mug))

    return(singw,singo,c)
    
# # Ponto Singular no nivel c   
def Psingular(c):
    muwsing = muw(c)
    singw = 2*muwsing/(muo + 2*muwsing + mug)
    singo = muo/(muo + 2*muwsing + mug)

    return(singw,singo)