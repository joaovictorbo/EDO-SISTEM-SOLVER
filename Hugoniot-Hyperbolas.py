# -*- coding: utf-8 -*-
"""
Hugoniots ao longo dos lados e das bifurcacoes secundarias
Formulas explicitas das dissertacoes da Joseane, Luciano e Patricio


Tentativa de melhorar a partir do arquivo .save

Verificar porque o ramo nao local não estah funcionado para v0 < 0.56
Fazer a intersecao do ramo com os lados mais preciso.

"""
import numpy as np  #funcoes matemáticas
#import sympy as sym #Funcoes simbolicas
import matplotlib.pyplot as plt #Para plotagem de graficos
#from mpl_toolkits.mplot3d import axes3d

import Functions as fun
import Inicia as ini
import FuncoesAuxiliares as fa


 
# Initialization

#Baricentric or retangular coordinates
mp = ini.baricentrica()

#Set the level
c0 = 0.0

# Domain
umin, umax, vmin, vmax, domain = ini.dominio()

#Reference viscosities
muw0, muo, mug = ini.viscosidades()

# update the water viscosity
muw0 = fun.muw(c0)

#Step h for sg in [0, 1] to pararetrize the hyperbolas branches
h = 0.00001#0.001
# Maximal number of steps along the hyperbola branches
N = 100000 #100000

# Umbilic point
mutotal = muw0 + muo + mug
u_umb = muw0 / mutotal
v_umb = muo / mutotal


## Input data to line EW
## U0 along the segment  EW: so/sg = muo/mug. 
## E = (0, muo/[muo + mug], mug/[muo + mug])
# w0 maximo = mug/(muo + mug)
#print('v0 maximo = ', muo/(muo + mug))
#v0 = muo/(muo + mug)#0.56
v0 = 0.5
w0 = (mug / muo) * v0
u0 = 1.0 - v0 - w0
print('u0=', u0, 'v0=', v0, 'w0=', w0,'\n')


# Speeds at U0
lbs0 = fun.lbdas(u0, v0, c0)
lbf0 = fun.lbdaf(u0, v0, c0)
lbc0 = fun.lbdc(u0, v0, c0)


## Implicit hyperbola equation for U0 along E-W
## (i) A*so^2  + B*sg^2 + C*so*sg  - D*so - E*sg + F =0, 
## Explicit hyperbola equations for U0 along E-W: so = so(sg)
## so = 1/2A * [(D - C*sg) +- sqrt((D - C*sg)^2 - 4*A*(B*sg^2 -E*sg + F)]
  
## Coeficientes
D0 = u0**2/muw0 + v0**2/muo + w0**2/mug

A = (muw0 + muo) / (muw0*muo) * v0**2 #coef de so^2 em (i)
B = (muw0 + mug) / (muw0*mug) * v0**2 #coef de sg*so em (i)
C = (muo/mug) * D0 + 2*v0**2 / muw0 # coef de sg^2 em (i)
D = (2*v0 / muw0 + D0) * v0 #coef de so em (i)
E = (2*v0 / muw0 + D0 * muo / mug) * v0 #coef de sg em (i)
F = v0**2 / muw0 #constante em (i)

# To control the sign  of sqrt(Delta) and to determin which is the local branch
theta = 1 # The local branch corresponds to the minus sign

# To store the branch with the sign "-" in the explict formulas
ulocal = []
vlocal = [] 
lbdaslocal = []
lbdaflocal = []
lbdaclocal = []
sigmalocal = []

# To store the branch with the sign "+" in the explict formulas
unlocal = []
vnlocal = []
lbdasnlocal = []
lbdafnlocal = []
lbdacnlocal = []
sigmanlocal = []

# To store the branch E-W of the hyperbola given by so/muo = sg/mug
vline = []
lbdasline = []
lbdafline = []
lbdacline = []
sigmaline = []


for k in range(2):       
    auxh = D - C*w0
    Delta = auxh**2 - 4*A*(B*w0**2 - E*w0 + F)
    
    if Delta >= 0.0:
        vaux1 = 1.0/(2*A)*(auxh - theta * np.sqrt(Delta))
        uaux1 = 1.0 - vaux1- w0       
        dist1 = abs(u0 - uaux1) + abs(v0 - vaux1)
        
        vaux2 = 1.0/(2*A)*(auxh + theta * np.sqrt(Delta))
        uaux2 = 1.0 - vaux2 - w0       
        dist2 = abs(u0 - uaux2) + abs(v0 - vaux2)
  
        if k == 0 and dist1 > dist2:
            uaux3 = uaux2
            vaux3 = vaux2
            uaux2 = uaux1
            vaux2 = vaux1
            uaux1 = uaux3
            vaux1 = vaux3            

            theta = - theta # The local branch corresponds to the plus sign
           
        for i in range(0,N):            
            # choose interior points
            interior1 = fa.triangletest(uaux1, vaux1, umin, umax, vmin, vmax, domain)
            if interior1 == 1:
                lbds = fun.lbdas(uaux1, vaux1, c0)
                lbdf = fun.lbdaf(uaux1, vaux1, c0)
                lbc = fun.lbdc(uaux1, vaux1, c0)
                sigma = fa.sigma(uaux1, vaux1, u0, v0, c0)

                ulocal.append(uaux1) 
                vlocal.append(vaux1)
                
                lbdaslocal.append(lbds)
                lbdaflocal.append(lbdf)
                lbdaclocal.append(lbc)
                sigmalocal.append(sigma)
                
            w1 = 1.0 - uaux1 - vaux1 + h
            auxh = D - C*w1
            Delta = auxh**2 - 4*A*(B*w1**2 - E*w1 + F) 

            if Delta >= 0.0:
                vaux1 = 1.0/(2*A)*(auxh - theta * np.sqrt(Delta))
                uaux1 = 1.0 - vaux1- w1
            
            interior2 = fa.triangletest(uaux2, vaux2, umin, umax, vmin, vmax, domain)
            if interior2 == 1:
                lbds = fun.lbdas(uaux2, vaux1, c0)
                lbdf = fun.lbdaf(uaux2, vaux1, c0)
                lbc = fun.lbdc(uaux2, vaux1, c0)
                sigma = fa.sigma(uaux2, vaux2, u0, v0, c0)

                
                unlocal.append(uaux2)
                vnlocal.append(vaux2)
                
                lbdasnlocal.append(lbds)
                lbdafnlocal.append(lbdf)
                lbdacnlocal.append(lbc)
                sigmanlocal.append(sigma)
                
            w2 = 1.0 - uaux2 - vaux2 + h
            auxh = D - C*w2
            Delta = auxh**2 - 4*A*(B*w2**2 - E*w2 + F) 
    
            if Delta >= 0.0:
                vaux2 = 1.0/(2*A)*(auxh + theta * np.sqrt(Delta))
                uaux2 = 1.0 - vaux2 - w2

# Reversing the sense of integration in each branch              
    h = -h
    if k == 0:
        #store the first part of the local branch
        u11 = ulocal
        v11 = vlocal
        
        lbdas11 = lbdaslocal
        lbdaf11 = lbdaflocal
        lbdac11 = lbdaclocal
        sigma11 = sigmalocal
        del(sigma11[0]) # Remove sigma(U0; U0)
        
        ulocal = [] # to store the second part of the local branch
        vlocal = []
        lbdaslocal = []
        lbdaflocal = []
        lbdaclocal = []
        sigmalocal = []

        #store first part of the detached branch
        u21 = unlocal
        v21 = vnlocal
        
        lbdas21 = lbdasnlocal
        lbdaf21 = lbdafnlocal
        lbdac21 = lbdacnlocal
        sigma21 = sigmanlocal
        
        unlocal = [] # to store the second part of the detached branch
        vnlocal = []
        lbdasnlocal = []
        lbdafnlocal = []
        lbdacnlocal = []
        sigmanlocal = []
        
    else: # k = 1
        u12 = ulocal
        v12 = vlocal
        
        lbdas12 = lbdaslocal
        lbdaf12 = lbdaflocal
        lbdac12 = lbdaclocal
        sigma12 = sigmalocal
        del(sigma12[0]) # Remove sigma(U0; U0)

        u22 = unlocal
        v22 = vnlocal
        
        lbdas22 = lbdasnlocal
        lbdaf22 = lbdafnlocal
        lbdac22 = lbdacnlocal
        sigma22 = sigmanlocal

# Points along the line EW: so/sg = muo / mug, or so = muo / (mu + mug)*sw
Nline = 100
uline = np.linspace(0, 1, Nline)
muo_over_muog = muo / (muo + mug) 
for k in range(Nline):
    vl = muo_over_muog * (1.0 - uline[k])
    
    lbds = fun.lbdas(uline[k], vl, c0)
    lbdf = fun.lbdaf(uline[k], vl, c0)
    lbc = fun.lbdc(uline[k], vl, c0)
    if (abs(uline[k] - u0) + abs(vl - v0)) >= 1 / Nline**2: # Testa se eh o ponto base U_0
        sigma = fa.sigma(uline[k], vl, u0,v0, c0)
    else:
        if u0 > u_umb:
            sigma = lbf0
        else:
            sigma = lbs0

    vline.append(vl)
    
    lbdasline.append(lbds)  
    lbdafline.append(lbdf)
    lbdacline.append(lbc)
    sigmaline.append(sigma)

N11 = len(u11)
for j11 in range(N11):
    u11[j11], v11[j11] = fun.map(u11[j11], v11[j11], mp)
    
N12 = len(ulocal)
for j12 in range(N12):
    u12[j12], v12[j12] = fun.map(u12[j12], v12[j12], mp)
    
N21 = len(u21)
for j21 in range(N21):
    u21[j21], v21[j21] = fun.map(u21[j21], v21[j21], mp)
    
N22 = len(u22)
for j22 in range(N22):
    u22[j22], v22[j22] = fun.map(u22[j22], v22[j22], mp)
    
#print(Nline, 'interior points in the line EW')
for jr in range(Nline):
    uline[jr], vline[jr] = fun.map(uline[jr], vline[jr], mp)


plt.figure('uv projection on saturation triangle')
u0, v0 = fun.map(u0, v0, mp)
plt.plot(u0, v0,'ko') # U0
u_umb, v_umb = fun.map(u_umb, v_umb, mp)
plt.plot(u_umb, v_umb,'go') # Umbilic
# Ordena os arrays e remove o primeiro item antes de concatenar
u11, v11 = sorted(u11[1:]), sorted(v11[1:])
u12, v12 = sorted(u12[1:]), sorted(v12[1:])

# Concatena os arrays
u_combined1 = u11 + u21
v_combined1 = v11 + v21
u_combined2 = u12 + u22
v_combined2 = v12 + v22

# Plot resultados
plt.figure('uv projection on saturation triangle')
u0, v0 = fun.map(u0, v0, mp)
plt.plot(u0, v0, 'ko')  # U0
u_umb, v_umb = fun.map(u_umb, v_umb, mp)
plt.plot(u_umb, v_umb, 'go')  # Umbilic

plt.plot(u_combined1, v_combined1, 'r-')  # Combined lower branches
plt.plot(u_combined2, v_combined2, 'y-')  # Combined upper branches
plt.plot(uline, vline, 'r-')
plt.show()


## Speeds along the local branch-11
plt.figure('speeds along branch-11')
# Map the indexes to the interval [0, 1]
Ns11 = len(lbdas11)
indexlocal = np.linspace(0, 1, Ns11)

plt.plot(indexlocal, lbdas11,'b-')
plt.plot(indexlocal, lbdaf11,'r-')

Nsigma11 = len(sigma11)
indexlocal = np.linspace(0, 1, Nsigma11)
plt.plot(indexlocal, sigma11,'r--')

## Speeds along the local branch
plt.figure('speeds along branch-12')
# Map the indexes to the interval [0, 1]
Ns12 = len(lbdas12)
indexlocal = np.linspace(0, 1, Ns12)

plt.plot(indexlocal, lbdas12,'b-')
plt.plot(indexlocal, lbdaf12,'r-')

Nsigma12 = len(sigma12)
indexlocal = np.linspace(0, 1, Nsigma12)
plt.plot(indexlocal, sigma12,'y--')

## Speeds along the nonlocal branch
plt.figure('speeds along branch-21')
# Map the indexes to the interval [0, 1]
Ns21 = len(lbdas21)
indexlocal = np.linspace(0, 1, Ns21)

plt.plot(indexlocal, lbdas21,'b-')
plt.plot(indexlocal, lbdaf21,'r-')

#Nsigma12 = len(sigma12)
#indexlocal = np.linspace(0, 1, Nsigma12)
plt.plot(indexlocal, sigma21,'b--')

## Speeds along the nonlocal branch
plt.figure('speeds along branch-22')
# Map the indexes to the interval [0, 1]
Ns22 = len(lbdas22)
indexlocal = np.linspace(0, 1, Ns22)

plt.plot(indexlocal, lbdas22,'b-')
plt.plot(indexlocal, lbdaf22,'r-')

#Nsigma12 = len(sigma12)
#indexlocal = np.linspace(0, 1, Nsigma12)
plt.plot(indexlocal, sigma22,'g--')

## Speeds along the line EW
plt.figure('speeds along EW')
# Map the indexes to the interval [0, 1]
indexlocal = np.linspace(0, 1, Nline)

plt.plot(indexlocal, lbdasline,'b-')
plt.plot(indexlocal, lbdafline,'r-')

#Nsigma12 = len(sigma12)
#indexlocal = np.linspace(0, 1, Nsigma12)
plt.plot(indexlocal, sigmaline,'r--')






