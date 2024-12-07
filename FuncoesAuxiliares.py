'''


Cuidado. Ver modificacoews em A-ProgramacaoPython. 06/01/24
Rever DiflbdclbdcZero


def prismdomain(): Desenha o Prisma e coloca os nomes dos eixos

def Diflbdc(u0,v0,c0, u, v): Diferenca lambdac(u, v,c0) por lambdac(u0,v0,c0)

Difcs(u,v,c) que devolve a diferença lambda^c(u,v,c)- lambda^s(u,v,c)

Difcf(u,v,c) que devolve a diferença lambda^c(u,v,c) - lambda^f(u,v,c)

Difcsigma(u,v,u0,v0,c0) que devolve a diferença lambda^c(u,v,c)- sigma((u,v,c0);(u0,v0,c0))

Difsigmaconst(u,v,u0,v0,c0,cnt) que devolve a diferença #calcula sigma((u,v,co);(u0,v0,c0)). cnt é qlq constante

lineBD(u0, v0, c0):  #Substitui (u,v) na equacao da reta BD e retorna o valor

Hugoniot(u0,v0,c0) que devolve os vetores u e v

sigma(u,v,u0,v0,c0) que devolve a velocidade de choque sigma

sigmacontact(u,v,c, u0,v0,c0) que devolve a velocidade de choque sigma ao longo d eum ramo de contato

coinc(fam, c0) calcula os pontos de T^s se fam=1 e de T^f se fam = 2

BoundaryPoints(u0,v0,c0): devolve quantos e quais pontos de H(U0) estao na fronteira

IntContact(u0, v0, c0, h, cmin, cmax, N):
    calcula o ramo de contato por um ponto U0=(u0, v0, c0)
    a partir de U0 com passo h via o Metodo de Euler
    Devolve os vetores u,v,c

SingleContatBranch(u0,v0,c0, cmin, cmax, updow, coinc): devolve os vetores u, v, c
                     de pontos comecando do nivel c0 para o nivel cmin, se updow = -1,
                     ou para o nivel cmax se updow = 1.
                     Se coinc = 1, o ramo para qdo atingir coincidencia

LinBissection(A, B, Ndim, f, TOL): Aproxima um zero da funcao f no segmento AB
                                   com Tolerancia TOL                                 

LinBissectionsigma(A, B, u0, v0, c0, cnt, f, TOL): Aproxima um zero da funcao f 
                   no segmento AB com Tolerancia TOL. Usar com f do tipo sigma - cnt

ExplCoinc(fam, c): fam = 1 =s low ou 2fam = 2 =f ast. c: nivel. 
    Ideia: Usar as formulas explictas das coincidencias (extensoes da fronteira GO)
    do paper Azevedo et al, SIAP-2014

Cextremos(c):    
    #Determinacao dos primeiros indices apos um valor extremo (minimo ou maximo)
    # em um vetor de pontos c.
        #Retorna a quantidade de mudancas no crescimento e os indices correspondentes
        
Coincidencias(a,b,c):
    Testa as diferencas lambdac - lambdas e lambdac - lambdaf no ponto (a,b,c)
    Devolve 1 se for lambdac = lambdas e 2 se lambdac = lambdaf
    
triangletest(u, v):
    Define se (u,v) pertence ao dominio de interesse interior ao triangulo de saturacoes
    
def Elipse(u0, v0, c0, N):
    # Funcao que obtem os pontos da elipse -curva de nivel lamndac = lambdac(U0) explicitamente
    # Devolve os vetores (up, vp):
    # N numero de subdivisoes do intervalo [umin umax], sendo umin e umax onde du/dv = 0
    
LambdacEqualSigma(u0, v0, c0, u, v):
    testa as diferencas lambdac - sigma ao longo do segmento de ramo de 
    Hugoniot (u, v), para c = constante
    para determinar pontos com sigma*((u0,v0); (u,v)) = lambdac(u0, v0, c0)

HugoniotContactIntersection(u0, v0, c0, u, v):
    determina pontos na interseccao de uma curva de Hugoniot por (u0, v0, c0) com 
    a elipse lambdac(u,v) = lambdac(u0,v0,c0) calculada por LambdacContourExpl(u0,v0,c0, N).
    Devolve os arrays u, v onde podem acontecer tais interseccoes
    
def cutoffTsTf(u,v,c0,coinc):
    # 2023-11-23
    # Cutoff the segment along the curve (u,v) at the level c = c0 before coincidences
    # u and v are arrays with the same length
    # if coinc = 1, lambda^c = lambda^s
    # if coinc = 2, lambda^c = lambda^f
    # return (unew, vnew)
    
def ball(a, b, c, r):
    # 2023-12-27
    # (a, b, c) center
    # r radius

'''
   
import numpy as np  #funcoes matemáticas
import math as mat #funcoes matemáticas

from scipy.integrate import solve_ivp

import Functions as fun
import Inicia as ini
    
# Define o dominio de calculos
umin, umax, vmin, vmax, triang = ini.dominio()
cmin, cmax = ini.concentrations()

#mapeia para baricentrica ou não
mp = ini.baricentrica()

def prismdomain():
  
    #mapeia para baricentrica ou não
    #mp = ini.baricentrica()
    
    # Define as concentracoes minima e maxima
    #cmin, cmax = ini.concentrations() 
        
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
    ax.set_xlabel('$u$')
    ax.set_ylabel('$v$')
    ax.set_zlabel('$z$')
    
          
    return(ax)

def Diflbdc(u0,v0,c0, u, v): # Melhorar para nao repetir
    #print('Diflbdc: dif0U =', fun.lbdc(u0,v0,c0), '\n')
    #print('Diflbdc: difU=', fun.lbdc(u,v,c0), '\n')
    dif = fun.lbdc(u0,v0,c0) - fun.lbdc(u,v,c0)
    #print('Diflbdc: dif =', dif, '\n')
    return(dif)

def Difcs(u,v,c):
    dif = fun.lbdc(u,v,c) - fun.lbdas(u,v,c)
    return(dif)
    
def Difcf(u,v,c):
    dif = fun.lbdc(u,v,c) - fun.lbdaf(u,v,c)
    return(dif)
    
def Difcsigma(u,v,u0,v0,c0):
    dif = fun.lbdc(u,v,c0) - sigma(u,v,u0,v0,c0)
    return(dif)


def lineBD(u0, v0, c0):  #Substtui (u,v) na equacao da reta BD
    muw0, muo, mug = ini.viscosidades()
    muw = fun.muw(c0)
    valor = v0 - (muo/muw)*((muw + mug)/(mug-muo))*u0 - muo/(muo-mug)
    return(valor)
  
def Hugoniot(u0,v0,c0):
     
    # Define os passos em u e em v
    du = 0.001
    dv = 0.01
    
    # Define a tolerancia para considerar uma igualdade
    tol = 10**(-5)
     
    # Para guardar os pontos da curva
    vectoru = []
    vectorv = []
    #vectorw = []
    
    v = vmin
    
    while v <=vmax:
        u1 = umin
        u2 = u1 + du
        #print('Print10: u1=', u1, 'u2=', u2,'\n')
        aux1 = fun.H(u1,v,u0,v0,c0)
        aux2 = fun.H(u2,v,u0,v0,c0)
        #print('Print20: aux1=', aux1, 'aux2=', aux2,'\n')

        while u2 <= umax and u2 + v <= 1.0:
            if aux1 == 0.0 and np.abs(aux2) > 0:            
                vectoru.append(u1)
                vectorv.append(v)
                #vectorw.append(1-u1-v)
                #print('Print 30: aux1 = 0, guardou (u1, v) \n')
                u1 = u2
                u2 = u1 + du
                aux1 = aux2
                aux2 = fun.H(u2,v,u0,v0,c0)
            elif aux2 == 0.0 and np.abs(aux1) > 0:           
                vectoru.append(u2)
                vectorv.append(v)
                #vectorw.append(1-u1-v)
                #print('Print 31: aux2 = 0, guardou (u2, v) \n')
                u1 = u2 + du
                u2 = u1 + du
                aux1 = fun.H(u1,v,u0,v0,c0)
                aux2 = fun.H(u2,v,u0,v0,c0)
            elif aux1 == 0.0 and aux2 == 0:           
                vectoru.append(u1)
                vectorv.append(v)
                #vectorw.append(1-u1-v)
                vectoru.append(u2)
                vectorv.append(v)
                #vectorw.append(1-u2-v)
                #print('Print 32: aux1 = aux2 = 0, guardou (u1, v) e (u2, v) \n')
                u1 = u2 + du
                u2 = u1 + du
                aux1 = fun.H(u1,v,u0,v0,c0)
                aux2 = fun.H(u2,v,u0,v0,c0)
            elif aux1*aux2 < 0:
                #print('Print 90: aux1*aux2 < 0. Há raiz. Bissecao \n')
                umedio = (u1+u2)/2.0
                auxmedio = fun.H(umedio,v,u0,v0,c0)
                u2save = u2
                #print('Print 100: u1=', u1, 'u2=', u2, 'umedio=', umedio, 'auxmedio=', auxmedio, '\n')
                while np.abs(auxmedio) >= tol:
                    #print('Print 101: u1=', u1, 'u2=', u2, 'auxmedio =', auxmedio,'\n')
                    if auxmedio*aux1 < 0:
                        u2 = umedio
                    else:
                        u1 = umedio
                    #print('Print 102: u1=', u1, 'u2=', u2, 'auxmedio =', auxmedio,'\n')
    
                    umedio = (u1 + u2)/2.0
                    auxmedio = fun.H(umedio,v,u0,v0,c0)
                    #print('Print 103: umedio=', umedio, 'auxmedio=', auxmedio, '\n')
                vectoru.append(umedio)
                vectorv.append(v)
                #vectorw.append(1-umedio-v)
                u2 = u2save
                #print('Print 105: guardou u=umedio=', umedio, '\n')
                u1 = u2
                u2 = u1 + du
                aux1 = aux2
                aux2 = fun.H(u2,v,u0,v0,c0)
            else:
                u1 = u2
                u2 = u1 + du
                aux1 = aux2
                aux2 = fun.H(u2,v,u0,v0,c0)
                #print('Print 106: aux1*aux2 > 0 \n')
                
            #print('Print 110: volta no loop em u, u1=', u1, 'u2=', u2, '\n',
             #     'aux1=', aux1, 'aux2=', aux2, '\n')
        v = v + dv
        #print('Print 120: volta no loop em v, v=', v, '\n')
        #print('v=',v,'\n')                  
            
    #return(vectoru, vectorv, vectorw)
    return(vectoru, vectorv)
'''  
#Para testar
u0 = 0.78441#0.699999
v0 = 0.199111#0.266669
c0 = 0.5
u, v= Hugoniot(u0, v0, c0)
#print('u      v  \n')
#for n in range(0,len(u)):
    #print(np.round(u[n],4), np.round(v[n],4), '\n')

import matplotlib.pyplot as plt #Para plotagem de graficos    
plt.figure('saturation triangle')
plt.plot(u, v, 'r.')
'''

# Funcao que calcula a velocidade sigma na Hugoniot de (u0,v0,c0)
def sigma(u,v,u0,v0,c0):
    denomu = u - u0
    numeru = fun.fw(u,v,c0) - fun.fw(u0,v0,c0)
    denomv = v - v0
    numerv = fun.fo(u,v,c0) - fun.fo(u0,v0,c0)
    if np.abs(denomu) > 0.0:
        valor = numeru/denomu
    elif np.abs(denomv) > 0.0:
        valor = numerv/denomv
    else:
        valor = -1.0 # sigma = 0, se R = L Consertar isto. Eliminar este else?
        print('\n Function sigma: U_R = U_L. Artifical sigma = ', -1, ' \n ')
    return valor
'''
#Para testar sigma, lambdas, lambdaf e lambdac ao longo dos ramos
u0 = 0.78441#0.699999
v0 = 0.199111#0.266669
c0 = 0.0
u, v = Hugoniot(u0, v0, c0)
N = len(u)
sigmas = []
lambdass = []
lambdasf = []
lambdasc = []
for n in range(0,N):
    sigmas.append(sigma(u[n],v[n], u0, v0, c0))
    lambdass.append(fun.lbdas(u[n],v[n], c0))
    lambdasf.append(fun.lbdaf(u[n],v[n], c0))
    lambdasc.append(fun.lbdc(u[n],v[n], c0))

#import matplotlib.pyplot as plt #Para plotagem de graficos    
plt.figure('velocidades')
plt.plot(v, sigmas, '.k')
#plt.plot(u, lambdass, '.b')
#plt.plot(u, lambdasf, '.r')
plt.plot(v, lambdasc, '.m')
plt.xlabel('v')
plt.ylabel('lambdac')
plt.grid()
print('labdac(U0)=', fun.lbdc(u0,v0,c0), '\n')
'''

# Funcao que calcula a velocidade sigma no ramo de contato da Hugoniot por (u0,v0,c0)
def sigmacontact(u,v,c, u0,v0,c0):
    denomu = u - u0
    numeru = fun.fw(u,v,c) - fun.fw(u0,v0,c0)
    denomv = v - v0
    numerv = fun.fo(u,v,c) - fun.fo(u0,v0,c0)
    if np.abs(denomu) > 0.0:
        valor = numeru/denomu
    elif np.abs(denomv) > 0.0:
        valor = numerv/denomv
    else:
        valor = -1.0 # sigma = 0, se R = L Consertar isto. Elimiar este else?
        print('\n Function sigma: U_R = U_L. Artifical sigma = ', -1, ' \n ')
    return valor


def coinc(fam,c): #fam: familia 1=slow ou 2=fast. c: nivel
# Usa as diferencas lambdac - lambdas e lambdac - lambdaf
# Faz uma malha no dominio e procura os pontos de coincidencia para muw = muw(c) 
# A funcao ExplCoinci abaixo usa formnulas explicitas
    
    # Define o passo nas direcoes u e v
    deltau = 0.01
    deltav = 0.001
    
    # Define a tolerancia para considerar uma igualdade
    tol = 10**(-5)
    
    # Para guardar os pontos da coincidencia lambdas = lambdac
    vectoru = []
    vectorv = []
    vectorw = []
   
    if fam == 1 and umin == 0 and vmin == 0:
           
        # Guarda o vertice G em Ts
        vectoru.append(0.0)
        vectorv.append(0.0)
        vectorw.append(1.0)
        
        v = deltav # Para iniciar o loop em v jah com v > 0 para a familia 1
                                             
    # Loop em v para pesquisar os pontos em cada reta v constante a partir de v=deltav
    if fam == 2:
        v = vmin # Inicia o loop em v jah em v = vmin para a familia 2
    
    u1 = umin
    u2 = u1 + deltau #Começa u já no interior
    
    sinal = 1.0 # Sabe-se que lambda_s - lambda_c = 0 no lado u = 0.
                 # e que lambda_s - lambda_c > 0 proximo ao lado u = 0.
                 # Sabe-se que lambda_f - lambda_c > no lado u = 0.
    
    v = vmin # Acrescentado em 02/09/23
    while v <= vmax: # Loop vertical em v
        if fam == 1:
            difsc2 = fun.lbdas(u2,v,c) - fun.lbdc(u2,v,c) # Testa se já houve interseccao com Ts
        if fam == 2:
            difsc2 = fun.lbdaf(u2,v,c) - fun.lbdc(u2,v,c) # Testa se já houve interseccao com Tf
        
        if np.abs(difsc2) <= tol:
            vectoru.append(u2)
            vectorv.append(v)
            vectorw.append(1.0-u2-v)
            sinal = -1.0 # Força não entrar no loop em u abaixo
            
        while sinal > 0.0 and u2 + v <= 1.0: # Loop horizontal em u ateh encontrar a interseccao              
            
            if difsc2 < 0.0: # Detectou uma intersecção. Use Bissecao.
                while np.abs(difsc2) > tol:
                    umedio = (u1+u2)/2.0
                    if fam == 1:
                        difsc2 = fun.lbdas(umedio,v,c) - fun.lbdc(umedio,v,c)
                    if fam == 2:
                        difsc2 = fun.lbdaf(umedio,v,c) - fun.lbdc(umedio,v,c)
                    if difsc2 < 0.0:
                        u2 = umedio
                    else:
                        u1 = umedio
        
                vectoru.append(umedio)
                vectorv.append(v)
                vectorw.append(1.0-umedio-v)
                sinal = -1.0
            else:
                u1 = u2    
                u2 = u2 + deltau #Para o while em u
                if fam == 1:
                    difsc2 = fun.lbdas(u2,v,c) - fun.lbdc(u2,v,c)
                if fam == 2:
                    difsc2 = fun.lbdaf(u2,v,c) - fun.lbdc(u2,v,c)
                if np.abs(difsc2) <= tol:
                    vectoru.append(u2)
                    vectorv.append(v)
                    vectorw.append(1.0-u2-v)
                    sinal = -1.0 # Força não entrar no loop em u
        
        v = v + deltav # para op while em v
        u1 = umin
        u2 = u1 + deltau # Pra recomecar no proximo nivel v
        sinal = 1.0
    # Guarda o vertice O em Ts
    if fam == 1 and umin == 0 and vmax == 1:
        vectoru.append(0.0)
        vectorv.append(1.0)
        vectorw.append(0.0)
    
    return(vectoru, vectorv, vectorw)
    
'''
#Para testar
u, v, w = coinc(1,0)
print('u      v      w \n')
for n in range(0,len(u)):
    print(np.round(u[n],4), np.round(v[n],4), np.round(w[n],4), '\n')

import matplotlib.pyplot as plt #Para plotagem de graficos    
plt.figure('saturation triangle')
plt.plot(u, v, 'r')
'''

## Pesquisa de pontos interseccao da Hugoniot de (u0,v0,c0) no lado GW em que v = 0
## Considerar (u0, v0) no interior do triangulo
def BoundaryPoints(u0,v0,c0):
    #Calcula os  pontos de interseccao da Hugoniot H(U0)
    # com o lados do dominio considerado.
    
    du = 0.1 # Passos para pesquisa
    dv = 0.1
    
    tol = 10**(-5) # tolerancia para considerar zero
    pontosu = [] # Para guardar os pontos de interseccao
    pontosv = [] # Para guardar os pontos de interseccao
    
    #muw0, muo, mug = ini.viscosidades()
    #muw0 = fun.muw(c0) # atualiza a viscosidade muw para o nivel c0
    
    fw0 = fun.fw(u0,v0,c0)
    fo0 = fun.fo(u0,v0,c0)
    
    ###########Intersecao com os lado v = vmin e v = vmax #########
    
    for i in range(0,2):
        if i == 0:
            cte = vmin
        if i == 1:
            cte = vmax       
        
        ua = umin
        ub = ua + du
        
        aux0 = u0*fo0 + (cte - v0)*fw0
        
        valora = aux0 - (cte - v0)*fun.fw(ua,cte,c0) +\
                (ua - u0)*fun.fo(ua,cte,c0) - ua*fo0
        
        usup = np.minimum(umax, 1.0 - cte) # Valor maximo de um interior ao triangulo
        while ub >= umin and ub <= usup:
            valorb = aux0 - (cte - v0)*fun.fw(ub,cte,c0) +\
                    (ub - u0)*fun.fo(ub,cte,c0) - ub*fo0        
            if valora == 0.0 and np.abs(valorb) > 0.0:
                pontosu.append(ua)
                pontosv.append(cte)
                valora = valorb
                ua = ub
                ub = ua + du
                valora = valorb
            elif valorb == 0.0 and np.abs(valora) > 0.0:
                pontosu.append(ub)
                pontosv.append(cte)
                ua = ub + du
                ub = ua + du
                valora = aux0 - (cte - v0)*fun.fw(ua,cte,c0) +\
                    (ua - u0)*fun.fo(ua,cte,c0) - ua*fo0
            elif valora == 0.0 and np.abs(valorb) == 0.0:
                pontosu.append(ua)
                pontosv.append(cte)
                pontosu.append(ub)
                pontosv.append(cte)
                ua = ub + du
                ub = ua + du
                valora = aux0 - (cte - v0)*fun.fw(ua,cte,c0) +\
                    (ua - u0)*fun.fo(ua,cte,c0) - ua*fo0
            elif valora*valorb < 0.0: #Aplique bissecao
                ubsave = ub #Salva o ponto a direita
                valorsave = valorb #Salva o ponto a direita
                um = (ua+ub)/2.0
                valorm = aux0 - (cte - v0)*fun.fw(um,cte,c0) +\
                    (um - u0)*fun.fo(um,cte,c0) - um*fo0
                while np.abs(valorm) > tol:
                    if valora*valorm < 0.0:
                        ub = um
                    else:
                        ua = um
                        valora = valorm
                    um = (ua + ub)/2.0
                    valorm = aux0 - (cte - v0)*fun.fw(um,cte,c0) +\
                        (um - u0)*fun.fo(um,cte,c0) - um*fo0
                pontosu.append(um)
                pontosv.append(cte)
                ua = ubsave
                ub = ua + du
                valora = valorsave
            elif valora*valorb > 0.0:
                ua = ub
                ub = ua + du
                valora = valorb

###########Intersecao com os lados u = umin e u = umax #########    
    for j in range(0, 2):
        if j == 0:
            cte = umin

        if j == 1:
            cte = umax

        vsup = np.minimum(vmax, 1.0 - cte) # Valor maximo de um interior ao triangulo

        va = vmin
        vb = va + dv
        
        aux0 = v0*fw0 + (cte - u0)*fo0

        valora = aux0 - (cte - u0)*fun.fo(cte,va,c0) +\
                (va - v0)*fun.fw(cte,va,c0) - va*fw0
        while vb >= vmin and vb <= vsup:
            valorb = aux0 - (cte - u0)*fun.fo(cte,vb,c0) +\
                (vb - v0)*fun.fw(cte,vb,c0) - vb*fw0
            if valora == 0.0 and np.abs(valorb) > 0.0:
                pontosu.append(cte)
                pontosv.append(va)
                valora = valorb
                va = vb
                vb = va + dv
                valora = valorb
            elif valorb == 0.0 and np.abs(valora) > 0.0:
                pontosu.append(cte)
                pontosv.append(vb)
                va = vb + dv
                vb = va + dv
                valora = aux0 - (cte - u0)*fun.fo(cte,va,c0) +\
                    (va - v0)*fun.fw(cte,va,c0) - va*fw0
            elif valora == 0.0 and np.abs(valorb) == 0.0:
                pontosu.append(cte)
                pontosv.append(va)
                pontosu.append(cte)
                pontosv.append(vb)
                va = vb + dv
                vb = va + dv
                valora = aux0 - (cte - u0)*fun.fo(cte,va,c0) +\
                    (va - v0)*fun.fw(cte,va,c0) - va*fw0
            elif valora*valorb < 0.0: #Aplique bissecao
                vbsave = vb #Salva o ponto a direita
                valorsave = valorb #Salva o ponto a direita
                vm = (va+vb)/2.0
                valorm = aux0 - (cte - u0)*fun.fo(cte,vm,c0) +\
                    (vm - v0)*fun.fw(cte,vm,c0) - vm*fw0
                while np.abs(valorm) > tol:
                    if valora*valorm < 0.0:
                        vb = vm
                    else:
                        va = vm
                        valora = valorm
                    vm = (va + vb)/2.0
                    valorm = aux0 - (cte - u0)*fun.fo(cte,vm,c0) +\
                        (vm - v0)*fun.fw(cte,vm,c0) - vm*fw0
                pontosu.append(cte)
                pontosv.append(vm)
                va = vbsave
                vb = va + dv
                valora = valorsave
            elif valora*valorb > 0.0:
                va = vb
                vb = va + dv
                valora = valorb
           
###########Intersecao com o lado u + v = 1 #########
    du = 0.1
    aux0 = u0*fo0 - v0*fw0
    uinf = np.maximum(umin, 1.0 - vmax) 
    usup = np.minimum(umax, 1.0 - vmin)

    ua = uinf
    ub = ua + du
    fwa = fun.fw(ua,1.0-ua,c0)
    foa = fun.fo(ua,1.0-ua,c0)
    valora = aux0 + (ua - u0)*foa -\
            (1.0 - ua)*(fwa - fw0) - ua*fo0 + v0*fwa
    while ub >= uinf and ub <= usup:
        fwb = fun.fw(ub,1.0-ub,c0)
        fob = fun.fo(ub,1.0-ub,c0)
        valorb = aux0 + (ub - u0)*fob -\
            (1.0 - ub)*(fwb - fw0) - ub*fo0 + v0*fwb
        if valora == 0.0 and np.abs(valorb) > 0.0:
            pontosu.append(ua)
            pontosv.append(1.0-ua)
            valora = valorb
            ua = ub
            fwa = fwb
            valora = valorb
            ub = ua + du
        elif valorb == 0.0 and np.abs(valora) > 0.0:
            pontosu.append(ub)
            pontosv.append(1.0-ub)
            ua = ub + du
            fwa = fun.fw(ua,1.0-ua,c0)
            foa = fun.fo(ua,1.0-ua,c0)
            ub = ua + du
            valora = aux0 + (ua - u0)*foa -\
                (1.0 - ua)*(fwa - fw0) - ua*fo0 + v0*fwa
        elif valora == 0.0 and np.abs(valorb) == 0.0:
            pontosu.append(ua)
            pontosv.append(1.0-ua)
            pontosu.append(ub)
            pontosv.append(1.0-ub)
            ua = ub + du
            fwa = fun.fw(ua,1.0-ua,c0)
            foa = fun.fo(ua,1.0-ua,c0)
            ub = ua + du
            valora = aux0 + (ua - u0)*foa -\
                (1.0 - ua)*(fwa - fw0) - ua*fo0 + v0*fwa
        elif valora*valorb < 0.0: #Aplique bissecao
            ubsave = ub #Salva o ponto a direita
            valorsave = valorb #Salva o ponto a direita
            um = (ua+ub)/2.0
            fwm = fun.fw(um,1.0-um,c0)
            fom = fun.fo(um,1.0-um,c0)
            valorm = aux0 + (um - u0)*fom -\
                (1.0 - um)*(fwm - fw0) - um*fo0 + v0*fwm
            while np.abs(valorm) > tol:
                if valora*valorm < 0.0:
                    ub = um
                else:
                    ua = um
                    valora = valorm
                    fwa = fwm
                    foa = fom
                um = (ua + ub)/2.0
                fwm = fun.fw(um,1.0-um,c0)
                fom = fun.fo(um,1.0-um,c0)
                valorm = aux0 + (um - u0)*fom -\
                   (1.0 - um)*(fwm - fw0) - um*fo0 + v0*fwm
            pontosu.append(um)
            pontosv.append(1.0-um)
            ua = ubsave
            ub = ua + du
            valora = valorsave
        elif valora*valorb > 0.0:
            ua = ub
            ub = ua + du
            valora = valorb

    NPGW = np.size(pontosu)

    return (pontosu, pontosv, NPGW) 

# Ramo da Hugoniot por um ponto U0, via EDOs/ Runge-Kutta Ordem 4
def HugBranch(uf, vf, u0, v0, c0):
    
    umin, umax, vmin, vmax = ini.dominio() #fixa o dominio de trabalho
    
    h = 0.001
    Nmax = 10000
    tol = 10**(-2)
    
    uhug = []
    vhug = []
    
    uhug.append(uf)
    vhug.append(vf)
    for n in range(1, Nmax):
        dHdu = fun.dHdu(uhug[n-1], vhug[n-1], u0, v0, c0)
        dHdv = fun.dHdv(uhug[n-1], vhug[n-1], u0, v0, c0)
        
        if np.abs(dHdv) > tol:
            uhug.append(uhug[n-1] + h)
            vhug.append(vhug[n-1] - h*dHdu/dHdv)
        
        elif np.abs(dHdu) > tol:
            uhug.append(uhug[n-1] - h*dHdv/dHdu)
            vhug.append(vhug[n-1] + h)
        
        else:
            print('Singularidade na Hugoniot por (',u0, v0, c0,')',
                  'Usingular = (', round(uhug[n-1],4), round(vhug[n-1],4), c0,')\n')
            break # Interrompe o loop for
            
        if fun.fronteira(uhug[n], vhug[n]) == True: # Se o ponto final estah fora do dominio
            # Faz a correcao para tomar U na fronteira e termina
            if uhug[n] < umin:
                print('fronteira u =', umin, '\n')
                if np.abs(dHdv) > tol:
                    hnew = umin - uhug[n-1]
                    uhug[n] = umin
                    vhug[n] = vhug[n-1] - hnew*dHdu/dHdv
                elif np.abs(dHdu) > tol:
                    hnew = (uhug[n-1] - umin)*dHdu/dHdv
                    uhug[n] = umin
                    vhug[n] = vhug[n-1] + hnew
                        
            if uhug[n] > umax:
                print('fronteira u =', umax, '\n')
                if np.abs(dHdv) > tol:
                    hnew = umax - uhug[n-1]
                    uhug[n] = umax
                    vhug[n] = vhug[n-1] - hnew*dHdu/dHdv
                elif np.abs(dHdu) > tol:
                    hnew = (uhug[n-1] - umax)*dHdu/dHdv
                    uhug[n] = umax
                    vhug[n] = vhug[n-1] + hnew
            
            if vhug[n] < vmin:
                print('fronteira v =', vmin, '\n')
                if np.abs(dHdv) > tol:
                    hnew = (vhug[n-1] - vmin)*dHdv/dHdu
                    uhug[n] = uhug[n-1] + hnew
                    vhug[n] = vmin
                elif np.abs(dHdu) > tol:
                    hnew = vmin - vhug[n-1]
                    uhug[n] = uhug[n-1] - hnew*dHdv/dHdu
                    vhug[n] = vmin
            
            if vhug[n] > vmax:
                print('fronteira v =', vmax, '\n')
                if np.abs(dHdv) > tol:
                    hnew = (vhug[n-1] - vmax)*dHdv/dHdu
                    uhug[n] = uhug[n-1] + hnew
                    vhug[n] = vmax
                elif np.abs(dHdu) > tol:
                    hnew = vmax - vhug[n-1]
                    uhug[n] = uhug[n-1] - hnew*dHdv/dHdu
                    vhug[n] = vmax
                
            break # Interrompe o loop for
        
    NPhugGW = np.size(uhug)
    
    return(uhug, vhug, NPhugGW)
                
'''
#Para testar
# Define por coordenadas baricentricas ou nao
mp = ini.baricentrica()
u0 = 0.076
v0 = 0.18
c0 = 0.5
print('U0 = (',u0, v0, c0,')')
Pontosu, Pontosv , NPGW = BoundaryPoints(u0,v0,c0)
print('Npontos=', NPGW, 'Pontosu=', Pontosu,'Pontosv=', Pontosv, '\n')
import matplotlib.pyplot as plt #Para plotagem de graficos    
plt.figure('saturation triangle') 

u0, v0 = fun.map(u0, v0, mp)
plt.plot(u0, v0, 'ro')
for k in range(0, NPGW):
    Pontosu[k], Pontosv[k] = fun.map(Pontosu[k], Pontosv[k], mp) 

plt.plot(Pontosu, Pontosv, 'bo')
'''

## Funcao que calcula o ramo de contato por um ponto U0=(u0, v0, c0)
## a partir de U0 com passo h pelo Metodo de Euler

def IntContact(u0, v0, c0, h, cmin, cmax, N):
    # U0 = (u0, v0, c0): Ponto inicial
    # h : Passo de integracao
    # cmin: Valor minimo de c fixado
    # cmax: valor maximo de c fixado
    # N: Numero maximo de pontosa serem calculados
        
    u = [] # Para guardar as aproximacoes com h > 0
    v = [] # Para guardar as aproximacoes 
    c = [] # Para guardar as aproximacoes
    
    u.append(u0) # Valor inicial
    v.append(v0)
    c.append(c0)
    
    # Sentido com passo h > 0
    #for n in range(0, N):
    n = 0
    auxc = c[0]
    while auxc >= cmin and auxc <= cmax and n < N:
        n = n+1
        # Integra
        u.append(u[n-1] + h*fun.efe(u[n-1], v[n-1],c[n-1]))
        v.append(v[n-1] + h*fun.ge(u[n-1], v[n-1],c[n-1]))
        c.append(c[n-1] + h*fun.aga(u[n-1], v[n-1],c[n-1]))
        
        auxc = c[n]
        
    hanterior = abs(fun.aga(u[n-1], v[n-1], c[n-1]))
    
    
    if auxc < cmin:
        #hanterior = abs(fun.aga(u[n-1], v[n-1], c[n-1]))
        if hanterior> 10**(-5):
            hnew = (cmin - auxc)/hanterior # Neste caso h > 0
            u[n] = u[n-1] + hnew*fun.efe(u[n-1], v[n-1], c[n-1])
            v[n] = v[n-1] + hnew*fun.ge(u[n-1], v[n-1], c[n-1]) 
            c[n] = cmin
      
    if auxc > cmax:
        #hanterior = abs(fun.aga(u[n-1], v[n-1], c[n-1]))
        if hanterior> 10**(-5):
            hnew = (auxc - cmax)/hanterior # Neste caso h > 0
            u[n] = u[n-1] + hnew*fun.efe(u[n-1], v[n-1], c[n-1])
            v[n] = v[n-1] + hnew*fun.ge(u[n-1], v[n-1], c[n-1]) 
            c[n] = cmax
          
    if n == N:
        print('\n Funcoes Auxiliares/IntContact usou o maximo de pontos N =', N,'com passo h =', h, '\n')
        
    return(u, v, c)
'''
# Para testar        
#Exibe ponts calculados num ramo de contato

cmin = 0.0
cmax = 1.0

u0 = 0.4
v0 = 0.4
c0 = 0.4

h = 0.01
N = 10

u, v, c = IntContact(u0, v0, c0, h, cmin, cmax, N)
print('\n Funcoes Auxiliares: Funcao IntContact. \n u=\n', u, '\n\n v =\n', v, '\n\n c =\n', c,'\n')

'''

## Funcao que calcula os pontos do campo de contato comecando em cmax e
## terminando em cmin, com cmin <= cmax

def SingleContatBranch(u0, v0, c0, cmin, cmax, updow, stopcoinc):
    #cmin deve nao ser superior a cmax e updow indica o
    #sentido de crescimento de c em c0. Se updow = -1 decresce e se updow = 1 cresce
    #Se stopcoinc = 1, entao forca a parada numa superficie de coincidencia    
    
    
    # Para integrar os contatos (por Euler)
    h = 0.01 # Passo
    Nmax = 10000 # Numero maximo de pontos
    
    # Para igualdades
    tol = 10**(-5)
    
    #exibe o ponto inicial U0 na tela
    #print('SingleContatBranch function. u0=', u0, 'v0=', v0, 'c0=',c0)
    
    # Calcula as diferencas do autovalor contato para os outros dois
    difcs = Difcs(u0,v0,c0) #fun.lbdc(u0, v0, c0) - fun.lbdas(u0, v0, c0)
    difcf = Difcf(u0,v0,c0) #fun.lbdc(u0, v0, c0) - fun.lbdaf(u0, v0, c0)
    
    if np.abs(difcs) <= tol:
        print('SingleContatBranch function. U0 sob a coincidencia lbdac = lbdas')
        Nmax = 0 #Calcula nada alem de U0
    if np.abs(difcf) <= tol:
        print('SingleContatBranch function. U0 sob a coincidencia lbdac = lbdaf')
        Nmax = 0 #Calcula nada alem de U0
    
    Npontos = 1 # Numero de pontos calculados. Ha pelo menos o ponto U_0
    
    u = [] # Para guardar as aproximacoes com c decrescente
    v = [] # Para guardar as aproximacoes 
    c = [] # Para guardar as aproximacoes
    
    # Fixa o sentido decrescente para c proximo de c0
    if updow == -1:       
        auxc = c0 + h*fun.aga(u0, v0, c0)
        if auxc > c0: # c crescente. Inverter o sentido de integracao
            h = -h
        auxc = c0 # Redefine caux para entrar no loop a seguir
        
        u.append(u0) # Valor inicial
        v.append(v0)
        c.append(c0)
        
        n = 0 # Para controlar o numero de pontos
        while auxc > cmin and n < Nmax:
            # Integra por Euler
            u.append(u[n] + h*fun.efe(u[n], v[n], c[n]))
            v.append(v[n] + h*fun.ge(u[n], v[n], c[n]))
            c.append(c[n] + h*fun.aga(u[n], v[n], c[n]))
           
            auxc = c[n+1] # Para testar o while
            n = n+1
        
        Npontos = np.size(c) # Fixa o n^o de pontos calculados
        
        if n == Nmax and c[n] > cmin:
            print('SingleContatBranch function. Nmax = ', Nmax,'not enough. \n')
            print('u =', np.round(u[Nmax], 6),
                      "v =", np.round(v[Nmax],6),
                      'c =', np.round(c[Nmax],6), '\n')
        elif c[n] <= cmin: # Faz a correcao no ultimo ponto para que c[n] = cmin
            hnew = (cmin - c[n-1]) / fun.aga(u[n-1], v[n-1],c[n-1])
            u[n] = u[n-1] + hnew*fun.efe(u[n-1], v[n-1], c[n-1])
            v[n] = v[n-1] + hnew*fun.ge(u[n-1], v[n-1], c[n-1])
            c[n] = cmin
        
            #Retorna o ponto de interseccao com o plano c = cmin
            #print("SingleContatBranch function. cminimum reached: u =", np.round(u[Npontos-1], 6),
                     #"v =", np.round(v[Npontos-1],6),
                      #'c =', np.round(c[Npontos-1],6), '\n')
            
    # Fixa o sentido crescente a partir de c0
    if updow == 1:
        auxc = c0 + h*fun.aga(u0, v0,c0) # Para testar o crescimento de c
        ccoinc = 1 # Para o caso de intersectar uma superficie de coincidencia
        
        if auxc < c0: # c decrescente. Inverter o sentido de integracao
            h = -h
          
        auxc = c0 # Redefine caux para entrar no loop a seguir
        
        u.append(u0) # Valores iniciais
        v.append(v0)
        c.append(c0)
        
        n = 0 # Para controlar o numero de pontos
        #lbdb = fun.lbdc(u[0],v[0],c[0])
        auxcsb = Difcs(u[0], v[0], c[0])
        auxcfb = Difcf(u[0], v[0], c[0]) 
        prodb = auxcsb*auxcfb
        #print('00  L779 Single: auxcsb=', auxcsb, 'auxcfb=', auxcfb, 'prodb=',prodb)
        while cmin <= auxc <= cmax and ccoinc == 1 and n < Nmax:
            # Integra por Euler
            #lbda = lbdb
            auxcsa = auxcsb
            auxcfa = auxcfb 
            proda = prodb
            #print('-01 proda=', proda, 'n=', n)
            
            u.append(u[n] + h*fun.efe(u[n], v[n], c[n]))
            v.append(v[n] + h*fun.ge(u[n], v[n], c[n]))
            c.append(c[n] + h*fun.aga(u[n], v[n], c[n]))
            
            #lbdb = fun.lbdc(u[n+1],v[n+1],c[n+1])
            auxcsb = Difcs(u[n+1],v[n+1],c[n+1]) #lbdb - fun.lbdas(u[n+1],v[n+1],c[n+1])
            auxcfb = Difcf(u[n+1],v[n+1],c[n+1]) #lbdb - fun.lbdaf(u[n+1],v[n+1],c[n+1]) 
            prodb = auxcsb*auxcfb
            #print('01 L796: auxcsb =', auxcsb, ' auxcfb=', auxcfb,'n=', n)
            # Para testar o final do while no caso de encontrar
            # uma superficie de coincidencia
            
            if proda*prodb <= 0.0 and stopcoinc ==1: #Deve parar se stopcoinc = 1
                ccoinc = -1 # caso cruze uma superficie de coincidencia
                #print('02 L802 n=', n, 'u[n]=',u[n], 'v[n]=',v[n], 'c[n]=',c[n], 'ccoinc=', ccoinc)
            auxc = c[n+1]
            n = n+1
        #print('06 n=', n)
        Npontos = np.size(c) # Fixa o n^o de pontos calculados
        #print('07 Npontos=', Npontos, 'n=', n)
        if n == Nmax:
            print('SingleContatBranch function. Nmax not enough. Singular point? \n')
            print('u =', np.round(u[n], 6),
                      "v =", np.round(v[n],6),
                      'c =', np.round(c[n],6), '\n')
        
        if auxc > cmax: # Faz a correcao no ultimo ponto para que c[n] = cmax
            n = Npontos - 1 # Porque foi acrescentado um no loop acima
            hnew = (cmax - c[n-1]) / fun.aga(u[n-1], v[n-1],c[n-1])
            u[n] = u[n-1] + hnew*fun.efe(u[n-1], v[n-1], c[n-1])
            v[n] = v[n-1] + hnew*fun.ge(u[n-1], v[n-1], c[n-1])
            c[n] = cmax
            
            #Retorna o ponto de interseccao com o plano c = cmax
            #print("SingleContatBranch function. cmaximum reached: u =", np.round(u[Npontos-1], 6),
                     #"v =", np.round(v[Npontos-1],6),
                      #'c =', np.round(c[Npontos-1],6), '\n')
        #print('L 825 ccoinc=', ccoinc, 'n=', n)    
        if ccoinc == -1:
            if Npontos > 1:
                n = Npontos-1 # Porque parou com o acrescimo n = n+1
                #rint('08 n=', n)
                ua = u[n-1] # Fixa o penultimo ponto calculado (na mesma regiao de U0)
                va = v[n-1]
                ca = c[n-1]
                #lbda = fun.lbdc(ua,va,ca)
                auxcsa = Difcs(ua,va,ca) #lbda - fun.lbdas(ua,va,ca)
                auxcfa = Difcf(ua,va,ca) #lbda - fun.lbdaf(ua,va,ca)          
                proda = auxcsa*auxcfa
                #print('L 837 SingleConhtact Branch auxcfa=', auxcfa, 'proda=',proda)
                    
                ub = u[n] # Fixa o ultimo ponto (na regiao oposta a U0)
                vb = v[n]
                cb = c[n]
                #lbdb = fun.lbdc(ub,vb,cb)
                auxcsb = Difcs(ub, vb, cb) #lbdb - fun.lbdas(ub, vb, cb)
                auxcfb = Difcf(ub, vb, cb) #lbdb - fun.lbdaf(ub, vb, cb)
                prodb = auxcsb*auxcfb
                #print('10 ub=', ub, 'vb=', vb, 'cb=',cb)
                #print('11 auxcsa=', auxcsa, 'auxcfa=', auxcfa)
                #print('L848 Single CB  auxcfb=', auxcfb, 'prodb=', prodb)
                #print('13 proda=', proda)
                if np.abs(auxcsa) < tol or np.abs(auxcfa) < tol:
                    print('SingleContatBranch function. Coincidence surface reached'
                       "u =", np.round(ua, 6),
                       "v =", np.round(va, 6),
                       "c =", np.round(ca, 6), '\n')
                else: # Aplicar Bissecao
                    hnew = h
                    # Para teste
                    p = 0
                    auxcsm = auxcsa
                    auxcfm = auxcfa                                       
                    while np.abs(auxcsm) >= tol and np.abs(auxcfm) >= tol: #Consertei aqui
                        p = p+1
                        #print('L862 SingleConhtact Branch p=', p)
                        #print('L863 SingleConhtact Branch tol =', tol, 'auxcsa=', auxcsa, 'auxcsb=', auxcsb)
                        #if p==20:
                            #break
                        #print('14 ua=', ua, 'va=', va, 'ca=',ca)
                        #print('14.5 ub=', ub, 'vb=', vb, 'cb=',cb)
                        hnew = hnew/2.0
                        #print('14.6 hnew=', hnew)
                        um = ua + hnew*fun.efe(ua, va, ca)
                        vm = va + hnew*fun.ge(ua, va, ca)
                        cm = ca + hnew*fun.aga(ua, va, ca)
                        #lbdm = fun.lbdc(um,vm,cm)
                        auxcsm = Difcs(um, vm, cm) #lbdm - fun.lbdas(um, vm, cm)
                        auxcfm = Difcf(um, vm, cm) #lbdm - fun.lbdaf(um, vm, cm)
                        prodm = auxcsm*auxcfm
                        #print('L877 Single CB auxcfm =', auxcfm)
                        #print('L878 SingleConhtact Branch um=', um, 'vm=', vm, 'cm=',cm, 'proda=', proda, 'prodm=', prodm)                    
                        if proda*prodm < 0.0:
                            ub = um
                            vb = vm
                            cb = cm
                            auxcsb = auxcsm
                            auxcfb = auxcfm
                            #print('L885 auxcfb=', auxcfb)
                        elif proda*prodm > 0.0:
                            ua = um
                            va = vm
                            ca = cm
                            auxcsa = auxcsm
                            auxcfa = auxcfm
                            proda = prodm
                            #print('L 893 auxcfa=auccfm', auxcfa)
                    u[n] = um
                    v[n] = vm
                    c[n] = cm
                    print("SingleContatBranch function. Coincidence reached: u =", np.round(um, 6),
                     "v =", np.round(vm,6),
                      'c =', np.round(cm,6), '\n')
                       
        if stopcoinc == 0 and c[n] <= cmin: # Faz a correcao no ultimo ponto para que c[n] = cmin
            hnew = (cmin - c[n-1]) / fun.aga(u[n-1], v[n-1],c[n-1])
            u[n] = u[n-1] + hnew*fun.efe(u[n-1], v[n-1], c[n-1])
            v[n] = v[n-1] + hnew*fun.ge(u[n-1], v[n-1], c[n-1])
            c[n] = cmin
            #Retorna o ponto de interseccao com o plano c = cmin
            #print("SingleContatBranch function. cminimum reached: u =", np.round(u[Npontos-1], 6),
                     #"v =", np.round(v[Npontos-1],6),
                      #'c =', np.round(c[Npontos-1],6), '\n')
    # Retorna os valores, tanto para o caso up como para o caso down
    #print('Npontos =', Npontos, 'n=', n,'\n')
    return(u, v, c, Npontos)

'''
# Para testar        
#Plotagem das projecoes do ramo de contato
# Define por coordenadas baricentricas ou nao
mp = ini.baricentrica()
#c0 = 0.5

c0 = 0.0
u0 = 0.4
v0 = 0.4

cmin = 0.0
cmax = 1.0

updow = 1
stopcoinc = 1

u, v, c, Npontos = SingleContatBranch(u0, v0, c0, cmin, cmax, updow, stopcoinc)
#print('u      v  \n')
#for n in range(0,len(u)):
#print(np.round(u[n],4), np.round(v[n],4), '\n')
import matplotlib.pyplot as plt #Para plotagem de graficos   
plt.figure('uc projection')
plt.plot(u0, c0, 'ro')
plt.plot(u,c,'m-')
plt.plot([umin, umax], [c0, c0], 'k--')
#plt.grid()
#plt.show()

# Projecao da curva de contato no plano vc
plt.figure('vc projection')
plt.plot(v0, c0, 'ro')
plt.plot(v,c,'m-')
plt.plot([umin, umax], [c0, c0], 'k--')
#plt.grid()
#plt.show()

# Faz o mapeamento do (u,v) para coordenadas baricentricas se mp = 1
for i in range(0, Npontos):
    u[i], v[i] = fun.map(u[i], v[i], mp)

# Projecao da curva de contato no triangulo de saturacoes c = c0
plt.figure('saturation triangle')
plt.plot(u[0], v[0], 'ro')   
plt.plot(u, v,'m-')
#plt.grid()
#plt.show() 
'''

## Funcao que calcula os pontos do campo de contato comecando num ponto
## U0 sobre uma superficie de coincidencia

def ContatBranchFromCoincidence(u0, v0, c0, cmin, side):
    # U0 = (u0, v0, c0) deve ser um ponto de coincidencia
    #cmin indica o menor nivel c permitido
    #side indica se o ramo qeu inicia numa superficie de coincidencia
    # seguirah para a esquerda (side = -1) ou para a direita (side = 1) de tal superficie
    
    tol = 10**(-3)
    Npontos = 0
    
    lambdac = fun.lbdc(u0, v0, c0)
    lambdas = fun.lbdas(u0, v0, c0)
    lambdaf = fun.lbdaf(u0, v0, c0)
    
    if abs(lambdac - lambdas) > tol and abs(lambdac - lambdaf) > tol:
        print('ContatBranchFromCoincidence function: U0 longe de coincidencias \n')
        print('u0 = ', u0, 'v0 = ', v0, 'c0 = ', c0, 'Npontos = ', Npontos, '\n')
        return(u0, v0, c0, Npontos)
    else:
    
        # Para integrar os contatos (por Euler)
        h = 0.01 # Passo
        Nmax = 10000 # Numero maximo de pontos
        
        Npontos = 1 # Numero de pontos calculados. Ha pelo menos o ponto U_0
        
        u = [] # Para guardar as aproximacoes com c decrescente
        v = [] # Para guardar as aproximacoes 
        c = [] # Para guardar as aproximacoes
    
        auxc = c0 # Define auxc para entrar no loop a seguir
        
        # Testa o sentido a partir do  ponto de coincidencia
        auxu = u0 + h*fun.efe(u0, v0, c0)    
        
        # Fixa o sentido decrescente para u proximo de u0 (esquerda da superfice)   
        if side == -1 and auxu > u0: # u crescente. Inverter o sentido de integracao
            h = -h
                
        if side == 1 and auxu < u0: # u decrescente. Inverter o sentido de integracao
            h = -h
                
        u.append(u0) # Valor inicial
        v.append(v0)
        c.append(c0)
        
        n = 0 # Para controlar o numero de pontos
        while auxc > cmin and n < Nmax:
            # Integra por Euler
            u.append(u[n] + h*fun.efe(u[n], v[n], c[n]))
            v.append(v[n] + h*fun.ge(u[n], v[n], c[n]))
            c.append(c[n] + h*fun.aga(u[n], v[n], c[n]))
           
            auxc = c[n+1] # Para testar o while
            n = n+1
        
        Npontos = np.size(c) # Fixa o n^o de pontos calculados
        if n == Nmax and c[n] > cmin:
            print('ContatBranchFromCoincience function. Nmax = ', Nmax,'insuficiente. \n')
            print('u =', np.round(u[Nmax], 6),
                      "v =", np.round(v[Nmax],6),
                      'c =', np.round(c[Nmax],6), '\n')
        elif c[n] < cmin: # Faz a correcao no ultimo ponto para que c[n] = cmin
            hnew = (cmin - c[n-1]) / fun.aga(u[n-1], v[n-1],c[n-1])
            u[n] = u[n-1] + hnew*fun.efe(u[n-1], v[n-1], c[n-1])
            v[n] = v[n-1] + hnew*fun.ge(u[n-1], v[n-1], c[n-1])
            c[n] = cmin
        
            #Retorna o ponto de interseccao com o plano c = cmin
            print('ContatBranchFromCoincience function. Extremidade: u =', np.round(u[Npontos-1], 6),
                     "v =", np.round(v[Npontos-1],6),
                      'c =', np.round(c[Npontos-1],6), '\n')
                
        #Retorna os valores, tanto para o caso up como para o caso down
        #print('Npontos =', Npontos, 'n=', n,'\n')
        return(u, v, c, Npontos)
'''  
# Para testar        
#Plotagem das projecoes do ramo de contato
# Define por coordenadas baricentricas ou nao
mp = ini.baricentrica()
c0 = 0.5

u0 = 0.645748
v0 = 0.302781

cmin = 0.0


u, v, c, Npontos = ContatBranchFromCoincidence(u0, v0, c0, cmin, -1)
if Npontos >= 1:
    
    import matplotlib.pyplot as plt #Para plotagem de graficos   
    plt.figure('uc projection')
    plt.plot(u0, c0, 'ro')
    plt.plot(u,c,'m-')
    plt.plot([umin, umax], [c0, c0], 'k--')
    #plt.grid()
    #plt.show()
    
    # Projecao da curva de contato no plano vc
    plt.figure('vc projection')
    plt.plot(v0, c0, 'ro')
    plt.plot(v,c,'m-')
    plt.plot([umin, umax], [c0, c0], 'k--')
    #plt.grid()
    #plt.show()
    
    # Faz o mapeamento do (u,v) para coordenadas baricentricas se mp = 1
    for i in range(0, Npontos):
        u[i], v[i] = fun.map(u[i], v[i], mp)
    
    # Projecao da curva de contato no triangulo de saturacoes c = c0
    plt.figure('saturation triangle')
    plt.plot(u[0], v[0], 'ro')   
    plt.plot(u, v,'m-')
    #plt.grid()
    #plt.show()
'''
#2024/01/17 Substituir as chamadas desta funcao por Bissection3
def LinBissection(A, B, c, f, TOL):
    # Find a zero of f(u,v,c) in the linear segment AB with
    # Tolerance TOL  in the plane c
    # A and B are arrays of Ndim components
    # f is a five real variable function
    # Raiz is an Ndim-array that is returned
    
    Ndim = 2
    Raiz = [] #To return the root
    vA = f(A[0],A[1],c)
    
    #vB = f(B[0],B[1],c)
    
    for i in range(0,Ndim):
        Raiz.append(0.5*(A[i] + B[i]))
    vM = f(Raiz[0],Raiz[1],c)
    
    MaxIteracoes = 6
    k = 0
    while abs(vM) > TOL and k <= MaxIteracoes:
        k = k + 1
        if vM*vA < 0.0:
            for i in range(0,Ndim):
                B[i] = Raiz[i]
                Raiz[i] = 0.5*(A[i] + B[i])
            vM = f(Raiz[0],Raiz[1],c)
        else:
            vA = f(Raiz[0],Raiz[1],c)
            for i in range(0,Ndim):
                A[i] = Raiz[i]
                Raiz[i] = 0.5*(A[i] + B[i])
            vM = f(Raiz[0],Raiz[1],c)

    return(Raiz)
    
'''
# Para testar
A = [0.0, 0.0]
B = [1.0, 1.0]
c = 0.0
N = 2
TOL = 0.001
def efe(x,y,c):
    valor = x**2 + y**2 - 1.0
    
    return(valor)

Raiz =  LinBissection(A, B, c, efe, TOL) 
print('Interseccao = ', Raiz)
print('Valor da funcao = ', efe(Raiz[0], Raiz[1],c))
  
'''
#2024/01/17 Substituir as chamadas desta funcao por Bissection5
def LinBissectionsigma(A, B, u0, v0, c0, cnt, f, TOL):  
    #calcula o zero de f(A,B,u0,v0,c0) - cnt
    # A e B devem ser arrays com Ndim componentes
    # f deve ser uma funcao de 5 variaveis reais
    # A raiz deve ser um array com Ndim componentes aramazenada em Raiz
    
    Ndim = 2
    Raiz = [] #To return the root
    vA = f(A[0],A[1],u0,v0,c0) - cnt
    
    #vB = f(B[0],B[1],c)
    
    for i in range(0,Ndim):
        Raiz.append(0.5*(A[i] + B[i]))
    vM = f(Raiz[0],Raiz[1],u0,v0,c0) - cnt
    
    while abs(vM) > TOL:
        if vM*vA < 0.0:
            for i in range(0,Ndim):
                B[i] = Raiz[i]
                Raiz[i] = 0.5*(A[i] + B[i])
            vM = f(Raiz[0],Raiz[1],u0,v0,c0) - cnt
        else:
            vA = f(Raiz[0],Raiz[1],u0,v0,c0) - cnt
            for i in range(0,Ndim):
                A[i] = Raiz[i]
                Raiz[i] = 0.5*(A[i] + B[i])
            vM = f(Raiz[0],Raiz[1],u0,v0,c0) - cnt

    return(Raiz)

'''   
# Para testar

def efe(x,y,u0,v0,c0):
    valor = x**2 + y**2
    
    return(valor)

A = [0.0, 0.0]
B = [1.0, 1.0]
u0 = 0.0
v0 = 0.0
c0 = 0.0
cnt = 1
N = 2
TOL = 0.001

Raiz =  LinBissectionsigma(A, B, u0,v0,c0, cnt, efe, TOL) 
print('Interseccao = ', Raiz)
print('Valor da funcao = ', efe(Raiz[0], Raiz[1],u0,v0,c0))
   
'''
#2024/01/17
def Bissection3(A, B, c, f, TOL, MaxIter):
    # Obtida por uma modificacao da funcao LinBissection(A, B, c, f, TOL):
    # Find a zero of f(u,v,c) in the linear segment AB with
    # Tolerance TOL  in the plane c
    # A and B are arrays of Ndim components
    # f is a three-real variable function
    # MaxIter is the maximum of iteractions
    # Raiz is an Ndim-array that is returned
    
    Ndim = 2
    Raiz = [] #To return the root
    vA = f(A[0],A[1],c)
    
    #vB = f(B[0],B[1],c)
    
    for i in range(0,Ndim):
        Raiz.append(0.5*(A[i] + B[i]))
    vM = f(Raiz[0],Raiz[1],c)
    
    k = 0
    while abs(vM) > TOL and k <= MaxIter:
        k = k + 1
        if vM*vA < 0.0:
            for i in range(0,Ndim):
                B[i] = Raiz[i]
                Raiz[i] = 0.5*(A[i] + B[i])
            vM = f(Raiz[0],Raiz[1],c)
        else:
            vA = f(Raiz[0],Raiz[1],c)
            for i in range(0,Ndim):
                A[i] = Raiz[i]
                Raiz[i] = 0.5*(A[i] + B[i])
            vM = f(Raiz[0],Raiz[1],c)
    #If necessary
    print('Funcoes Auxiliares. Bissection3: Number of Iteractions =', k-1, '\n')
    print('Funcoes Auxiliares. Bissection3: Maximum of Iteractions =', MaxIter, '\n')
    print('Funcoes Auxiliares. Bissection3: Tolerance =', TOL, '\n')

    return(Raiz)
    
'''
# Para testar
A = [0.0, 0.0]
B = [1.0, 1.0]
c = 0.0
N = 2
TOL = 0.001
MaxIter = 6
def efe(x,y,c):
    valor = x**2 + y**2 - 1.0
    
    return(valor)

Raiz =  Bissection3(A, B, c, efe, TOL, MaxIter) 
print('Interseccao = ', Raiz)
print('Valor da funcao = ', efe(Raiz[0], Raiz[1],c))
  
'''
#2024/01/17
def Bissection5(A, B, u0, v0, c0, cnt, f, TOL, MaxIter):  
    # Find a root of f(u,v,u0,v0,c0) - cnt in the linear segment AB
    # A and B are arrays of Ndim components
    # f is a five-real variable function 
    # TOL  is the tolerance
    # MaxIter is the maximum of iteractions
    # Root is an Ndim-array that is returned
    
    Ndim = 2
    Root = [] #To return the root
    vA = f(A[0],A[1],u0,v0,c0) - cnt
    
    #vB = f(B[0],B[1],c)
    
    for i in range(0,Ndim):
        Root.append(0.5*(A[i] + B[i]))
    vM = f(Root[0],Root[1],u0,v0,c0) - cnt
    
    k = 0
    while abs(vM) > TOL and k <= MaxIter:
        k = k + 1
        if vM*vA < 0.0:
            for i in range(0,Ndim):
                B[i] = Root[i]
                Root[i] = 0.5*(A[i] + B[i])
            vM = f(Root[0],Root[1],u0,v0,c0) - cnt
        else:
            vA = f(Root[0],Root[1],u0,v0,c0) - cnt
            for i in range(0,Ndim):
                A[i] = Root[i]
                Root[i] = 0.5*(A[i] + B[i])
            vM = f(Root[0],Root[1],u0,v0,c0) - cnt
    #If necessary
    #print('Funcoes Auxiliares. Bissection5: Number of Iteractions =', k-1, '\n')
    #print('Funcoes Auxiliares. Bissection5: Maximum of Iteractions =', MaxIter, '\n')
    #print('Funcoes Auxiliares. Bissection5: Tolerance =', TOL, '\n')
    
    return(Root)

''' 
# Para testar

def efe(x,y,u0,v0,c0):
    valor = x**2 + y**2
    
    return(valor)

A = [0.0, 0.0]
B = [1.0, 1.0]
u0 = 0.0
v0 = 0.0
c0 = 0.0
cnt = 1.0
N = 2
TOL = 0.001
MaxIter = 6

Raiz =  Bissection5(A, B, u0,v0,c0, cnt, efe, TOL, MaxIter) 
print('Raiz = ', Raiz)
print('f - cnt = ', efe(Raiz[0], Raiz[1],u0,v0,c0) - cnt)    

'''
def ExplCoinc(fam, c): #fam: familia 1=slow ou 2=fast. c: nivel. 
    '''
    Ideia: Usar as formulas explictas das coincidencias (extensoes da fronteira GO)
    do paper Azevedo et al, SIAP-2014
    
    A funcao Coincid acima usa as formulas das diferencas lambdac - lambdas e lambdac - lambdaf
    ''' 
    
    # Define o numero de pontos em cada segmento da curva
    N = 1000
       
    #Viscosidades 
    muw0, muo, mug = ini.viscosidades() 
        
    # Razoes de viscosidades no nivel c
    muwdec = fun.muw(c)
    rw = muwdec / muo
    rg = mug / muo
    ro = 1.0
    rwg = rw + rg
    rtot = rw + rg + ro
    
    # Para guardar os pontos da coincidencia 
    vectoruplus = [] # Sinal positivo na formula (5.3) do paper
    vectorvplus = []    
    
    #Formulas 5.5
    swE1 = rw / mat.sqrt((rw + ro)*rwg)
    swE2 = mat.sqrt(rw/rtot)
    #print('swE1=', swE1, 'swE2=', swE2)
           
    if fam == 1:
        #print('Ts')
        # Passos da parametrizacao na variavel sw
        h1 = swE1 / N #Para sw no intervalo [0, swE1]
       

        #Primeiro Loop cuidando do sinal +        
        for i in range(0, N): # i = N eh considerado a parte para não ter dissc < 0
            sw = i*h1 #Parametrizacao por sg
            
            #Formula 5.4
            Discr = 4.0*(rw*rtot*rwg*(1.0 + rw)*sw**4 - (2*rw*rtot + rg)*rw**2*sw**2 + rw**4)
            #Formula 5.3
            numsg1 = -2*rw*((rw+ro)*sw**2 - rg*sw - rw)
            numsg2 = mat.sqrt(Discr)
            denomsg = 2*rw*((rw + rtot)*sw + 2*rw)
            
            sgplus = (numsg1 + numsg2) / denomsg
            vectoruplus.append(sw)
            vectorvplus.append(1.0 - sgplus - sw)
            
        # Para i = N, como sw = swE1 pode acontecer de Discr < 0. Então fazemos Discr = 0
        #Formula 5.3
        sw = swE1
        Discr = 0.0
        numsg1 = -2*rw*((rw+1.0)*sw**2 - rg*sw - rw)
        numsg2 = mat.sqrt(Discr)
        denomsg = 2*rw*((rw + rtot)*sw + 2*rw)
        
        sgplus = (numsg1 + numsg2) / denomsg
        vectoruplus.append(sw)
        soplus = 1.0 - sgplus - sw
        vectorvplus.append(soplus)
        
        # Segundo Loop considerando o sinal -
        for i in range(1, N): # i = N eh considerado a parte pra não ter dissc < 0
            sw = swE1 - i*h1 #Parametrizacao por sw
            
            #Formula 5.4
            Discr = 4.0*(rw*rtot*rwg*(1.0 + rw)*sw**4 - (2*rw*rtot + rg)*rw**2*sw**2 + rw**4)

            #Formula 5.3
            numsg1 = -2*rw*((rw+ro)*sw**2 - rg*sw - rw)
            numsg2 = mat.sqrt(Discr)
            denomsg = 2*rw*((rw + rtot)*sw + 2*rw)
            
            sgplus = (numsg1 - numsg2) / denomsg
            vectoruplus.append(sw)
            vectorvplus.append(1.0 - sgplus - sw)
         
                                                
    if fam == 2:
        #Primeiro faz a parte do ramo que cruza o lado so = 0 comecando com so = 0
        #Segundo faz o ramo que cruza o lado sg = 0 terminando com sg = 0
        #print('Tf')
        # Para determinar o ponto final na face so = 0
        coefa = rw*rtot - rw**2 - 2*rw - rw*rg - rg
        coefb = 2*rw**2 + 2*rg*rw + 2*rw - 2*rw*rtot
        coefc = rw*rtot - rg*rw - rw**2
        Delta = coefb**2 - 4*coefa*coefc
        RaizDelta = mat.sqrt(Delta)

        
        #Limite superior para sw no ramo que intercepta o lado so = 0
        swminus = (-coefb - RaizDelta) / (2*coefa) 
        sgminus = 1.0 - swminus
        #print('swminus=', swminus)
        h2 = (swminus - swE2) / N #Para sw no intervalo [swE2, swminus]
        
        
        # Para i de 0 ateh N-1 use a formula de Discr evitando Discr < 0
        for i in range(0, N):
            sw = swminus - i*h2 #Parametrizacao por sg
            
            #Formula 5.4
            Discr = 4.0*(rw*rtot*rwg*(1.0 + rw)*sw**4 - (2*rw*rtot + rg)*rw**2*sw**2 + rw**4)
            #Formula 5.3
            numsg1 = -2*rw*((rw+1.0)*sw**2 - rg*sw - rw)
            numsg2 = mat.sqrt(Discr)
            denomsg = 2*rw*((rw + rtot)*sw + 2*rw)

            sgplus = (numsg1 + numsg2) / denomsg
            soplus = 1.0 - sgplus - sw

            vectoruplus.append(sw)
            vectorvplus.append(soplus)
        
        # Para sw = swE2 , ou i = N, facca a parte para evitar que Discr < 0
        sw = swE2 #Parametrizacao por sw
        
        #Formula 5.4
        Discr = 0.0
        #Formula 5.3
        numsg1 = -2*rw*((rw+1.0)*sw**2 - rg*sw - rw)
        numsg2 = mat.sqrt(Discr)
        denomsg = 2*rw*((rw + rtot)*sw + 2*rw) 
        sgplus = (numsg1 + numsg2) / denomsg
        soplus = 1.0 - sgplus - sw        #print ('sgplus=', swglus, 'soplus=', soplus, '\n')

        vectoruplus.append(sw)
        vectorvplus.append(soplus)
        ###########################
        
        # Limite superior para sw no ramo que intercepta o lado sg = 0            
        swsg0 = mat.sqrt(rw/(rw+1)) 
        h2 = (swsg0 - swE2)/N
               
        # Para j de 1 a N use a formula com Discr > 0. j = 0 jah foi feito
        for j in range(1,N+1):
            sw = swE2 + j*h2 #Parametrizacao por sw
            
            #Formula 5.4
            Discr = 4.0*(rw*rtot*rwg*(1.0 + rw)*sw**4 - (2*rw*rtot + rg)*rw**2*sw**2 + rw**4)
            #Formula 5.3
            numsg1 = -2*rw*((rw+1.0)*sw**2 - rg*sw - rw)
            numsg2 = mat.sqrt(Discr)
            denomsg = 2*rw*((rw + rtot)*sw + 2*rw)
            
            sgminus = (numsg1 - numsg2) / denomsg
            sominus = 1.0 - sgminus - sw
            vectoruplus.append(sw)
            vectorvplus.append(sominus)
            
    # if save == 1:
    #     print('See File Coincidence')
    #     Npoints = np.size(vectoruplus)
    #     f = open('Coincidence','w')
    #     for i in range(Npoints):
    #         f.write(str(vectoruplus[i])+'  '+str(vectorvplus[i])+'\n')
    #     f.close()
    
    return(vectoruplus, vectorvplus)

'''
#Para testar
import matplotlib.pyplot as plt #Para plotagem de graficos    
plt.figure('saturation triangle', figsize=(8,6), dpi=80)

fam = 1
c = 0.5
u, v = ExplCoinc(fam, c)

N = len(u)
print('option = ', fam, 'len(u) =', N)
for i in range(N):
    u[i], v[i] = fun.map(u[i], v[i], mp)

if fam == 1:
    plt.plot(u, v, 'b-')
if fam == 2:
    plt.plot(u, v, 'r')

'''
  
def Cextremos(c):
    # c deve ser uma lista com as coordenadas c de um segmento de contato, por exemplo
    
    #Determinacao dos primeiros indices apos um valor extremo (minimo ou maximo)
    # em um vetor de pontos c.
    #Retorna a quantidade de mudancas no crescimento e os indices correspondentes
   
    NP = len(c)
    if NP > 1:
        Coinc = []
        aux1 = 1.0 # Para testar crescimento de c e pontos de coincidencia
        if c[1] - c[0] < 0:
            aux1 = -aux1 
     
        for i in range(1,NP):
            aux2 = (c[i] - c[i-1]) * aux1
            if aux2 < 0: # Houve inversao no crescimento de c
                aux1 = -aux1
                Coinc.append(i)
    Qindex = len(Coinc) # Quantidade de pontos de coicnidencia

    return(Qindex, Coinc) #Retorna a quantidade e os indices dos pontos de coincidencai 

'''
#Para testar
v = [-1,2,3, 2, 4, -1]
N, indices = Cextremos(v)
print('\n No de extremos', N, '\n indices=', indices)
'''      
                
def Coincidencias(a,b,c):
    # Testa as diferencas lambdac - lambdas e lambdac - lambdaf
    # para saber qual a coincidencia
    dif1 = abs(Difcs(a, b, c))
    dif2 = abs(Difcf(a, b, c))
    
    if dif1 < dif2:
        coinc = 1
    else:
        coinc = 2
            
    return(coinc)

'''
#Para testar

a, b, c = 0.26398, 0.70329, 0.0
Coincid = Coincidencias(a,b,c)
print('\n coinc = ', Coincid)
''' 

def triangletest(u, v, umin, umax, vmin, vmax, domain):
# Define se (u,v) pertence ao dominio deinteresse interior ao triangulo de saturacoes
    #umin, umax, vmin, vmax, domain = ini.dominio()
    
    if umin <= u <= umax and vmin <= v <= vmax and 0 <= u + v <= 1.0:
        return(1)
    else:
        return(0)

'''
#Para testar

u0 = 0.5
v0 = 0.2
dominio = triangletest(u0,v0)
print('\n Teste 1 sim, ', dominio)

u0 = 0.5
v0 = 0.5
dominio = triangletest(u0,v0)
print('\n Teste 2 sim, ', dominio)

u0 = 0.5
v0 = 0.7
dominio = triangletest(u0,v0)
print('\n Teste 3 nao ', dominio)
'''

#2024-01-17
def find_closest_point(matrix, x0, y0):
    # Convert the matrix to a numpy array for easier manipulation
    matrix = np.array(matrix)
    
    # Calculate the distances between each point in the matrix and (x0, y0)
    distances = (matrix[0] - x0)**2 + (matrix[1] - y0)**2
    
    # Find the index of the point with the minimum distance
    min_distance_index = np.argmin(distances)
    
    # Get the coordinates of the closest point
    closest_point = (matrix[0][min_distance_index], matrix[1][min_distance_index])
    
    # Remove the closest point from the matrix
    matrix = np.delete(matrix, min_distance_index, axis=1)
    
    return matrix, closest_point, min_distance_index

'''
# To test
line0 = [1,-2,3,4,5]
line1 = [7,8,-9,10,11]

A = []

A.append(line0)
A.append(line1)

print('Original Matrix = ', A, '\n')

B, Min = find_closest_point(A, -1.0, 1.0)
print('Final Matrix = ', B, 'Minimum Point =', Min)
'''

# Funcao que obtem os pontos da elipse -curva de nivel lamndac = lambdac(U0) explicitamente
def Ellipse(u0, v0, c0, N):
    # Devolve os vetores (eliu, eliv)

    # N numero de subdivisoes do intervalo [umin umax], sendo umin e umax onde du/dv = 0
    
    
    # Ver a possibiidade de não definir un e usar algo do tipo "extend ou concatenate" para acrescentar
    # os pontos correspondentes a raiz negativa
   
    # Arrays para as saidas
    up = []
    vp = []
    un = []
    vn = []
    
    muw0, muo, mug = ini.viscosidades()
    muw0 = fun.muw(c0)
    
    ro = muw0 / muo
    rg = muw0 / mug
    K = fun.lbdc(u0, v0, c0)
    
    a = K*(1.0 + rg)
    b = K*(ro + rg)
    c = 2*K*rg
    d = -(1.0 + 2*K*rg)
    e = -2*K*rg
    f = K*rg
    
# Calcula o intervalo u que define a elipse. Pontos onde du/dv = 0   
      
    aux1 = 4*b*d - 2*c*e
    aux2 = (4*a*b - c**2)
    aux3 = 4*b*f - e**2
    deltinha = aux1**2 - 4 * aux2 * aux3
    aux4 = np.sqrt(deltinha)
    
    umin = (-aux1 - aux4 ) / (2*aux2)
    umax = (-aux1 + aux4 ) / (2*aux2)
    #print('\n umin = ', umin)
    #print('\n umax= ', umax)
    
    vmin = -(c*umin + e) / (2*b)
    vmax = -(c*umax + e) / (2*b)
    #print('\n vmin = ', vmin)
    #print('\n vmax = ', vmax)
        
    A = b
    
    # 
    u  = np.linspace(umin, umax, N)
    
    up.append(umin)
    vp.append(vmin)  #Primeiro ponto com raiz positiva derivada em v nula
      
    for n in range(1, N-1): # Evita a divisao por zero onde du/dv = 0
        up.append(u[n])
        un.append(u[n])
        #up = np.append(u[n])
        #un = np.append(u[n])
        
        B = c*up[n] + e
        C = a*up[n]**2 + d*up[n] + f
        Delta = B**2 - 4*A*C
        vp.append((-B + np.sqrt(Delta)) / (2*A))
        vn.append((-B - np.sqrt(Delta)) / (2*A))
        #vp = np.append((-B + np.sqrt(Delta)) / (2*A))
        #vn = np.append((-B - np.sqrt(Delta)) / (2*A))
    
    up.append(umax)
    vp.append(vmax)    #Ultimo ponto com raiz positiva derivada em v nula
    #up = np.append(umax)
    #vp = np.append(vmax)
    
    # Inverte pontos com raiz negativa e junta ao da raiz positiva
    for n in range(1, N-1):
        up.append(un[N-2-n])
        vp.append(vn[N-2-n])
        #up = np.append(un[N-2-n])
        #vp = np.append(vn[N-2-n])
        
        
    # To draw the complete elipse
    up.append(umin)
    vp.append(vmin)
    
    # Exchange the closest point to the initial (u0, v0) by (u0, v0)
    #matrix = []
    #matrix.append(up)
    #matrix.append(vp)
    NewElipse, ClosestPoint, index = find_closest_point([up, vp], u0, v0)

    up[index] = u0
    vp[index] = v0
    
    # To start the picture from (u0, v0)
    upnew = []
    vpnew = []

    for k in range(0, len(up) - index):
        upnew.append(up[index + k])
        vpnew.append(vp[index + k])
    
    for m in range(0, index):
        upnew.append(up[m])
        vpnew.append(vp[m])
    

    return(upnew, vpnew)
'''
# Only to test
#import matplotlib as mpl
#from mpl_toolkits.mplot3d import Axes3D
#import numpy as np
import matplotlib.pyplot as plt

u0 = 0.4246216873843006#0.46
v0 = 0.319988910113042#0.5
c0 = 0.0#0.5928522166980628
N = 200
UP, VP = Ellipse(u0, v0, c0, N)

# Para fechar a Elipse na plotagem
UP.append(UP[0])
VP.append(VP[0])

#print('\n\n\n UP', UP, '\n\n\n VP', VP)
plt.figure('saturation triangle')
plt.plot(UP, VP, 'm-')
plt.plot([u0, 0.4255072525169559],[v0, 0.31891115030911066],'k-')

plt.grid()
'''


def LambdacEqualSigma(u0, v0, c0, u, v):
    # Testa as diferencas lambdac - sigma ao longo do da curva definida por (u, v)
    # para c = constante
    # para determinar pontos com sigma*((u0,v0); (u,v)) = lambdac(u0, v0, c0)
    
    # Em processo de melhoria. Trocar teste de velocidadds por teste de pontos
    # || Hug(u0, v0, c0, u,v) - (u, v) || < 10**(-5)
    #def H(u,v,u0,v0,c0):
        #return (fw(u,v,c0) - fw(u0,v0,c0))*(v-v0) - (fo(u,v,c0) - fo(u0,v0,c0))*(u-u0)
    ucoinc = []
    vcoinc = []
    N = len(u)
    i = 0
    
    print('FunctionLambdacEqualSigma: lbdc(u0,v0,c0)=', fun.lbdc(u0,v0,c0), 'sigma(u0,v0, u[0], v[0], c0)=', sigma(u0,v0,u[0],v[0],c0), '\n')
    dif = Difcsigma(u0, v0, u[0], v[0], c0) # verifica a diferença lambda^c- sigma no primeiro ponto
    print('Function LambdacEqualSigma: i =', i,'\n')
    print('Function LambdacEqualSigma: u[0]= ', u[0], '\n')
    print('Function LambdacEqualSigma: u0= ', u0, '\n')

    if abs(dif) <= 10**(-5):
        ucoinc.append(u[0])  # Valor inicial
        vcoinc.append(v[0])
        print('Function LambdacEqualSigma: dif1=', dif,'\n')
        print('Function LambdacEqualSigma: ucoinc1=', ucoinc,'\n')
    
    while abs(dif) <= 10**(-5) and i < N-1:
        i = i + 1
        dif = Difcsigma(u0, v0, u[i], v[i], c0)
        print('Function LambdacEqualSigma: i =', i,'\n')
        print('FunctionLambdacEqualSigma: lbdc(u0,v0,c0)=', fun.lbdc(u0,v0,c0), 'sigma(u0,v0, u[i], v[i], c0)=', sigma(u0,v0,u[i],v[i],c0), '\n')

        print('Function LambdacEqualSigma: dif2=', dif,'\n')
          
    i0 = i + 1
    
    if dif > 0:
        difaux = 1
    else:
        difaux = -1
    
    for i in range(i0, N):
        dif = Difcsigma(u0, v0, u[i], v[i], c0) # verifica a diferença lambda^c(u,v,c0)- sigma((u,v,c0);(u0,v0,c0))
        print('Function LambdacEqualSigma: dif3=', dif,'\n')
        teste = difaux*dif
        if teste <= 0:
            if i0 == 1:
                ucoinc = np.array(u[i])
                vcoinc = np.array(v[i])
            else:
                ucoinc = np.append(ucoinc, u[i])
                vcoinc = np.append(vcoinc, v[i])
            print('Function LambdacEqualSigma: ucoinc2=', ucoinc,'\n')
            while abs(dif) <= 10^(-5) and i < N:
                i = i + 1
                dif = Difcsigma(u0, v0, u[i], v[i], c0)
                print('Function LambdacEqualSigma: dif4=', dif,'\n')
                
            #precaucao para o caso de dif ter sido zero    
            if dif > 0:
                difaux = 1
            else:
                difaux = -1
                
            ucoinc = np.append(ucoinc,0.5)
            vcoinc = np.append(vcoinc,0.2)
            print('Function LambdacEqualSigma: ucoinc3=', ucoinc,'\n')
                
    return(ucoinc, vcoinc)

'''
#Para testar

u0 = 0.4246216873843006#0.46
v0 = 0.319988910113042#0.5
c0 = 0.0#0.5928522166980628
# u0 = 0.455
# v0 = 0.4951
# c0 = 0.4
N = 10
UP, VP = Elipse(u0, v0, c0, N)

U, V = LambdacEqualSigma(u0, v0, c0, UP, VP)
print('\n Teste Function LambdacEqualSigma: Ucoinc = ', U)
print('\n Teste unction LambdacEqualSigma: Vcoinc = ', V)

h = 0.01
N = 10000
cmin = 0.0
cmax = 1.0
updow = 1
stopcoinc = 1


u1, v1, c1, Npontos1 = SingleContatBranch(U[0], V[0], c0, cmin, cmax, updow, stopcoinc)
#u2, v2, c2, Npontos2 = SingleContatBranch(U[1], V[1], c0, cmin, cmax, updow, stopcoinc)


import matplotlib.pyplot as plt #Para plotagem de graficos 

mp = 1  
plt.figure('uc projection')
plt.plot(u0, c0, 'ro')
plt.plot(u1,c1,'m-')
#plt.plot(u2,c2,'b-')
plt.grid()


# Projecao da curva de contato no plano vc
plt.figure('vc projection')
plt.plot(v0, c0, 'ro')
plt.plot(v1, c1,'m-')
#plt.plot(v2, c2,'b-')
plt.grid()
#plt.show()

# Faz o mapeamento do (u,v) para coordenadas baricentricas se mp = 1
for i in range(0, Npontos1):
    u1[i], v1[i] = fun.map(u1[i], v1[i], mp)

# Projecao da curva de contato no triangulo de saturacoes c = c0
plt.figure('saturation triangle')
plt.plot(u1[0], v1[0], 'ro')   
plt.plot(u1, v1,'m-')

# for j in range(0, Npontos2):
#     u2[j], v2[j] = fun.map(u2[j], v1[j], mp)

# plt.plot(u2[0], v2[0], 'ro')   
# plt.plot(u2, v2,'b-')
#plt.grid()
#plt.show()
'''

##################

# Consertar e completar esta funcao abaixo e usar mesmas ideias em outras. 06/01/24.# 
# Retornar apenas pontos que estão no dominio
def LambdacEqualLambdacZero(u0, v0, c0, u, v):
    # Testa as diferencas lambdac(U) - lambda(U0) ao longo da curva (u,v)
    # contendo o ponto U0=(u0, v0) no nivel c = c0 e
    # retorna os pontos (ucoinc, vcoinc) em que lambdac = lambdac(U0)
    # Em processo de melhoria 06/01/14

    ucoinc = []
    vcoinc = []
    ucoinc.append(u0) # O ponto U0 deve satisfazer a igualdade
    vcoinc.append(v0)
    K = 1
    print('LambdacEqualLambdacZero1: K = len(ucoinc) = ', K, '\n')
    print('LambdacEqualLambdacZero1: ucoinc = ', ucoinc, '\n')
    
    N = len(u)
    i = 0
    
    difaux = 1
    
    if ( (u0 - u[0])**2 + (v0 - v[0])**2 ) >= 10*(-1):
        dif = Diflbdc(u0, v0, c0, u[0], v[0])
        print('LambdacEqualLambdacZero2: i = ', i, 'dif = ', dif,'\n')
        difaux = dif*difaux
        if difaux > 0:
            difaux = 1
        else:
            difaux = -1
        
        # if abs(dif) <= 10**(-5):
        #     ucoinc = np.append(ucoinc, u[0])
        #     vcoinc = np.append(vcoinc, v[0])
        #     print('i = ', i,'vcoinc = ', vcoinc)
        #     K = K + 1
        #     print('LambdacEqualLambdacZero3: K = len(ucoinc) = ', K, '\n')
        #     i = 1
        
        while i < N: # Procura troca de sinal
            dif = Diflbdc(u0, v0, c0, u[i], v[i])
            print('LambdacEqualLambdacZero4: i = ', i, 'dif = ', dif,'\n')
            
            if np.abs(dif) <= 10**(-5):
                teste = 0 # Testar se eh um dos pontos anteriores
                for k in range(K):
                    dist = (ucoinc[k] - u[i])**2 + (vcoinc[k] - v[i])**2
                    if dist < 10**(-1):
                        teste = 1
                if teste == 0:
                    ucoinc = np.append(ucoinc, u[i])
                    vcoinc = np.append(vcoinc, v[i])
                    K = K + 1
                    print('LambdacEqualLambdacZero5: K = len(ucoinc) = ', K, '\n')
                    print('LambdacEqualLambdacZero5: ucoinc = ', ucoinc, '\n')
            i = i + 1
        
    return(ucoinc, vcoinc)


'''
#Para testar

u0 = 0.4246216873843006#0.46
v0 = 0.319988910113042#0.5
c0 = 0.0#0.5928522166980628
# u0 = 0.455
# v0 = 0.4951
# c0 = 0.4
N = 20
UP, VP = Elipse(u0, v0, c0, N)

U, V = LambdacEqualSigma(u0, v0, c0, UP, VP)
print('\n Ucoinc = ', U)
print('\n Vcoinc = ', V)

h = 0.01
N = 10000
cmin = 0.0
cmax = 1.0
updow = 1
stopcoinc = 1


u1, v1, c1, Npontos1 = SingleContatBranch(U, V, c0, cmin, cmax, updow, stopcoinc)
#u2, v2, c2, Npontos2 = SingleContatBranch(U[1], V[1], c0, cmin, cmax, updow, stopcoinc)


import matplotlib.pyplot as plt #Para plotagem de graficos 

mp = 1  
plt.figure('uc projection')
plt.plot(u0, c0, 'ro')
plt.plot(u1,c1,'m-')
#plt.plot(u2,c2,'b-')
plt.grid()


# Projecao da curva de contato no plano vc
plt.figure('vc projection')
plt.plot(v0, c0, 'ro')
plt.plot(v1, c1,'m-')
#plt.plot(v2, c2,'b-')
plt.grid()
#plt.show()

# Faz o mapeamento do (u,v) para coordenadas baricentricas se mp = 1
for i in range(0, Npontos1):
    u1[i], v1[i] = fun.map(u1[i], v1[i], mp)

# Projecao da curva de contato no triangulo de saturacoes c = c0
plt.figure('saturation triangle')
plt.plot(u1[0], v1[0], 'ro')   
plt.plot(u1, v1,'m-')

# for j in range(0, Npontos2):
#     u2[j], v2[j] = fun.map(u2[j], v1[j], mp)

# plt.plot(u2[0], v2[0], 'ro')   
# plt.plot(u2, v2,'b-')
#plt.grid()
#plt.show()
'''

def HugoniotContactIntersection(u0, v0, c0, u, v):
    # (u0, v0): Ponto base da Hugoniot no plano c = c0
    # Testa se H(u0,v0,c0) intersecta o segmento definido pelos arrays (u, v)
    # que definem alguma curva

    # Devolve os arrays (uint, vint) de pontos de interseccao
    
    # Em processo de melhoria
    uint = []
    vint = []
    
    
    Nelup = len(u)
    
    if Nelup > 0: # Se houver pontos da parte da raiz positiva na elipse
        dif = fun.H(u[0],v[0],u0,v0,c0) # Testa sinal
        if abs(dif) < 10**(-5): # Estah na Hugoniot
            uint.append(u[0])
            vint.append(v[0])
        
       
        i = 1 # Para testar outros possiveis pontos
        while abs(dif) < 10**(-5) and i < Nelup:
            dif = fun.H(u[i],v[i],u0,v0,c0) # Testa se (ulup[i], elvp[i]) estah
            # na Hugoniuot de (u0, v0, c0)
            i = i + 1
        # Identificou o primeiro ponto com dif \ne 0    
        
        # Define um sinal caso dif \ne 0
        if dif < 0:
            difaux = -1
        if dif > 0:
            difaux = 1
        
        while i < Nelup:
            #print('\n\n i =', i)
            dif = fun.H(u[i],v[i],u0,v0,c0) # Testa se (ulup[i], elvp[i]) estah
            # na Hugoniuot de (u0, v0, c0)
            difaux = dif*difaux
                   
            if difaux < 0:
                uint.append(u[i])
                vint.append(v[i])
            difaux = dif
            #rint('\n LambdacEqualSigma: i = ', i, 'dif=', dif, '\n')
            i = i + 1
                
    return(uint, vint)

'''
#Para testar

u0 = 0.4
v0 = 0.5
c0 = 0.0
N = 500
up, vp, un, vn = LambdacContourExpl(u0, v0, c0, N)
#print('\n UElipse= ', u)
#print('\n VElipse= ', v)

Ucoincid, Vcoincid = HugoniotContactIntersection(u0, v0, c0, up, vp)
#print('\n\n Ucoincid = ', Ucoincid, '\n\n Vcoincid = ', Vcoincid, '\n\n')
a = fun.lbdc(Ucoincid[0], Vcoincid[0],c0)

Ucoincidn, Vcoincidn = HugoniotContactIntersection(u0, v0, c0, un, vn)

b = fun.lbdc(Ucoincidn[0], Vcoincidn[0],c0)

c = fun.lbdc(u0, v0,c0)

print('\n lambdac1 =', a, '\n lambdac2 = ', b, '\n lambdc(U0) =', c, '\n')
'''
def cutoffTsTf(u,v,c0,coinc):
    # 2023-11-23
    # Cutoff the segment along the curve (u,v) at the level c = c0 before coincidences
    # u and v are arrays with the same length
    # if coinc = 1, lambda^c = lambda^s
    # if coinc = 2, lambda^c = lambda^f
    # return (unew, vnew) 

    tol = 10**(-3)
    
    N = len(u)
    
    unew = [] # to save the "new" segment
    vnew = []
    
    q = 0 # to define unew, vnew
    
    if coinc == 1:
        teste = Difcs(u[q],v[q], c0) #lambdac - lambdas
        
        if teste > 0:
            sinal = 1
            prod = 1.0
        elif teste < 0:
            sinal = -1
            prod = 1.0
        else:
            sinal = 0
            prod = 0.0
            unew.append(u[0])
            vnew.append(v[0])
        
        while prod > 0.0 and q < N-1:
            unew.append(u[q])
            vnew.append(v[q])      
            q = q + 1

            a, b = u[q],v[q]
            teste = Difcs(a, b, c0)
            prod = teste*sinal

        # Bissection
        if q > 0:
            [unew[q-1],vnew[q-1]] = LinBissection([unew[q-1],vnew[q-1]], [a, b], c0, Difcs, tol)

    else:
        teste = Difcf(u[q],v[q], c0) #lambdac - lambdaf
         
        if teste > 0:
            sinal = 1
            prod = 1.0
        elif teste < 0:
            sinal = -1
            prod = 1.0
        else:
            sinal = 0
            prod = 0.0
            unew.append(u[0])
            vnew.append(v[0])
            
            
        while prod > 0 and q < N-1:
            unew.append(u[q])
            vnew.append(v[q])      
            q = q + 1

            a, b = u[q],v[q]
            teste = Difcs(a, b, c0)
            prod = teste*sinal
            
        # Bissection
        if q > 0:
            [unew[q-1],vnew[q-1]] = LinBissection([unew[q-1],vnew[q-1]], [a, b], c0, Difcs, tol)
    
    return(unew, vnew)

def ball(a, b, c, r):
    # (a, b, c) center
    # r radius
    N = 15
    M = 2*N
    theta = np.linspace(-np.pi, np.pi, N)
    phi = np.linspace(-np.pi/2.0, np.pi/2.0, M)

    theta, phi = np.meshgrid(theta, phi)
    u = a + r*np.cos(theta)*np.cos(phi)
    v = b + r*np.sin(theta)*np.cos(phi)
    w = c + r*np.sin(phi)

    return (u, v, w)

'''
#Para testar

import matplotlib.pyplot as plt

a = 0.3
b = 0.3
c = 0.3
r = 0.0085

ub, vb, cb = ball(a, b, c, r)

#ax = plt.axes(projection ='3d')
ax = ini.ambiente3d()
ax.plot_surface(ub, vb, cb) 
plt.show()
'''
#2024-01-17
def HugContactOrigins(u0, v0, c0, u, v):
    # Obtain the intersection points of the Hugoniot curve from (u0,v0,c0) with
    # the curve parametrized by (u, v) on the plane c = c0
    
    ucoinc = []
    vcoinc = []
    N = len(u)
    i = 0
    
    #print('Function HugContactBranches: lbdc(u0,v0,c0)=', fun.lbdc(u0,v0,c0), '\n')
    value = fun.H(u[0], v[0], u0, v0, c0)
   
    while abs(value) <= 10**(-12) and i < N-1:
        ucoinc.append(u[i])
        vcoinc.append(v[i])
        i = i + 1
        value = fun.H(u[i], v[i], u0, v0, c0)
        
    i0 = i
    
    if value > 0:
        difaux = 1
    else:
        difaux = -1
    
    for i in range(i0, N):
        value = fun.H(u[i], v[i], u0, v0, c0)
        teste = difaux*value
        dist = (u[i] - u0)**2 + (v[i] - v0)**2
        
        if teste < 0:
            difaux = - difaux
            A = [u[i-1], v[i-1]]
            B = [u[i], v[i]]
            TOL = 10**(-12)
            MaxIter = 10

            U = Bissection5(A, B, u0, v0, c0, 0.0, fun.H, TOL, MaxIter)
            if dist >=10**(-12): # To avoid repetitions
                ucoinc.append(U[0])
                vcoinc.append(U[1])
    
    ucoinc = np.array(ucoinc) 
    vcoinc = np.array(vcoinc)
                     
    return(ucoinc, vcoinc)
'''
# Only to test

import matplotlib.pyplot as plt #Para plotagem de graficos 

#mp = 0
h = 0.01
N = 10000
cmin = 0.0
cmax = 1.0
updow = 1
stopcoinc = 1


u0 = 0.39
v0 = 0.56
c0 = 0.0
N = 1000
UP, VP = Ellipse(u0, v0, c0, N)


U, V = HugContactOrigins(u0, v0, c0, UP, VP)
print('Ucoinc=', U, '\n')
print('Vcoinc=', V, '\n')

plt.figure('saturation triangle')
plt.plot(UP, VP, 'g-') # Plot the elipse labdc = labdc(U0) on the plane c = c0

#uc, vc, cc, Npontosc = SingleContatBranch(u0, v0, c0, cmin, cmax, updow, stopcoinc)

plt.plot(u0, v0, 'ko') # Plot the contact branch for U0, from c = c0 to a concidence surface


for m in range(len(U)):

    plt.plot(U[m], V[m], 'ro') # Plot the origin points for contact branches


for n in range(len(U)):
    test = triangletest(U[n], V[n]) # verify if (U[n], V[n]) is in the sat triangle

    if test == 1:
        u1, v1, c1, Npontos1 = SingleContatBranch(U[n], V[n], c0, cmin, cmax, updow, stopcoinc)
        
        # Projection of the contact branch on plane uc
        plt.figure('uc projection')
        plt.plot(U[n], c0, 'ro')
        plt.plot(u1,c1,'m-')
        #plt.plot(uc, cc, 'k-')
        #plt.plot(u0, c0, 'ko')
        #plt.plot(u2,c2,'b-')
        plt.grid()
        
        
        # Projection of the contact branch on plane vc
        plt.figure('vc projection')
        plt.plot(V[n], c0, 'ro')
        plt.plot(v1, c1,'m-')
        #plt.plot(vc, cc,'k-')
        #plt.plot(v0, c0, 'ko')
        plt.grid()
        #plt.show()
        
        # Map (u,v) to baricentric coordinates, if  mp = 1
        for i in range(0, Npontos1):
            u1[i], v1[i] = fun.map(u1[i], v1[i], mp)
        
        # Projection of the contact branch on the plane c = c0
        plt.figure('saturation triangle')
        plt.plot(U[n], V[n], 'ro')   
        plt.plot(u1, v1,'m-')
'''


# ###2024-02-13 Eliminar repeticoes com funcoes acima
# def ContactIntegration(u0, v0, c0, cmin, cmax, h, tol, Nmax):
#     #cmin deve nao ser superior a cmax
#     #Se stopcoinc = 1, entao forca a parada nas superficies de coincidencia    
    
    
#     #print('SingleContatBranch function. u0=', u0, 'v0=', v0, 'c0=',c0)
    
#     # Adjust  c0, if c0 is near cmin or near cmax
#     aga0 = fun.aga(u0, v0, c0)
    
#     if abs( (c0 - cmin) / aga0 ) <= abs(h):
#         c0 = cmin
        
#     if abs( (c0 - cmax) / aga0 ) <= abs(h):
#         c0 = cmax

#     n = 0 # To count the number of storaged points in a while
    
#     # To store the approximations
#     u = []
#     v = []
#     c = []
    
#     # Store the initial point       
#     u = np.array([u0])
#     v = np.array([v0])
#     c = np.array([c0])
    
#     un = u0
#     vn = v0
#     cn = c0

        
#     difcs = Difcs(u0,v0,c0) #fun.lbdc(u0, v0, c0) - fun.lbdas(u0, v0, c0)
#     difcf = Difcf(u0,v0,c0) #fun.lbdc(u0, v0, c0) - fun.lbdaf(u0, v0, c0)
#     #print('SingleContatBranch function. difcs(U0), difcf(U0) =', difcs0, difcf0, '\n')
    
#     # To prevent coincidence point
#     u1 = u0
#     v1 = v0
#     c1 = c0
    
#     # To define signs
#     if np.abs(difcs) <= tol or np.abs(difcf) <= tol:
#         n = 1

#         while np.abs(difcs) <= tol:
#             u1 = u1 + h*fun.efe(u1, v1, c1)
#             v1 = v1 + h*fun.ge(u1, v1, c1)
#             c1 = c1 + h*fun.aga(u1, v1, c1)
#             difcs = Difcs(u1,v1,c1)
            
           
#         while np.abs(difcf) <= tol:
#             u1 = u1 + h*fun.efe(u1, v1, c1)
#             v1 = v1 + h*fun.ge(u1, v1, c1)
#             c1 = c1 + h*fun.aga(u1, v1, c1)
#             difcf = Difcf(u1,v1,c1)
            
#         u = np.append(u, u1)
#         v = np.append(v, v1)
#         c = np.append(c, c1)
          
#     # The difference signs are defined
#     if difcs > 0.0:
#         signcs = 1
#     elif difcs < 0.0:
#         signcs = -1
          
#     if  difcf > 0.0:
#         signcf = 1
#     elif difcf < 0.0:
#         signcf = -1
             
#     un = u1
#     vn = v1
#     cn = c1
     
#     auxc = cn            
#     while cmin <= auxc <= cmax and n < Nmax:
#         # Integration by the Euler Method
#         un = un + h*fun.efe(un, vn, cn)
#         vn = vn + h*fun.ge(un, vn, cn)
#         cn = cn + h*fun.aga(un, vn, cn)
       
#         u = np.append(u, un)
#         v = np.append(v, vn)
#         c = np.append(c, cn)

#         auxc = cn # Para testar o while        
#         n = n + 1
        
#         #Test if it is near a singular point
#         dist = abs(un) + abs(vn) + abs(cn)
#         if dist < tol:
#             singular = 1
#             print('Near a singular point',  'dist = ', dist, '\n')
#             break
#         else:
#            singular = 0

#         difcs = Difcs(un,vn,cn) #lbdb - fun.lbdas(u[n+1],v[n+1],c[n+1])
#         #print('Integration: dif_cs =', difcs, '\n')
#         signcs = difcs * signcs
#         if signcs <= 0.0:
#             #print('Coincidence State lambda^c = lambda^s.')
#             break               
#         signcs = difcs
    
#         difcf = Difcf(un,vn,cn) #lbdb - fun.lbdaf(u[n+1],v[n+1],c[n+1]) 
#         #print('Integration: dif_cf =', difcf, '\n')
#         signcf = difcf * signcf

#         if signcf <= 0.0:
#             #print('Coincidence State lambda^c = lambda^f.')
#             break
#         signcf = difcf
          
#     if cmin <= auxc <= cmax and n == Nmax:
#         if Nmax > 2:
#             print('ContactIntegration function. Nmax = ', Nmax,'not enough. \n')
    
#     elif cn < cmin: # Correction if c[n] = cmin
#         hnew = (cmin - c[n-1]) / fun.aga(u[n-1], v[n-1], c[n-1])
#         u[n] = u[n-1] + hnew*fun.efe(u[n-1], v[n-1], c[n-1])
#         v[n] = v[n-1] + hnew*fun.ge(u[n-1], v[n-1], c[n-1])
#         c[n] = cmin
        
#     elif cn > cmax: # Correction if c[n] > cmax
#         hnew = (cmax - c[n-1]) / fun.aga(u[n-1], v[n-1],c[n-1])
#         u[n] = u[n-1] + hnew*fun.efe(u[n-1], v[n-1], c[n-1])
#         v[n] = v[n-1] + hnew*fun.ge(u[n-1], v[n-1], c[n-1])
#         c[n] = cmax

       
#     return(u, v, c, singular)

## 2024-02-23. Copiando de Teste 6
#Campo de contato
def F(t,y):
    eq1 = fun.efe(y[0], y[1], y[2])
    eq2 = fun.ge(y[0], y[1], y[2])
    eq3 = fun.aga(y[0], y[1], y[2])
    
    return [eq1, eq2, eq3]

## 2024-02-23. Substituindo a partir de Teste6 usando solve_ivp
def ContactIntegration(u0, v0, c0, cmin, cmax, h, tol, Nmax):
    #cmin deve nao ser superior a cmax
    #Se stopcoinc = 1, entao forca a parada nas superficies de coincidencia    
    
    
    #print('SingleContatBranch function. u0=', u0, 'v0=', v0, 'c0=',c0)
    
    # Adjust  c0, if c0 is near cmin or near cmax
    aga0 = fun.aga(u0, v0, c0)
    
    if abs( (c0 - cmin) / aga0 ) <= abs(h):
        c0 = cmin
        
    if abs( (c0 - cmax) / aga0 ) <= abs(h):
        c0 = cmax

    n = 0 # To count the number of storaged points in a while
    
    # To store the approximations
    u = []
    v = []
    c = []
    
    # Store the initial point       
    u = np.array([u0])
    v = np.array([v0])
    c = np.array([c0])
    
    un = u0
    vn = v0
    cn = c0

        
    difcs = Difcs(u0,v0,c0) #fun.lbdc(u0, v0, c0) - fun.lbdas(u0, v0, c0)
    difcf = Difcf(u0,v0,c0) #fun.lbdc(u0, v0, c0) - fun.lbdaf(u0, v0, c0)
    #print('SingleContatBranch function. difcs(U0), difcf(U0) =', difcs0, difcf0, '\n')
    
    # To prevent coincidence point
    u1 = u0
    v1 = v0
    c1 = c0
    
    t = np.arange(0, 2*h, h) # Malha em t
    
    # To define signs
    if np.abs(difcs) <= tol or np.abs(difcf) <= tol:
        n = 1

        while np.abs(difcs) <= tol:
            #u1 = u1 + h*fun.efe(u1, v1, c1)
            #v1 = v1 + h*fun.ge(u1, v1, c1)
            #c1 = c1 + h*fun.aga(u1, v1, c1)
            sol = solve_ivp(F, [0, h], [u1, v1, c1], t_eval = t)
            u1 = sol.y[0]
            v1 = sol.y[1]
            c1 = sol.y[2]
            difcs = Difcs(u1,v1,c1)
            
           
        while np.abs(difcf) <= tol:
            #u1 = u1 + h*fun.efe(u1, v1, c1)
            #v1 = v1 + h*fun.ge(u1, v1, c1)
            #c1 = c1 + h*fun.aga(u1, v1, c1)
            sol = solve_ivp(F, [0, h], [u1, v1, c1], t_eval = t)
            u1 = sol.y[0]
            v1 = sol.y[1]
            c1 = sol.y[2]
            
            difcf = Difcf(u1,v1,c1)
            
        u = np.append(u, u1)
        v = np.append(v, v1)
        c = np.append(c, c1)
          
    # The difference signs are defined
    if difcs > 0.0:
        signcs = 1
    elif difcs < 0.0:
        signcs = -1
          
    if  difcf > 0.0:
        signcf = 1
    elif difcf < 0.0:
        signcf = -1
             
    un = u1
    vn = v1
    cn = c1
     
    auxc = cn          
    while cmin <= auxc <= cmax and n < Nmax:
        # Integration by the Euler Method
        #un = un + h*fun.efe(un, vn, cn)
        #vn = vn + h*fun.ge(un, vn, cn)
        #cn = cn + h*fun.aga(un, vn, cn)
        sol = solve_ivp(F, [0, 2*h], [un, vn, cn], t_eval = t)
        un = sol.y[0][1]
        vn = sol.y[1][1]
        cn = sol.y[2][1]        
       
        u = np.append(u, un)
        v = np.append(v, vn)
        c = np.append(c, cn)

        auxc = cn # Para testar o while        
        n = n + 1
        
        #Test if it is near a singular point
        dist = abs(un) + abs(vn) + abs(cn)
        if dist < tol:
            singular = 1
            print('Near a singular point',  'dist = ', dist, '\n')
            break
        else:
           singular = 0

        difcs = Difcs(un,vn,cn) #lbdb - fun.lbdas(u[n+1],v[n+1],c[n+1])
        #print('Integration: dif_cs =', difcs, '\n')
        signcs = difcs * signcs
        if signcs <= 0.0:
            #print('Coincidence State lambda^c = lambda^s.')
            break               
        signcs = difcs
    
        difcf = Difcf(un,vn,cn) #lbdb - fun.lbdaf(u[n+1],v[n+1],c[n+1]) 
        #print('Integration: dif_cf =', difcf, '\n')
        signcf = difcf * signcf

        if signcf <= 0.0:
            #print('Coincidence State lambda^c = lambda^f.')
            break
        signcf = difcf
          
    if cmin <= auxc <= cmax and n == Nmax:
        if Nmax > 2:
            print('ContactIntegration function. Nmax = ', Nmax,'not enough. \n')
    
    elif cn < cmin: # Correction if c[n] = cmin
        hnew = (cmin - c[n-1]) / fun.aga(u[n-1], v[n-1], c[n-1])
        u[n] = u[n-1] + hnew*fun.efe(u[n-1], v[n-1], c[n-1])
        v[n] = v[n-1] + hnew*fun.ge(u[n-1], v[n-1], c[n-1])
        c[n] = cmin
        
    elif cn > cmax: # Correction if c[n] > cmax
        hnew = (cmax - c[n-1]) / fun.aga(u[n-1], v[n-1],c[n-1])
        u[n] = u[n-1] + hnew*fun.efe(u[n-1], v[n-1], c[n-1])
        v[n] = v[n-1] + hnew*fun.ge(u[n-1], v[n-1], c[n-1])
        c[n] = cmax

       
    return(u, v, c, singular)
'''
# Para testar        
#Plotagem das projecoes do ramo de contato
# Define por coordenadas baricentricas ou nao
mp = 0# 1 # to baricentricas
#c0 = 0.5

c0 = 0.0 #0.4#0.47633542532979983
u0 = 0.694# 0.4#0.19842524359890665
v0 = 0.278#0.22013451204279894

cmin = 0.0
cmax = 1.0

umin = 0.0
umax = 1.0
vmin = 0.0
vmax = 1.0


h = 0.1
tol = 10**(-5)
Nmax = 10000


u, v, c, sing = ContactIntegration(u0, v0, c0, cmin, cmax, h, tol, Nmax)
#print('u      v  \n')
#for n in range(0,len(u)):
#print(np.round(u[n],4), np.round(v[n],4), '\n')
import matplotlib.pyplot as plt #Para plotagem de graficos   
plt.figure('uc projection')
plt.plot(u0, c0, 'ro')
plt.plot(u,c,'m-')
plt.plot([umin, umax], [c0, c0], 'k--')
#plt.grid()
#plt.show()

# Projecao da curva de contato no plano vc
plt.figure('vc projection')
plt.plot(v0, c0, 'ro')
plt.plot(v,c,'m-')
plt.plot([umin, umax], [c0, c0], 'k--')
#plt.grid()
#plt.show()

# Faz o mapeamento do (u,v) para coordenadas baricentricas se mp = 1
Npontos = len(u)
print('teste singleContactBranch :NPontos = ', Npontos, "\n")
for i in range(0, Npontos):
    u[i], v[i] = fun.map(u[i], v[i], mp)

# Projecao da curva de contato no triangulo de saturacoes c = c0
plt.figure('saturation triangle')
plt.plot(u0, v0, 'ro')   
plt.plot(u, v,'m-')
#plt.grid()
#plt.show() 

'''
#2024-02-13 Cido.
def ContactBranches(laxclass, u0, v0, c0, cmin, cmax, h, tol, Nmax):
    #laxclass is the Lax classification of the contact
       
    #Matrix = [] # To store the branch
    
    # Local branch
    difcs = Difcs(u0, v0, c0)
    difcf = Difcf(u0, v0, c0)
    
    if difcs < 0:
        testclass = 1
    elif difcs > 0.0 and difcf < 0.0:
        testclass = 2
    elif difcf > 0:
        testclass = 3
    elif difcs == 0.0:
        testclass = 12
        print('U0 on a coincidence lambda^c = lambda^s')
    elif difcf == 0.0:
        testclass = 23
        print('U0 on a coincidence lambda^c = lambda^f')
    if abs(testclass - laxclass) >  0:    
        print('Not a', laxclass, ', but a', testclass,'contact branch from U0 \n')
        u1 = [u0]
        v1 = [v0]
        c1 = [c0]
                     
    if laxclass == 1 and difcs <= 0.0:
        if h > 0.0:
            h = - h
        # To test if the shock classification is correct
        Naux = 2
        uaux, vaux, caux, singular = ContactIntegration(u0, v0, c0, cmin, cmax, h, tol, Naux)
        testeaux = Difcs(uaux[1], vaux[1], caux[1])

        if testeaux > 0:
            h = - h # Probably U0 lies on T^s with c0 minimum.
            
        if c0 >= cmin:
            u1, v1, c1, singular = ContactIntegration(u0, v0, c0, cmin, cmax, h, tol, Nmax)
            
        if cmin < c0 < cmax and difcs < 0.0: # there is a second contact segment to be calculated
            h = - h #revert the integration sense
            u2, v2, c2, singular = ContactIntegration(u0, v0, c0, cmin, cmax, h, tol, Nmax)
            
            ## order inversion of U1
            unew = u1[::-1]
            vnew = v1[::-1]
            cnew = c1[::-1]
            
            u1 = np.concatenate((unew, u2[1:]))
            v1 = np.concatenate((vnew, v2[1:]))
            c1 = np.concatenate((cnew, c2[1:]))

    if laxclass == 3 and difcf >= 0.0:
        if h > 0.0:
            h = - h
            
        if c0 >= cmin:
            u1, v1, c1, singular = ContactIntegration(u0, v0, c0, cmin, cmax, h, tol, Nmax)
            
        if cmin < c0 < cmax and difcf > 0.0: # there is a second contact segment to be calculated
            h = - h #revert the integration sense
            u2, v2, c2, singular = ContactIntegration(u0, v0, c0, cmin, cmax, h, tol, Nmax)
            
            ## order inversion of U1
            unew = u1[::-1]
            vnew = v1[::-1]
            cnew = c1[::-1]
            
            u1 = np.concatenate((unew, u2[1:]))
            v1 = np.concatenate((vnew, v2[1:]))
            c1 = np.concatenate((cnew, c2[1:]))

    if laxclass == 2 and difcs >= 0.0 and difcf <= 0.0:
        if h < 0.0:
            h = - h
            
        # To test if the shock classification is correct
        Naux = 2
        uaux, vaux, caux, singular = ContactIntegration(u0, v0, c0, cmin, cmax, h, tol, Naux)
        if caux[1] < cmin:
            print('U0 with c0 = cmin') # TO DO: Check it
        else:
            u1, v1, c1, singular = ContactIntegration(u0, v0, c0, cmin, cmax, h, tol, Nmax)
                        
        if cmin < c0 < cmax and difcs > 0.0 and difcf < 0.0: # there is a second contact segment to be calculated
            h = - h #revert the integration sense
            u2, v2, c2, singular = ContactIntegration(u0, v0, c0, cmin, cmax, h, tol, Nmax)
            
            ## order inversion of U1
            unew = u1[::-1]
            vnew = v1[::-1]
            cnew = c1[::-1]
            
            u1 = np.concatenate((unew, u2[1:]))
            v1 = np.concatenate((vnew, v2[1:]))
            c1 = np.concatenate((cnew, c2[1:]))

        #### TO DO: implement the case with a nonlocal contact branch

        #Matrix = np.append([laxclass, u1, v1, c1])
            
    #return(np.array(Matrix))
    return(u1, v1, c1)     


'''
 Matrix [[0,[1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1]], [1,[2,2,2,2,2,2],[3,3,3,3,3,3],[4,4,4,4,4,4]]
    Matrix[1][1] = [2,2,2,2,2,2]
    Matrix[1][2] = [3,3,3,3,3,3]
    Matrix[1][3] = [4,4,4,4,4,4]
    Matrix[0][0] = 0
    Matrix[1][0] = 1
'''

'''

# Para teste

u0 = 0.1#0.47633542532979983
v0 = 0.19#0.2#0.19842524359890665
c0 = 0.4#0.22013451204279894

cmin = 0.0
cmax = 1.0

umin = 0.0
umax = 1.0
vmin = 0.0
vmax = 1.0

h = 0.1
tol = 10**(-9)
Nmax = 10000

laxclass = 1


u3, v3, c3 = ContactBranches(laxclass, u0, v0, c0, cmin, cmax, h, tol, Nmax)

import matplotlib.pyplot as plt #Para plotagem de graficos  

if laxclass == 1:
    color = 'b-'
elif laxclass == 2:
    color = 'm-'
else:
    color = 'r'
    
plt.figure('uc projection')
plt.plot(u0, c0, 'ko')
plt.plot(u3, c3, color)
plt.plot([umin, umax], [c0, c0], 'k--')

# Projecao da curva de contato no triangulo de saturacoes c = c0
plt.figure('saturation triangle')
plt.plot(u0, v0, 'ko')   
plt.plot(u3, v3, color)
#plt.grid()
#plt.show() 


plt.figure('vc projection')
plt.plot(v0, c0, 'ko')
plt.plot(v3, c3, color)
plt.plot([vmin, vmax], [c0, c0], 'k--')


print('ufirst = ', u3[0], 'vfirst = ', v3[0], 'cfirst = ', c3[0])
print('ulast = ', u3[len(u3)-1], 'vlast = ', v3[len(u3)-1], 'clast = ', c3[len(u3)-1])
'''