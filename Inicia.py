# -*- coding: utf-8 -*-
"""
dominio(): fixa o dominio retangular no plano c = c0
concentrations(): fixa as concentracoes cmin e cmax
viscosidades(): fixa as viscosidades muw0(inicial da agua), muo e mug
baricentrica(): fixa se vai usar coordenadas baricentricas ou nao
DadosIntegracao(): fixa o comprimento do passo de integracao 
                   e os numero de pontos com passo axes = figure.add_subplot(projection='3d')h
ambiente3d(): Seta o ambiente 3d e devolve "axes" para plotagem 3d

"""
import matplotlib.pyplot as plt

def dominio():
# Define uma regiao R retangular para calculos.
# Quando tal regiao R ultrapassar os limites do triangulo, as
# funcoes especificas devem fazer esta restricao adicional de considerar
# pontos apenas na interseccao da regiao R com o triangulo
    
    
    triang = 1 # 1 caso seja calculado em todo o triangulo.
    
    if triang == 1: # Calculos no triangulo todo
        umin = 0.0
        umax = 1.0
        vmin = 0.0
        vmax = 1.0
    else: # Calculo em regiao retangular especifica para "zoom"
        umin = 0.0#0.315#0.2
        umax = 1.0#0.415#0.6
        vmin = 0.0#0.0#0.515#0.35
        vmax = 1.0#1.0#0.630#0.75
        
    return(umin, umax, vmin, vmax, triang)
    

def concentrations(): #Viscosidade da agua
    cmin = 0.0#0.5928522166980628
    #cmax = 1.0#0.666616 #0.284124  # #1.0 Normal# L1 0.28017 #L2 e L3 0.284124
    cmax = 1.0#0.66#0.28017 #L1
    return(cmin,cmax)

def viscosidades(): #Viscosidades da agua, do oleo e do gas
    muw0 = 1.0 #Viscosidade inicial sem polimero
    muo = 4.0#9.5
    mug = 0.25#0.45
    return(muw0, muo, mug)
    
def baricentrica(): #Mapeamento ou nao?
    x = 0# x = 1 se quiser coordenadas baricentricas
    return(x)
    
def DadosIntegracao(): #Para integracao dos contatos
    N = 200000 # Numero de pontos na curva integral do contato para h > 0
    h = 0.02#  # Passo de integração
    return(N, h)

def bib():
    lib = 1 # 0 LapTop,  1 DeskTop
    return(lib)

def ambiente3d():
    figure = plt.figure('PrismDomain')
    lib = bib()
    print('\n \n Inicia.bib(): lib = 1 if matplotlib3.2-Laptop')
    print('               lib = 0 if not matplotlib3.2-desktop.')
    print('               Actual value of lib =', lib, '\n\n')
    if lib == 0:
        axes = figure.gca(projection='3d')
    else:
        axes = figure.add_subplot(projection='3d')  
    return(axes)
    
####################
# Definir funcao que armazena estilos d elinhas, dashed, densed dashed , solid,etc
# from collections import OrderedDict

# linestyles_dict = OrderedDict(
#     [('solid',               (0, ())),
#      ('loosely dotted',      (0, (1, 10))),
#      ('dotted',              (0, (1, 5))),
#      ('densely dotted',      (0, (1, 1))),

#      ('loosely dashed',      (0, (5, 10))),
#      ('dashed',              (0, (5, 5))),
#      ('densely dashed',      (0, (2, 1))),

#      ('loosely dashdotted',  (0, (3, 10, 1, 10))),
#      ('dashdotted',          (0, (3, 5, 1, 5))),
#      ('densely dashdotted',  (0, (3, 1, 1, 1))),

#      ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
#      ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
#      ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])

############################
