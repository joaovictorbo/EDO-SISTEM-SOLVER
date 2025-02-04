# newton.py
import newton_funcoes as funcoes
import newton_plot as plot
import sistema12_funcoes as sistema12

def run_newton(u0, v0, c0):
    print("Executando Newton com condições iniciais:", u0, v0, c0)
    
    # Calcular a trajetória com as condições iniciais
    sol, sol2 = funcoes.resolver_trajetoria(u0, v0, c0)
    sol3, sol4 = sistema12.resolver_trajetoria(u0, v0, c0)
    # Dividir a trajetória para obter apenas as partes dentro do triângulo
    trajetorias_dentro = funcoes.dividir_trajetorias(sol)
    trajetorias_dentro2 = funcoes.dividir_trajetorias(sol2)
    # trajetorias_dentro3 = sistema12.dividir_trajetorias(sol3)
    # trajetorias_dentro4 = sistema12.dividir_trajetorias(sol4)
    
    # Plotar apenas as trajetórias dentro do triângulo
    plot.plotar_trajetorias(trajetorias_dentro + trajetorias_dentro2)
