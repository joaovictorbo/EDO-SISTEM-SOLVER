# sistema12.py
import sistema12_funcoes as funcoes
import sistema12_plot as plot

def run(u0, v0, c0):
    """
    Executa o sistema 12 com as condições iniciais fornecidas.
    
    Parâmetros:
        u0, v0, c0: Condições iniciais para u, v, c.
    """
    # Resolver as trajetórias
    sol, sol2 = funcoes.resolver_trajetoria(u0, v0, c0)
    
    # Dividir as trajetórias
    trajetorias1 = funcoes.dividir_trajetorias(sol)
    trajetorias2 = funcoes.dividir_trajetorias(sol2)
    
    # Plotar as trajetórias
    plot.plotar_trajetorias(trajetorias1, trajetorias2, u0, v0, c0)
