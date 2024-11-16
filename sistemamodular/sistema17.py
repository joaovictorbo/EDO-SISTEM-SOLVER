import sistema17_funcoes as funcoes
import sistema17_plot as plot

def run(u0, v0, c0):
    # Resolve as trajetórias
    sol, sol2 = funcoes.resolver_trajetoria(u0, v0, c0)

    # Divide as trajetórias em segmentos
    trajetorias1 = funcoes.dividir_trajetorias(sol)
    trajetorias2 = funcoes.dividir_trajetorias(sol2)

    # Plota as trajetórias
    plot.plotar_trajetorias(trajetorias1, trajetorias2, u0, v0, c0)
