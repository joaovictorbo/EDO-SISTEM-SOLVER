# main.py
import sistemaadmissibilidade as adm
import sistemaadmissibilidadeplot as plot

def run_admissibilidade(u0, v0, z0, u_R, v_R, z_R, f_R, g_R):
    print("Executando sistema de admissibilidade com condições iniciais:")
    print(f"u0={u0}, v0={v0}, z0={z0}")
    
    # Resolver a trajetória com as condições iniciais
    t_span = (0, 10)
    t_values, y_values = adm.resolver_trajetoria(u0, v0, z0, u_R, v_R, z_R, f_R, g_R, t_span)
    
    # Filtrar trajetórias (opcional, dependendo da aplicação)
    trajetorias = [y_values.T]  # Adapte para incluir apenas partes desejadas
    
    # Plotar as trajetórias
    plot.plotar_trajetorias(trajetorias)

if __name__ == "__main__":
    # Condições iniciais
    u0, v0, z0 = 0.5, 0.2, 0.1
    u_R, v_R, z_R = 0.7, 0.5, 0.4
    f_R, g_R = 1.0, 1.0  # Exemplos de valores de referência

    # Executar o sistema
    run_admissibilidade(u0, v0, z0, u_R, v_R, z_R, f_R, g_R)