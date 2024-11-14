# main.py
import sistema12
import newton

def main():
    # Apresentar opções ao usuário
    print("Escolha o sistema para executar:")
    print("1 - Executar Sistema 12")
    print("2 - Executar Newton")
    print("3 - Executar ambos")
    escolha = input("Digite o número da sua escolha: ")

    # Obter as condições iniciais do usuário
    u0 = float(input("Digite o valor inicial para u: "))
    v0 = float(input("Digite o valor inicial para v: "))
    c0 = float(input("Digite o valor inicial para c: "))

    # Executar de acordo com a escolha do usuário
    if escolha == "1":
        print("Executando Sistema 12 com as condições iniciais fornecidas...")
        sistema12.run(u0, v0, c0)
    elif escolha == "2":
        print("Executando Newton com as condições iniciais fornecidas...")
        newton.run_newton(u0, v0, c0)
    elif escolha == "3":
        print("Executando Sistema 12 e Newton com as condições iniciais fornecidas...")
        sistema12.run(u0, v0, c0)
        newton.run_newton(u0, v0, c0)
    else:
        print("Escolha inválida. Por favor, tente novamente.")

if __name__ == "__main__":
    main()
