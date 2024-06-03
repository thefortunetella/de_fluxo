include("cfb.jl")
using Main.CFB

# Função para medir e imprimir o tempo de execução
function measure_time(f, args...)
    start = time()
    result = f(args...)
    println("Tempo de execução: ", time() - start, " segundos")
    return result
end

# Escolha do arquivo a ser utilizado
sistema = "ieee300"
arquivo = "ieee300.txt"

# Carregar o grafo a partir do arquivo
g = read_edgelist(arquivo)

# Definir os parâmetros para o DE
k_values = [1, 2, 3]
pop_size = 100
crossover_rate = 0.7
beta_min = 0.2
beta_max = 0.8
iter_num = 10

# Executar o DE para diferentes valores de k e medir o tempo de execução
for k in k_values
    println("Executando DE para k = $k")
    measure_time(cfb_de, g, k, pop_size, crossover_rate, beta_min, beta_max, iter_num, sistema)
end

# script para chamar as funções com os novos parâmetros definir os valores de tensão e impedância e passar esses parâmetros para a função cfb_de

include("cfb.jl")
using Main.CFB

# Função para medir e imprimir o tempo de execução
function measure_time(f, args...)
    start = time()
    result = f(args...)
    println("Tempo de execução: ", time() - start, " segundos")
    return result
end

# Escolha do arquivo a ser utilizado
sistema = "ieee300"
arquivo = "ieee300.txt"

# Carregar o grafo a partir do arquivo
g = read_edgelist(arquivo)

# Definir os parâmetros para o DE
k_values = [1, 2, 3]
pop_size = 100
crossover_rate = 0.7
beta_min = 0.2
beta_max = 0.8
iter_num = 10

# Definir tensões e impedâncias (exemplos, você pode ajustá-los conforme necessário)
voltages = Complex{Float64}[1.0 + 0.0im for i in 1:nv(g)]
impedances = Complex{Float64}[0.01 + 0.05im for i in 1:ne(g)]

# Executar o DE para diferentes valores de k e medir o tempo de execução
for k in k_values
    println("Executando DE para k = $k")
    measure_time(cfb_de, g, k, pop_size, crossover_rate, beta_min, beta_max, iter_num, voltages, impedances, sistema)
end

