
module CFB

using Graphs: Graph, SimpleGraph, nv, ne, adjacency_matrix
using Graphs: incidence_matrix, rem_edge!, add_edge!, add_vertices!
using Graphs: is_connected, edges, src, dst
using LinearAlgebra: transpose, diagm, pinv
using IterTools: subsets
using StatsBase: sample, sortperm
using CSV, DataFrames, Tables
using Random: seed!

export read_edgelist
export cfb_de
export current_flow_betweenness

seed!(2021)

function read_edgelist(s::String)
    f = open(s, "r")
    g = SimpleGraph(Int64)
    nodes = Set([])
    edges = Set([])

    while !eof(f)
        line = readline(f)
        numbers = [parse(Int64, m.match) for m in eachmatch(r"\d+", line)]
        if length(numbers) >= 2
            src, dst = numbers[1], numbers[2]
            push!(nodes, src)
            push!(nodes, dst)
            push!(edges, [src, dst])
        end
    end
    close(f)

    add_vertices!(g, length(nodes))
    for e in edges
        add_edge!(g, e[1], e[2])
    end

    return g
end

function current_flow_betweenness(g::Graph)
    n = nv(g)
    m = ne(g)
    betweenness = zeros(Float64, n)
    A = adjacency_matrix(g)
    b = -transpose(incidence_matrix(g, oriented=true))
    n_b = (n - 1) * (n - 2)
    L = diagm(0 => reduce(+, A, dims=[1])[1:n]) - A
    L_tilde = L[2:n, 2:n]
    C = vcat(zeros(1, n), hcat(zeros(n - 1, 1), pinv(L_tilde)))
    F = b * C
    for j in 1:m
        pos = sortperm(sortperm(-F[j, :]))
        e = findall(b[j, :] .!= 0)
        s = e[1]
        d = e[2]
        for i in 1:n
            betweenness[s] += (i - pos[i]) * F[j, i]
            betweenness[d] += (n + 1 - i - pos[i]) * F[j, i]
        end
    end
    norm_vec = 1:n
    betweenness -= norm_vec
    betweenness .+= 1
    betweenness .*= 2 / n_b
    return betweenness
end

function edge_indices_k_tuples(m::Integer, k::Integer)
    subs = subsets(1:m, k)
    edge_index_matrix = zeros(Integer, length(subs), k)
    for (i, s) in enumerate(subs)
        edge_index_matrix[i, :] = s
    end
    return edge_index_matrix
end

function edges_k_tuples(g::Graph, indices::Matrix{Integer})
    edgs = [e for e in edges(g)]
    n_subsets, size_subsets = size(indices)
    edge_matrix = zeros(Integer, n_subsets, 2 * size_subsets)
    for i in 1:n_subsets
        for (j, e) in enumerate(indices[i, :])
            edg = edgs[e]
            edge_matrix[i, 2*j-1] = src(edg)
            edge_matrix[i, 2*j] = dst(edg)
        end
    end
    return edge_matrix
end

function edit_graph(g::Graph, edgs::Vector{Integer})
    g_bkp = Graph(g)
    n_edges = length(edgs)
    n_edges = Int64(n_edges / 2)
    for i in 1:n_edges
        rem_edge!(g_bkp, edgs[2*i-1], edgs[2*i])
    end
    return g_bkp
end

function check_valid_removals(g::Graph, tuples::Matrix{Integer})
    n_subsets, size_subsets = size(tuples)
    valids = [true for i in 1:n_subsets]
    for i in 1:n_subsets
        g_edit = edit_graph(g, tuples[i, :])
        valids[i] = is_connected(g_edit)
    end
    return valids
end

function filter_valid_removals(tuples::Matrix{Integer}, valids::Vector{Bool})
    n_subsets, size_subsets = size(tuples)
    valid_indices = [i for i in 1:n_subsets if valids[i]]
    num_valids = length(valid_indices)
    edge_tuples_matrix = zeros(Integer, num_valids, size_subsets)
    for (i, ti) in enumerate(valid_indices)
        edge_tuples_matrix[i, :] = tuples[ti, :]
    end
    return edge_tuples_matrix
end

function filter_invalid_removals(tuples::Matrix{Integer}, valids::Vector{Bool})
    n_subsets, size_subsets = size(tuples)
    invalid_indices = [i for i in 1:n_subsets if !valids[i]]
    num_invalids = length(invalid_indices)
    edge_tuples_matrix = zeros(Integer, num_invalids, size_subsets)
    for (i, ti) in enumerate(invalid_indices)
        edge_tuples_matrix[i, :] = tuples[ti, :]
    end
    return edge_tuples_matrix
end

function betweenness_from_removals(g::Graph, tuples::Matrix{Integer})
    n = nv(g)
    n_subsets, _ = size(tuples)
    betweenness = zeros(Float64, n_subsets, n)
    for i in 1:n_subsets
        g_edit = edit_graph(g, tuples[i, :])
        betweenness[i, :] = current_flow_betweenness(g_edit)
    end
    return betweenness
end

function betweenness_deltas(betweenness::Matrix{Float64}, reference::Vector{Float64})
    n_subsets, n = size(betweenness)
    deltas = zeros(Float64, n_subsets, n)
    for i in 1:n_subsets
        deltas[i, :] = abs.(betweenness[i, :] .- reference)
    end
    return deltas
end

function de_cost_function(deltas::Matrix{Float64})
    n_subsets, _ = size(deltas)
    costs = vec(sum(deltas, dims=2))
    return costs
end

function de_mutation(m::Integer, k::Integer, pop_size::Integer, indices_pop::Matrix{Integer}, beta_min=0.2::Float64, beta_max=0.8::Float64)
    coefs = rand(-m:m, pop_size, k)
    beta = rand(pop_size, k) .* (beta_max - beta_min) .+ beta_min
    indices_pop_mut = round.(Integer, indices_pop .+ (beta .* coefs))
    indices_pop_mut = clamp.(indices_pop_mut, 1, m)
    return indices_pop_mut
end

function de_crossover(pop_size::Integer, pop_indices::Matrix{Integer}, pop_indices_mut::Matrix{Int64}, crossover_rate=0.5::Float64)
    flags = rand(pop_size) .< crossover_rate
    for c in 1:pop_size
        if flags[c]
            pop_indices[c, :] = pop_indices_mut[c, :]
        end
    end
    return pop_indices
end

function de_iter!(g::Graph, pop_indices::Matrix{Integer}, ref_bets::Vector{Float64}, crossover_rate::Float64, beta_min::Float64, beta_max::Float64, ultima=false::Bool)
    m = ne(g)
    edgs = edges_k_tuples(g, edge_indices_k_tuples(m, 1))
    e = edges_k_tuples(g, pop_indices)
    n_pop, k = size(pop_indices)
    v = check_valid_removals(g, e)
    while sum(v) == 0
        pop_indices = sample_initial_pop(m, n_pop, k)
        e = edges_k_tuples(g, pop_indices)
        v = check_valid_removals(g, e)
    end
    ve = filter_valid_removals(e, v)
    ive = filter_invalid_removals(e, v)
    bets = betweenness_from_removals(g, ve)
    deltas = betweenness_deltas(bets, ref_bets)
    costs = de_cost_function(deltas)
    p = sortperm(-costs)
    costs = costs[p]
    pop_indices = pop_indices[p, :]
    if !ultima
        best_edg = pop_indices[1, :]
        pop_indices_mut = de_mutation(m, k, size

(pop_indices, 1), pop_indices, beta_min, beta_max)
        pop_indices_cross = de_crossover(size(pop_indices, 1), pop_indices, pop_indices_mut, crossover_rate)
        pop_indices = pop_indices_cross
        pop_indices[1, :] = best_edg
        return nothing, nothing, nothing, nothing
    else
        global_norms = sum(abs.(deltas), dims=1)
        deltas_by_tuple = sum(abs.(deltas), dims=2)
        deltas_by_edge = zeros(Float64, m)
        d = Dict([0, 0] => 0)
        for i in 1:m
            d[edges(g)[i]] = i
        end
        for i in 1:size(ve, 1)
            for j in 1:k
                edg = ve[i, 2*j-1:2*j]
                idx = d[edg]
                deltas_by_edge[idx] += deltas_by_tuple[i]
            end
        end
        return ve, ive, DataFrame(V=1:length(global_norms), DELTA=global_norms), DataFrame(DELTA=deltas_by_edge), vec(deltas_by_tuple)
    end
end

function cfb_de(g::Graph, k::Integer, pop_size::Integer, crossover_rate::Float64, beta_min::Float64, beta_max::Float64, iter_num::Integer = 10, arquivo_saida::String = "result")
    ref_cfb = current_flow_betweenness(g)
    m = ne(g)
    edge_index_pop = sample_initial_pop(m, pop_size, k)
    for i in 1:iter_num-1
        de_iter!(g, edge_index_pop, ref_cfb, crossover_rate, beta_min, beta_max)
    end
    ve, ive, global_norms_df, edge_norms_df, local_norms = de_iter!(g, edge_index_pop, ref_cfb, crossover_rate, beta_min, beta_max, true)
    global_norms_vector = Vector{Float64}(global_norms_df[:, "DELTA"])
    edge_norms_vector = Vector{Float64}(edge_norms_df[:, "DELTA"])
    export_de_results(g, ve, global_norms_vector, edge_norms_vector, local_norms, ive, k, arquivo_saida, "de")
end

function export_de_results(g::Graph, valid_tuples::Matrix{Integer}, globals::Vector{Float64}, edge_globals::Vector{Float64}, locals::Vector{Float64}, disconnects::Matrix{Integer}, k::Integer, filename::String, metodo::String)
    dir = string(metodo, "_", filename, "_", string(k))
    dir_bkp = pwd()
    if !isdir(dir)
        mkdir(dir)
    end
    cd(dir)
    CSV.write("valid_tuples.csv", Tables.table(valid_tuples), writeheader=false)
    CSV.write("disconnects.csv", Tables.table(disconnects), writeheader=false)
    CSV.write("local_deltas.csv", Tables.table(locals), writeheader=false)
    CSV.write("vertex_global_deltas.csv", Tables.table(globals), writeheader=false)
    CSV.write("edge_global_deltas.csv", Tables.table(edge_globals), writeheader=false)
    cd(dir_bkp)
end

end

# FUNÇAO FLUXO DE POTENCIA:

function calc_power_flow(g::Graph, voltages::Vector{Complex{Float64}}, impedances::Vector{Complex{Float64}})
    m = ne(g)
    power_flows = zeros(Complex{Float64}, m)
    edges_list = collect(edges(g))

    for (i, edge) in enumerate(edges_list)
        src_node = src(edge)
        dst_node = dst(edge)
        voltage_diff = voltages[src_node] - voltages[dst_node]
        power_flows[i] = voltage_diff / impedances[i]
    end

    return power_flows
end

#Integrando o cálculo de fluxo de potência no algoritmo de evolução diferencial(modificaçao necessaria na função de_iter!

function de_iter!(g::Graph, pop_indices::Matrix{Integer}, ref_bets::Vector{Float64}, crossover_rate::Float64, beta_min::Float64, beta_max::Float64, impedances::Vector{Complex{Float64}}, voltages::Vector{Complex{Float64}}, ultima=false::Bool)
    m = ne(g)
    edgs = edges_k_tuples(g, edge_indices_k_tuples(m, 1))
    e = edges_k_tuples(g, pop_indices)
    n_pop, k = size(pop_indices)
    v = check_valid_removals(g, e)
    while sum(v) == 0
        pop_indices = sample_initial_pop(m, n_pop, k)
        e = edges_k_tuples(g, pop_indices)
        v = check_valid_removals(g, e)
    end
    ve = filter_valid_removals(e, v)
    ive = filter_invalid_removals(e, v)
    bets = betweenness_from_removals(g, ve)
    deltas = betweenness_deltas(bets, ref_bets)
    costs = de_cost_function(deltas)

    # Calcular fluxo de potência
    power_flows = calc_power_flow(g, voltages, impedances)

    p = sortperm(-costs)
    costs = costs[p]
    pop_indices = pop_indices[p, :]

    if !ultima
        best_edg = pop_indices[1, :]
        pop_indices_mut = de_mutation(m, k, size(pop_indices, 1), pop_indices, beta_min, beta_max)
        pop_indices_cross = de_crossover(size(pop_indices, 1), pop_indices, pop_indices_mut, crossover_rate)
        pop_indices = pop_indices_cross
        pop_indices[1, :] = best_edg
        return nothing, nothing, nothing, nothing
    else
        global_norms = sum(abs.(deltas), dims=1)
        deltas_by_tuple = sum(abs.(deltas), dims=2)
        deltas_by_edge = zeros(Float64, m)
        d = Dict([0, 0] => 0)
        for i in 1:m
            d[edges(g)[i]] = i
        end
        for i in 1:size(ve, 1)
            for j in 1:k
                edg = ve[i, 2*j-1:2*j]
                idx = d[edg]
                deltas_by_edge[idx] += deltas_by_tuple[i]
            end
        end
        return ve, ive, DataFrame(V=1:length(global_norms), DELTA=global_norms), DataFrame(DELTA=deltas_by_edge), vec(deltas_by_tuple), power_flows
    end
end

#  Atualizando a função principal cfb_de incluir novos parametros de impedncia e tensao
function cfb_de(g::Graph, k::Integer, pop_size::Integer, crossover_rate::Float64, beta_min::Float64, beta_max::Float64, iter_num::Integer, voltages::Vector{Complex{Float64}}, impedances::Vector{Complex{Float64}}, arquivo_saida::String)
    ref_cfb = current_flow_betweenness(g)
    m = ne(g)
    edge_index_pop = sample_initial_pop(m, pop_size, k)
    for i in 1:iter_num-1
        de_iter!(g, edge_index_pop, ref_cfb, crossover_rate, beta_min, beta_max, impedances, voltages)
    end
    ve, ive, global_norms_df, edge_norms_df, local_norms, power_flows = de_iter!(g, edge_index_pop, ref_cfb, crossover_rate, beta_min, beta_max, impedances, voltages, true)
    global_norms_vector = Vector{Float64}(global_norms_df[:, "DELTA"])
    edge_norms_vector = Vector{Float64}(edge_norms_df[:, "DELTA"])
    export_de_results(g, ve, global_norms_vector, edge_norms_vector, local_norms, ive, k, arquivo_saida, "de")

    # Exportar fluxos de potência
    CSV.write("power_flows.csv", Tables.table(power_flows), writeheader=false)
end


