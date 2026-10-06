#!/usr/bin/env julia
# FlashWeave bootstrap worker used by enbi_analysis.R.
# Input: samples x taxa TSV. Output: bootstrap, rho, n_edges, n_samples.

using FlashWeave
using DelimitedFiles
using Random

function parse_args()
    length(ARGS) == 5 || error("expected 5 args: in.tsv out.tsv B fraction seed")
    return (
        ARGS[1],
        ARGS[2],
        parse(Int, ARGS[3]),
        parse(Float64, ARGS[4]),
        parse(Int, ARGS[5])
    )
end

function rho_from_edgelist(path)
    (isfile(path) && filesize(path) > 0) || return (NaN, 0)
    edges = readdlm(path, '\t')
    size(edges, 2) >= 3 || return (NaN, 0)

    weights = Float64[]
    for value in edges[:, 3]
        if value isa Number
            push!(weights, Float64(value))
        else
            parsed = tryparse(Float64, strip(string(value)))
            parsed === nothing || push!(weights, parsed)
        end
    end

    isempty(weights) && return (NaN, 0)
    denominator = sum(abs.(weights))
    return (
        denominator > 0 ? sum(weights) / denominator : NaN,
        length(weights)
    )
end

function main()
    infile, outfile, n_bootstrap, fraction, seed = parse_args()
    raw = readdlm(infile, '\t')
    data = Float64.(raw[2:end, 2:end])
    n_samples = size(data, 1)
    subsample_n = min(max(3, round(Int, fraction * n_samples)), n_samples)
    rng = MersenneTwister(seed)

    println(
        "group table: $infile | samples=$n_samples taxa=$(size(data, 2)) " *
        "| B=$n_bootstrap fraction=$fraction subsample_n=$subsample_n"
    )

    open(outfile, "w") do io
        println(io, "bootstrap\trho\tn_edges\tn_samples")
        for bootstrap in 1:n_bootstrap
            selected = randperm(rng, n_samples)[1:subsample_n]
            subset = data[selected, :]
            present = vec(sum(subset, dims = 1) .> 0)
            subset = subset[:, present]

            rho = NaN
            n_edges = 0
            try
                network = learn_network(
                    subset;
                    sensitive = true,
                    heterogeneous = false,
                    verbose = false
                )
                edge_path = tempname() * ".edgelist"
                save_network(edge_path, network)
                rho, n_edges = rho_from_edgelist(edge_path)
                rm(edge_path, force = true)
            catch error
                @warn "bootstrap failed" bootstrap error
            end

            println(io, "$bootstrap\t$rho\t$n_edges\t$subsample_n")
            flush(io)
        end
    end

    println("wrote $outfile")
end

main()
