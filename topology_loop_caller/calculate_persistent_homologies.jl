using ArgParse
using CSV
using DataFrames
using Eirene
using Logging
using NPZ

function get_file_extention(filename)
    return filename[findlast(isequal('.'), filename):end]
end

"""
    parse_commandline()

Parses arguments for CLI Topological Data Analysis step.
"""
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--matrices-path"
        help = "A path to a folder with distance matrix/matrices (only)"
        arg_type = String
        default = "../results/distance_matrices/"
        "--results-path"
        help = "A path to a folder where resulting CSV will be saved"
        arg_type = String
        default = "../results/persistent_homology_results/"
        "--maxdim"
        help = "Compute persistent homology in dimensions 0, ..., k."
        arg_type = Int
        default = 2
        "--minrad"
        help = "Compute homology from time t onward."
        arg_type = Float64
        default = -Inf
        "--maxrad"
        help = "Stop computing homology after time t."
        arg_type = Float64
        default = Inf
        "--numrad"
        help = "Int or Inf, divide the interval from minrad to maxrad into N equally spaced steps, and compute the homology of each step. If the value of numrad is set to Inf, then homology will computed at every time point."
        arg_type = Float64
        default = Inf
        "--model"
        help = "Used Eirene model, 'pc' (point cloud), 'vr' (vietoris-rips), or 'complex'."
        arg_type = String
        default = "vr"
        "--zero-order-homologies-skip"
        help = "Whether to skip zero order homologies."
        arg_type = Bool
        default = true
    end
    return parse_args(s)
end

@info "Finished importing."

function main()
    @show parsed_args = parse_commandline()

    matrices_path = parsed_args["matrices-path"]
    results_path = parsed_args["results-path"]
    maxdim = parsed_args["maxdim"]
    minrad = parsed_args["minrad"]
    maxrad = parsed_args["maxrad"]
    numrad = parsed_args["numrad"]
    model = parsed_args["model"]
    zero_order_homologies_skip = parsed_args["zero-order-homologies-skip"]

    # Creating a folder with results:
    if .!isdir(results_path)
        mkpath(results_path)
        @info "Created path $results_path for results."
    end

    if zero_order_homologies_skip
        dimensions_range = 1:maxdim
    else
        dimensions_range = 0:maxdim
    end

    # Parsing all files inside the folder
    paths = readdir(matrices_path, join=true)
    # Filter filenames, find .npy extentions:
    paths = filter(x -> get_file_extention(x) == ".npy", paths)
    @info "Starting to parse the following files: $paths"
    for path in paths
        m = npzread(path)
        @info "NPZ read is done for $path"
        C = eirene(m, model=model, maxdim=maxdim, minrad=minrad, maxrad=maxrad, numrad=numrad)
        @info "Eirene done for $path"
        # Form a DataFrame
        df = DataFrame(Replica=String[],
            Class=Int64[], Dim=Int64[],
            Birth=Float64[], Death=Float64[], Lifetime=Float64[],
            Numvert=Int64[], Range=Int64[], Vertices=Array{Int64,1}[],
            Unsorted_Vertices=Array{Int64,1}[])

        # Filename / Replica (taken from filename)
        rep = String(match(r"[^\/]+(?=\.[^\/.]*$)", path).match)
        for d in dimensions_range
            bar = barcode(C, dim=d)
            dimsize = size(bar)[1]
            #  Dataframe: chr | class | dim | birth | death | lifetime | cycle rep
            replicas = [rep for i in 1:dimsize]
            classes = [i for i in 1:dimsize]
            dims = [d for i in 1:dimsize]
            births = bar[:, 1]
            deaths = bar[:, 2]
            lifetimes = bar[:, 2] - bar[:, 1]
            # Extract only unique vertices of representative cycles
            unique_vertices = [sort(unique(vec(classrep(C, class=cl, dim=d)))) for cl in 1:dimsize]
            vertices = [vec(classrep(C, class=cl, dim=d)) for cl in 1:dimsize]
            numvert = [size(u)[1] for u in unique_vertices]
            range = [maximum(u) - minimum(u) for u in unique_vertices]
            tmpdf = DataFrame(Replica=replicas,
                Class=classes, Dim=dims,
                Birth=births, Death=deaths, Lifetime=lifetimes,
                Numvert=numvert, Range=range, Vertices=unique_vertices, Unsorted_Vertices=vertices)
            @info "tmpdf created for $path and $d"
            append!(df, tmpdf)
        end
        # Save results to DataFrame
        @info "DF for $rep is created."
        savename = joinpath(results_path, "$rep.csv")
        CSV.write(savename, df)
        @info " Saved file: $savename"
    end
end

main()
