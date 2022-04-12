using Eirene
using NPZ
using CSV
using DataFrames
print("Finished importing. \n")
paths = ["0A.npy", "0B.npy", "2A.npy", "2B.npy", "5A.npy", "5B.npy", "8A.npy", "8B.npy"]
paths_copy = []
for path in paths
    m = npzread("/gss/home/a.kuznetsov/NoNAN_DM_new/$path");
        print("NPZ read is done for $path \n")
    C = eirene(m, model="vr", maxdim=2, maxrad=0.9)
    print("Eirene done for $path \n")
    #form DataFrame
    df = DataFrame(Replica = String[],
        Class = Int64[], Dim = Int64[], 
        Birth = Float64[], Death = Float64[], Lifetime = Float64[], 
        Numvert = Int64[], Range = Int64[],  Vertices = Array{Int64,1}[],
        Unsorted_Vertices = Array{Int64,1}[])
    print("dataframe is  done for $path \n")

#Replica number (take from filename)
    rep = String(match(r"\d+.", path).match)
    
    maxdim=2
    for d in 1:maxdim
            bar = barcode(C, dim=d)
                    print("barcode done for $path \n")
            dimsize = size(bar)[1]
            #n_cyclenumber = size(C["cyclerep"][d+2])[1]
            print("dimsize done for $path and $d \n")

            #  Dataframe: chr | class | dim | birth | death | lifetime | cycle rep
            #chrs = [chr for i in 1:dimsize]
            #segments = [segment for i in 1:dimsize]
            replicas = [rep for i in 1:dimsize]
            print("replicas done for $path and $d \n")
            classes = [i for i in 1:dimsize]
            print("classes done for $path and $d \n")
            dims = [d for i in 1:dimsize]
            print("dims done for $path and $d \n")
            births = bar[:,1]
            print("births done for $path and $d \n")
            deaths = bar[:,2]
            print("deaths done for $path and $d \n")
            lifetimes = bar[:,2] - bar[:,1]
            print("lifetime done for $path and $d \n")
            #extract only unique vertices of representative cycles
            unique_vertices = [sort(unique(vec(classrep(C, class=cl, dim=d)))) for cl in 1:dimsize]
            print("unique vert is  done for $path and $d \n")
            vertices = [vec(classrep(C, class=cl, dim=d)) for cl in 1:dimsize]

            numvert = [size(u)[1] for u in unique_vertices]
            print("numvert done for $path and $d \n")

            range = [maximum(u) - minimum(u) for u in unique_vertices]
            print("range done for $path and $d \n")


            tmpdf = DataFrame(  Replica = replicas,
                                Class = classes, Dim = dims, 
                                Birth = births, Death = deaths, Lifetime = lifetimes, 
                                Numvert = numvert, Range=range, Vertices = unique_vertices, Unsorted_Vertices = vertices)
            print("tmpdf done for $path and $d \n")

            append!(df, tmpdf)
            print("append is  done for $path and $d \n")
        end
#save results to DataFrame
    savename = "/gss/home/a.kuznetsov/results_corr/$rep.csv"
    CSV.write(savename, df)
    print("save file: $savename \n")
    end