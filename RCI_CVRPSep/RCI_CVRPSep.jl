using CVRPSEP, Test, CSV, DataFrames, JSON

function read_Kth_data(data,k)
    customers::Int64 = data[k,"customers"]
    
    K::Int64 = data[k,"vehicles"]
    
    C::Int64 = data[k,"capacity"]
     
    edges_raw::Vector{Vector{Int64}} = Vector{Vector{Int64}}(JSON.parse(data[k,"edges"]))
    
    xij_raw::Vector{Float64} = Vector{Float64}(JSON.parse(data[k,"x_bar"]))
    
    demands::Vector{Int64} = Vector{Int64}(JSON.parse(data[k,"demand"]))

    edges::Array{Int64, 2} = Array{Int64, 2}(undef,Int64(length(edges_raw)/2), 2)

    xij:: Vector{Float64} = Vector{Float64}(undef,Int64(length(xij_raw)/2))

    index::Int64 = 1

    for i in 1:2:length(edges_raw)  # get undirected support graph s.t Xij >0
        edges[index,1] = edges_raw[i][1]
        edges[index,2] = edges_raw[i][2]
        xij[index] = xij_raw[i]
        index = index + 1
    end

    input = Vector{Any}([customers, K, C, demands, edges, xij])
    
    return input
end

#=function testdata_provider(input)
    n_customers = input[1] + 1
    demand = input[4]
    capacity = input[3]
    edge_tail = input[5][:,2]
    edge_head = input[5][:,1]
    edge_x = input[6]
    n_edges = length(edge_head)
    max_n_cuts = 10
    integrality_tolerance = 1e-6


    return demand, capacity, edge_head, edge_tail, edge_x
end=#

data = CSV.read("test.csv",DataFrame)

@testset "CVRPSEP.jl" begin
    # Write your tests here.
    fail = 0
    for i in 1:size(data)[1]
        input = read_Kth_data(data,i)
        demand, capacity, edge_head, edge_tail, edge_x = input[4], input[3], input[5][:,1], input[5][:,2],input[6]
        cut_manager = CutManager()

        idx = findall(v -> v > 0.0, edge_x)
        edge_tail = edge_tail[idx]
        edge_head = edge_head[idx]
        edge_x = edge_x[idx]

        S, RHS = rounded_capacity_inequalities!(
            cut_manager, 
            demand, 
            capacity, 
            edge_tail, 
            edge_head, 
            edge_x,
            integrality_tolerance = 1e-6,
            max_n_cuts = 6000
        )
        if isempty(S)
            fail = fail + 1
        end

        @show S
        @show RHS
    end
    #println(fail/100)

end
