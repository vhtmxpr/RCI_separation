using CSV, DataFrames, JSON, SparseArrays, CVRPSEP

function read_Kth_data(data, k)

    customers::Int64 = data[k,"customers"]
    
    K::Int64 = data[k,"vehicles"]
    
    C::Int64 = data[k,"capacity"]
     
    edges_raw::Vector{Vector{Int64}} = Vector{Vector{Int64}}(JSON.parse(data[k,"edges"]))
    
    xij_raw::Vector{Float64} = Vector{Float64}(JSON.parse(data[k,"x_bar"]))
    
    demands::Vector{Int64} = Vector{Int64}(JSON.parse(data[k,"demand"]))

    edges::Array{Int64, 2} = Array{Int64, 2}(undef,Int64(length(edges_raw)/2), 2)

    xij::Array{Float64, 2} = zeros(Int64,data[k,"customers"] + 1, data[k,"customers"] + 1)
    
    xij_for_CVRPSep::Vector{Float64} = Vector{Float64}(undef,Int64(length(xij_raw)/2))

    index::Int64 = 1
    
    for i in 1:2:length(edges_raw)  # get undirected support graph s.t Xij >0
        if xij_raw[i] > 0.0 # Xij = 0인 edge는 가져올 필요가 없음
            edges[index,1] = edges_raw[i][1]
            edges[index,2] = edges_raw[i][2]
            xij[edges_raw[i][1], edges_raw[i][2]] = xij_raw[i]
            xij_for_CVRPSep[index] = xij_raw[i]
            index = index + 1
        end
    end

    input = Vector{Any}([customers, K, C, demands, edges, sparse(xij),xij_for_CVRPSep])
    
    return input
end

function support_graph(input::Vector{Any})
    
    g::Dict{String, Vector{Int64}} = Dict{String, Vector{Int64}}()
    
    for i in 1:input[1]+1
        g["$i"] = []
    end
    
    for i in 1:size(input[5])[1]
        push!(g["$(input[5][i,1])"],input[5][i,2])
        push!(g["$(input[5][i,2])"],input[5][i,1])
    end

    return g
end

function initial_sol_random(customer::Int64, demand::Vector{Int64}, edges::Array{Int64, 2} , xij::SparseMatrixCSC{Float64, Int64}, C::Int64)
    sol_initial_y::Array{Int64,2} = Array{Int64,2}(undef,customer+1,2)
    sum_demand::Int64 = 0

    sol_initial_y[1,:] = [0,0]
    num::Int64 = 0

    for i in 2:customer + 1 
        num = rand(0:1)
        sol_initial_y[i,:] = [num, demand[i]]
        if num == 1
            sum_demand = sum_demand + demand[i]
        end
    end
    
    sol_initial_w::Array{Int64,2} = zeros(Int64,customer+1,customer+1)
    value::Int64 = 0
    Obj_value::Float64 = 0.0
    neighbor_set::Dict{String, Int64} =Dict{String, Int64}()

    for i in 1:size(edges)[1]
        if sol_initial_y[edges[i,1],1] == 1 && sol_initial_y[edges[i,2],1] == 0
            value = 1
        elseif sol_initial_y[edges[i,1],1] == 0 && sol_initial_y[edges[i,2],1] == 1
            value = 1
        else
            value = -1
        end
        sol_initial_w[edges[i,1],edges[i,2]] = value
        if value == 1
            if edges[i,1] != 1
                if !haskey(neighbor_set, "$(edges[i,1])")
                    neighbor_set["$(edges[i,1])"] = 1
                else
                    neighbor_set["$(edges[i,1])"] = neighbor_set["$(edges[i,1])"] + 1
                end
            end
            if edges[i,2] != 1
                if !haskey(neighbor_set, "$(edges[i,2])") 
                    neighbor_set["$(edges[i,2])"] = 1
                else
                    neighbor_set["$(edges[i,2])"] = neighbor_set["$(edges[i,2])"]+ 1
                end
            end 
            Obj_value = Obj_value + xij[edges[i,1],edges[i,2]]*value
        end
    end

    Obj_value = Obj_value - 2*ceil(sum_demand/C)
    
    return Vector{Any}([sol_initial_y,sparse(sol_initial_w), sum_demand, Obj_value,neighbor_set])
end

function initial_sol_max_subset(customer::Int64, demand::Vector{Int64}, edges::Array{Int64, 2} , xij::SparseMatrixCSC{Float64, Int64}, C::Int64)

    sol_initial_y::Array{Int64,2} = Array{Int64,2}(undef,customer+1,2)
    sum_demand::Int64 = sum(demand)
    
    sol_initial_y[1,:] = [0,0]
    for i in 2:customer + 1 
        sol_initial_y[i,:] = [1, demand[i]]
    end

    sol_initial_w::Array{Int64,2} = zeros(Int64,customer+1,customer+1)
    value::Int64 = 0
    Obj_value::Float64 = 0.0
    neighbor_set::Dict{String, Int64} =Dict{String, Int64}()

    for i in 1:size(edges)[1]
        if sol_initial_y[edges[i,1],1] == 1 && sol_initial_y[edges[i,2],1] == 0
            value = 1
        elseif sol_initial_y[edges[i,1],1] == 0 && sol_initial_y[edges[i,2],1] == 1
            value = 1
        else
            value = -1
        end
        sol_initial_w[edges[i,1],edges[i,2]] = value
        if value == 1
            if edges[i,1] != 1
                if !haskey(neighbor_set, "$(edges[i,1])")
                    neighbor_set["$(edges[i,1])"] = 1
                else
                    neighbor_set["$(edges[i,1])"] = neighbor_set["$(edges[i,1])"] + 1
                end
            end
            if edges[i,2] != 1
                if !haskey(neighbor_set, "$(edges[i,2])") 
                    neighbor_set["$(edges[i,2])"] = 1
                else
                    neighbor_set["$(edges[i,2])"] = neighbor_set["$(edges[i,2])"]+ 1
                end
            end 
            Obj_value = Obj_value + xij[edges[i,1],edges[i,2]]*value
        end
    end

    Obj_value = Obj_value - 2*ceil(sum_demand/C)

    return Vector{Any}([sol_initial_y,sparse(sol_initial_w), sum_demand, Obj_value,neighbor_set])
 
end

function initial_sol_one_subset(customer::Int64, demand::Vector{Int64}, edges::Array{Int64, 2} , xij::SparseMatrixCSC{Float64, Int64}, C::Int64, g::Dict{String, Vector{Int64}})
    
    sol_initial_y::Array{Int64,2} = Array{Int64,2}(undef,customer+1,2)
    sol_initial_y[1,:] = [0,0]
    for i in 2:customer + 1 
        sol_initial_y[i,:] = [0, demand[i]]
    end
    index_max_neighbor::Int64 = 1
    max_neighbor::Int64 = -1
    for (key, value) in g
        if key == "1"
            continue
        end
        num_neighbor = length(value)
        if num_neighbor > max_neighbor
            max_neighbor = num_neighbor
            index_max_neighbor = parse(Int64,key)
        end
    end
    sol_initial_y[index_max_neighbor,1] = 1
    sum_demand::Int64 = demand[index_max_neighbor]
   
    sol_initial_w::Array{Int64,2} = zeros(Int64,customer+1,customer+1)
    neighbor_set::Dict{String, Int64} =Dict{String, Int64}()
    value::Int64 = 0
    Obj_value::Float64 = 0.0

    for i in 1:size(edges)[1]
        if sol_initial_y[edges[i,1],1] == 1 && sol_initial_y[edges[i,2],1] == 0
            value = 1
        elseif sol_initial_y[edges[i,1],1] == 0 && sol_initial_y[edges[i,2],1] == 1
            value = 1
        else
            value = -1
        end
        sol_initial_w[edges[i,1],edges[i,2]] = value
        if value == 1
            if edges[i,1] != 1
                if !haskey(neighbor_set, "$(edges[i,1])")
                    neighbor_set["$(edges[i,1])"] = 1
                else
                    neighbor_set["$(edges[i,1])"] = neighbor_set["$(edges[i,1])"] + 1
                end
            end
            if edges[i,2] != 1
                if !haskey(neighbor_set, "$(edges[i,2])") 
                    neighbor_set["$(edges[i,2])"] = 1
                else
                    neighbor_set["$(edges[i,2])"] = neighbor_set["$(edges[i,2])"]+ 1
                end
            end 
            Obj_value = Obj_value + xij[edges[i,1],edges[i,2]]*value
        end
    end

    Obj_value = Obj_value - 2*ceil(sum_demand/C)

    return Vector{Any}([sol_initial_y,sparse(sol_initial_w), sum_demand, Obj_value,neighbor_set])
 
end

function initial_sol_CVRPSep(demands::Vector{Int64}, edges::Array{Int64, 2}, xij_for_CVRPSep::Vector{Float64}, C::Int64, customer::Int64, xij::SparseMatrixCSC{Float64, Int64})

    demand = demands
    capacity = C
    edge_head = edges[:,1]
    edge_tail = edges[:,2]
    edge_x =  xij_for_CVRPSep
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
        max_n_cuts = 1
    )
    #@show S[1]
    if !isempty(S)
        sol_cvrp::Vector{Int64} = S[1]
        sol_initial_y::Array{Int64,2} = Array{Int64,2}(undef,customer+1,2)
        sum_demand::Int64 = 0
        index::Int64 = 1
        len::Int64 = length(sol_cvrp)

        for i in 1:length(demands)
            if index <= len
                if i == sol_cvrp[index] 
                    sol_initial_y[i,:] = [1,demands[i]]
                    sum_demand = sum_demand + demands[i]
                    index = index + 1
                    continue
                end
            end
            sol_initial_y[i,:] = [0,demands[i]]

        end

        sol_initial_w::Array{Int64,2} = zeros(Int64,customer+1,customer+1)
        value::Int64 = 0
        Obj_value::Float64 = 0.0
        neighbor_set::Dict{String, Int64} =Dict{String, Int64}()
        
        for i in 1:size(edges)[1]
            if sol_initial_y[edges[i,1],1] == 1 && sol_initial_y[edges[i,2],1] == 0
                value = 1
            elseif sol_initial_y[edges[i,1],1] == 0 && sol_initial_y[edges[i,2],1] == 1
                value = 1
            else
                value = -1
            end
            sol_initial_w[edges[i,1],edges[i,2]] = value
            if value == 1
                if edges[i,1] != 1
                    if !haskey(neighbor_set, "$(edges[i,1])")
                        neighbor_set["$(edges[i,1])"] = 1
                    else
                        neighbor_set["$(edges[i,1])"] = neighbor_set["$(edges[i,1])"] + 1
                    end
                end
                if edges[i,2] != 1
                    if !haskey(neighbor_set, "$(edges[i,2])") 
                        neighbor_set["$(edges[i,2])"] = 1
                    else
                        neighbor_set["$(edges[i,2])"] = neighbor_set["$(edges[i,2])"] + 1
                    end
                end 
                Obj_value = Obj_value + xij[edges[i,1],edges[i,2]]*value
            end
        end

        Obj_value = Obj_value - 2*ceil(sum_demand/C)

        return Vector{Any}([sol_initial_y,sparse(sol_initial_w), sum_demand, Obj_value,neighbor_set])
    else
        return 0
    end
end

function initial_sol_CVRPSep_best(demands::Vector{Int64}, edges::Array{Int64, 2}, xij_for_CVRPSep::Vector{Float64}, C::Int64, customer::Int64, xij::SparseMatrixCSC{Float64, Int64})

    demand = demands
    capacity = C
    edge_head = edges[:,1]
    edge_tail = edges[:,2]
    edge_x =  xij_for_CVRPSep
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
    #@show S

    if !isempty(S)
        Obj_value_best::Float64 = 100.0
        S_best::Vector{Int64} = []
        for j in 1:length(S)
            demand_raw_sum::Int64 = 0
            Obj_value_raw::Float64 = 0.0
            for k in 1: length(S[j])
                demand_raw_sum = demand_raw_sum + demands[S[j][k]]
            end
            for i in 1:size(edges)[1]
                if (edges[i,1] in S[j]) && !(edges[i,2] in S[j])
                    value = 1
                elseif !(edges[i,1] in S[j]) && (edges[i,2] in S[j])
                    value = 1
                else
                    value = 0
                end
                Obj_value_raw = Obj_value_raw + xij[edges[i,1],edges[i,2]]*value
            end
            Obj_value_raw = Obj_value_raw - 2*ceil(demand_raw_sum/C)
            if Obj_value_best > Obj_value_raw 
                Obj_value_best = Obj_value_raw
                S_best = S[j]
            end
        end
    end

    if !isempty(S)
        sol_cvrp::Vector{Int64} = S_best
        sol_initial_y::Array{Int64,2} = Array{Int64,2}(undef,customer+1,2)
        sum_demand::Int64 = 0
        index::Int64 = 1
        len::Int64 = length(sol_cvrp)

        for i in 1:length(demands)
            if index <= len
                if i == sol_cvrp[index] 
                    sol_initial_y[i,:] = [1,demands[i]]
                    sum_demand = sum_demand + demands[i]
                    index = index + 1
                    continue
                end
            end
            sol_initial_y[i,:] = [0,demands[i]]
        end

        sol_initial_w::Array{Int64,2} = zeros(Int64,customer+1,customer+1)
        value::Int64 = 0
        Obj_value::Float64 = 0.0
        neighbor_set::Dict{String, Int64} =Dict{String, Int64}()
        
        for i in 1:size(edges)[1]
            if sol_initial_y[edges[i,1],1] == 1 && sol_initial_y[edges[i,2],1] == 0
                value = 1
            elseif sol_initial_y[edges[i,1],1] == 0 && sol_initial_y[edges[i,2],1] == 1
                value = 1
            else
                value = -1
            end
            sol_initial_w[edges[i,1],edges[i,2]] = value
            if value == 1
                if edges[i,1] != 1
                    if !haskey(neighbor_set, "$(edges[i,1])")
                        neighbor_set["$(edges[i,1])"] = 1
                    else
                        neighbor_set["$(edges[i,1])"] = neighbor_set["$(edges[i,1])"] + 1
                    end
                end
                if edges[i,2] != 1
                    if !haskey(neighbor_set, "$(edges[i,2])") 
                        neighbor_set["$(edges[i,2])"] = 1
                    else
                        neighbor_set["$(edges[i,2])"] = neighbor_set["$(edges[i,2])"] + 1
                    end
                end 
                Obj_value = Obj_value + xij[edges[i,1],edges[i,2]]*value
            end
        end

        Obj_value = Obj_value - 2*ceil(sum_demand/C)

        return Vector{Any}([sol_initial_y,sparse(sol_initial_w), sum_demand, Obj_value,neighbor_set])
    else
        return 0
    end
end

function stop(T::Float64)

    stop::Int64 = 0
    
    if(T < 0.01)
        stop = 1
    end

    return stop
end

function neighbor_sol_demand(y_cur::Array{Int64,2}, y_cur_copy::Array{Int64,2}, w_cur::SparseMatrixCSC{Int64, Int64}, w_cur_copy::SparseMatrixCSC{Int64, Int64}, g::Dict{String, Vector{Int64}}, demand_cur::Int64, Obj_cur::Float64, xij::SparseMatrixCSC{Float64, Int64}, C::Int64, T_cur::Float64,neighbor_set::Dict{String, Int64},in_delta_S::Int64) 
    neighbor_candidate::Vector{Int64} = []
    neighbor_candidate2::Vector{Int64} = []
    neighbor_candidate3::Vector{Int64} = []
    for i in 1: size(y_cur)[1]
        if y_cur[i,1]==0
            if ceil((demand_cur+ y_cur[i,2])/C) > ceil((demand_cur)/C)
                push!(neighbor_candidate,i)
            elseif  ceil((demand_cur+ 2*y_cur[i,2])/C) > ceil((demand_cur)/C)
                push!(neighbor_candidate2,i)
            elseif ceil((demand_cur+ 3*y_cur[i,2])/C) > ceil((demand_cur)/C)
                push!(neighbor_candidate3,i)
            else
            end
        end
    end

    index_fliped1::Int64 = -1
    index_fliped2::Int64 = -1
    index_fliped3::Int64 = -1
    result::Int64 = 0
    Obj_new::Float64 = Obj_cur
    demand_new::Int64 = demand_cur
    neighbor_list::Vector{Int64} = []
    neighbor_list2::Vector{Int64} = []
    neighbor_list3::Vector{Int64} = []

 ##############
    #=case::Int64 = rand(1:3)
    check1::Int64 = 0
    check2::Int64 = 0
    check3::Int64 = 0
    while check1 == 1 && check2 == 1 && check3 == 1
        if case == 1
            if !isempty(neighbor_candidate)
                index_fliped1 = neighbor_candidate[rand(1:length(neighbor_candidate))]
                result = 1
                break
            else
                case = rand(2:3)
                check1 = 1
            end
        elseif case == 2
            if length(neighbor_candidate2) >= 2
                index_fliped1 = neighbor_candidate2[rand(1:length(neighbor_candidate2))]
                index_fliped2 = index_fliped1
                while index_fliped2 == index_fliped1
                    index_fliped2 = neighbor_candidate2[rand(1:length(neighbor_candidate2))]
                end
                result = 2
                break
            else
                case = rand([1,3])
                check2 = 1
            end
        else
            if length(neighbor_candidate3) >= 3
                index_fliped1 = neighbor_candidate3[rand(1:length(neighbor_candidate3))]
                index_fliped2 = index_fliped1
                while index_fliped2 == index_fliped1
                    index_fliped2 = neighbor_candidate3[rand(1:length(neighbor_candidate3))]
                end
                index_fliped3 = index_fliped2
                while index_fliped3 == index_fliped2 || index_fliped3 == index_fliped1
                    index_fliped3 = neighbor_candidate3[rand(1:length(neighbor_candidate3))]
                end
                result = 3
                break
            else
                case = rand(1:2)
                check3 = 1
            end
        end
    end=#
 ################
    if !isempty(neighbor_candidate)
        index_fliped1 = neighbor_candidate[rand(1:length(neighbor_candidate))]
        result = 1
    else
        #return [0.0,0,0]
        if length(neighbor_candidate2) >= 2
            index_fliped1 = neighbor_candidate2[rand(1:length(neighbor_candidate2))]
            index_fliped2 = index_fliped1
            while index_fliped2 == index_fliped1
                index_fliped2 = neighbor_candidate2[rand(1:length(neighbor_candidate2))]
            end
            result = 2
        else
            return [0.0,0,0]
            if length(neighbor_candidate3) >= 3
                index_fliped1 = neighbor_candidate3[rand(1:length(neighbor_candidate3))]
                index_fliped2 = index_fliped1
                while index_fliped2 == index_fliped1
                    index_fliped2 = neighbor_candidate3[rand(1:length(neighbor_candidate3))]
                end
                index_fliped3 = index_fliped2
                while index_fliped3 == index_fliped2 || index_fliped3 == index_fliped1
                    index_fliped3 = neighbor_candidate3[rand(1:length(neighbor_candidate3))]
                end
                result = 3
            else
                result = 4
            end
        end
    end
 ####################
    if result == 1||result ==2||result == 3
        if result == 1
            y_cur_copy[index_fliped1,1] = 1
            demand_new = demand_new + y_cur_copy[index_fliped1,2]
            neighbor_list = g["$index_fliped1"]
        elseif result == 2
            y_cur_copy[index_fliped1,1] = 1
            y_cur_copy[index_fliped2,1] = 1
            demand_new = demand_new + y_cur_copy[index_fliped1,2] + y_cur_copy[index_fliped2,2]
            neighbor_list = g["$index_fliped1"]
            neighbor_list2 = g["$index_fliped2"]
        else
            y_cur_copy[index_fliped1,1] = 1
            y_cur_copy[index_fliped2,1] = 1
            y_cur_copy[index_fliped3,1] = 1
            demand_new = demand_new + y_cur_copy[index_fliped1,2] + y_cur_copy[index_fliped2,2] +  y_cur_copy[index_fliped3,2]
            neighbor_list = g["$index_fliped1"]
            neighbor_list2 = g["$index_fliped2"]
            neighbor_list3 = g["$index_fliped3"]
        end
        Obj_new = Obj_cur + 2*ceil(demand_cur/C) - 2*ceil(demand_new/C)


        len1::Int64 = length(neighbor_list)
        len2::Int64 = length(neighbor_list2)
        len3::Int64 = length(neighbor_list3)
        len::Int64 = len1 + len2 + len3
        edges_changed::Array{Int64,2} = Array{Int64,2}(undef,len,2)

        for i in 1:len1
            if w_cur_copy[index_fliped1,neighbor_list[i]] != 0
                if  w_cur_copy[index_fliped1,neighbor_list[i]] == 1
                    w_cur_copy[index_fliped1,neighbor_list[i]] = -1
                    Obj_new = Obj_new - xij[index_fliped1,neighbor_list[i]]
                    edges_changed[i,:] = [index_fliped1,neighbor_list[i]]
                else
                    w_cur_copy[index_fliped1,neighbor_list[i]] = 1
                    Obj_new = Obj_new + xij[index_fliped1,neighbor_list[i]]
                    edges_changed[i,:] = [index_fliped1,neighbor_list[i]]
                end
            else
                if w_cur_copy[neighbor_list[i],index_fliped1] == 1
                    w_cur_copy[neighbor_list[i],index_fliped1] = -1
                    Obj_new = Obj_new - xij[neighbor_list[i],index_fliped1]
                    edges_changed[i,:] = [neighbor_list[i],index_fliped1]
                    
                else
                    w_cur_copy[neighbor_list[i],index_fliped1] = 1
                    Obj_new = Obj_new + xij[neighbor_list[i],index_fliped1]
                    edges_changed[i,:] = [neighbor_list[i],index_fliped1]
                end
            end
        end
        if result != 1
            for i in 1:len2
                if w_cur_copy[index_fliped2,neighbor_list2[i]] != 0
                    if  w_cur_copy[index_fliped2,neighbor_list2[i]] == 1
                        w_cur_copy[index_fliped2,neighbor_list2[i]] = -1
                        Obj_new = Obj_new - xij[index_fliped2,neighbor_list2[i]]
                        edges_changed[i+len1,:] = [index_fliped2,neighbor_list2[i]]
                    else
                        w_cur_copy[index_fliped2,neighbor_list2[i]] = 1
                        Obj_new = Obj_new + xij[index_fliped2,neighbor_list2[i]]
                        edges_changed[i+len1,:] = [index_fliped2,neighbor_list2[i]]
                    end
                else
                    if w_cur_copy[neighbor_list2[i],index_fliped2] == 1
                        w_cur_copy[neighbor_list2[i],index_fliped2] = -1
                        Obj_new = Obj_new - xij[neighbor_list2[i],index_fliped2]
                        edges_changed[i+len1,:] = [neighbor_list2[i],index_fliped2]
                        
                    else
                        w_cur_copy[neighbor_list2[i],index_fliped2] = 1
                        Obj_new = Obj_new + xij[neighbor_list2[i],index_fliped2]
                        edges_changed[i+len1,:] = [neighbor_list2[i],index_fliped2]
                    end
                end
            end
        end
        if result == 3
            for i in 1:len3
                if w_cur_copy[index_fliped3,neighbor_list3[i]] != 0
                    if  w_cur_copy[index_fliped3,neighbor_list3[i]] == 1
                        w_cur_copy[index_fliped3,neighbor_list3[i]] = -1
                        Obj_new = Obj_new - xij[index_fliped3,neighbor_list3[i]]
                        edges_changed[i+len1+len2,:] = [index_fliped3,neighbor_list3[i]]
                    else
                        w_cur_copy[index_fliped3,neighbor_list3[i]] = 1
                        Obj_new = Obj_new + xij[index_fliped3,neighbor_list3[i]]
                        edges_changed[i+len1+len2,:] = [index_fliped3,neighbor_list3[i]]
                    end
                else
                    if w_cur_copy[neighbor_list3[i],index_fliped3] == 1
                        w_cur_copy[neighbor_list3[i],index_fliped3] = -1
                        Obj_new = Obj_new - xij[neighbor_list3[i],index_fliped3]
                        edges_changed[i+len1+len2,:] = [neighbor_list3[i],index_fliped3]
                        
                    else
                        w_cur_copy[neighbor_list3[i],index_fliped3] = 1
                        Obj_new = Obj_new + xij[neighbor_list3[i],index_fliped3]
                        edges_changed[i+len1+len2,:] = [neighbor_list3[i],index_fliped3]
                    end
                end
            end
        end
        ΔObj_value::Float64 = Obj_new - Obj_cur
    
        if ΔObj_value <= 0 || rand() < exp(-ΔObj_value/T_cur)
            y_cur[index_fliped1,1] = y_cur_copy[index_fliped1,1]
            if result != 1
                y_cur[index_fliped2,1] = y_cur_copy[index_fliped2,1]
            end
            if result == 3
                y_cur[index_fliped3,1] = y_cur_copy[index_fliped3,1]
            end
            for i in 1:len
                if w_cur[edges_changed[i,1],edges_changed[i,2]] == w_cur_copy[edges_changed[i,1],edges_changed[i,2]]  
                    continue
                end
                w_cur[edges_changed[i,1],edges_changed[i,2]] = w_cur_copy[edges_changed[i,1],edges_changed[i,2]]
                ##neighbor_sol_2flip_in_delat_S와 함께 쓸 경우에만 사용##
                if in_delta_S == 1 
                    if w_cur[edges_changed[i,1],edges_changed[i,2]] == 1
                        if edges_changed[i,1] != 1
                            if !haskey(neighbor_set, "$(edges_changed[i,1])")
                                neighbor_set["$(edges_changed[i,1])"] = 1
                            else
                                neighbor_set["$(edges_changed[i,1])"] = neighbor_set["$(edges_changed[i,1])"] + 1
                            end
                        end
                        if edges_changed[i,2] != 1
                            if !haskey(neighbor_set, "$(edges_changed[i,2])")
                                neighbor_set["$(edges_changed[i,2])"] =  1
                            else
                                neighbor_set["$(edges_changed[i,2])"] = neighbor_set["$(edges_changed[i,2])"] + 1
                            end
                        end
                    else
                        if edges_changed[i,1] != 1
                            if neighbor_set["$(edges_changed[i,1])"] == 1
                                delete!(neighbor_set, "$(edges_changed[i,1])")
                            else
                                neighbor_set["$(edges_changed[i,1])"] = neighbor_set["$(edges_changed[i,1])"] - 1
                            end
                        end
                        if edges_changed[i,2] != 1
                            if neighbor_set["$(edges_changed[i,2])"] == 1
                                delete!(neighbor_set, "$(edges_changed[i,2])")
                            else
                                neighbor_set["$(edges_changed[i,2])"] = neighbor_set["$(edges_changed[i,2])"] - 1
                            end
                        end
                    end
                end
                ##-----------------------------------------------------##
            end
            demand_cur = demand_new
            Obj_cur = Obj_new
        else
            y_cur_copy[index_fliped1,1] = y_cur[index_fliped1,1]
            if result != 1
                y_cur_copy[index_fliped2,1] = y_cur[index_fliped2,1]
            end
            if result == 3
                y_cur_copy[index_fliped3,1] = y_cur[index_fliped3,1]
            end
            for i in 1:len
                w_cur_copy[edges_changed[i,1],edges_changed[i,2]] = w_cur[edges_changed[i,1],edges_changed[i,2]]
            end
            demand_new = demand_cur
            Obj_new = Obj_cur
        end

        return Vector{Any}([Obj_cur, demand_cur,result])
    else
        return Vector{Any}([0.0,0,result])
    end


end

function neighbor_sol_hamming_ball(y_cur::Array{Int64,2}, y_cur_copy::Array{Int64,2}, w_cur::SparseMatrixCSC{Int64, Int64}, w_cur_copy::SparseMatrixCSC{Int64, Int64}, g::Dict{String, Vector{Int64}}, demand_cur::Int64, Obj_cur::Float64, xij::SparseMatrixCSC{Float64, Int64}, C::Int64, T_cur::Float64,mod::Int64,num::Int64,neighbor_set::Dict{String, Int64},in_delta_S::Int64)
    num_ball::Int64 = num
    num_element::Int64 = ceil(size(y_cur)[1]/num_ball)
    index_fliped_set::Array{Int64,2} = Array{Int64,2}(undef,num_ball,2)
    i::Int64 = 1
    index_start::Int64 = 1
    index_end::Int64 = num_element
    index_fliped::Int64 = 0
    Obj_new::Float64 = Obj_cur
    demand_new::Int64 = demand_cur
    neighbor_list::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef,num_ball)
    len::Vector{Int64} = Vector{Int64}(undef,num_ball)
    sum_len::Int64 = 0
    while i <= num_ball
        index_fliped = rand(index_start:index_end)
        if index_fliped == 1
            continue
        else
            index_fliped_set[i,1] = index_fliped
            index_fliped_set[i,2] = rand(0:1)
            if index_fliped_set[i,2] == 1
                if y_cur_copy[index_fliped,1] == 1
                    y_cur_copy[index_fliped,1] = 0 
                    demand_new = demand_new - y_cur_copy[index_fliped,2]
                else
                    y_cur_copy[index_fliped,1] = 1
                    demand_new = demand_new + y_cur_copy[index_fliped,2]
                end
                neighbor_list[i] = g["$index_fliped"]
                len[i] = length(neighbor_list[i])
                sum_len = sum_len + len[i]
            else
               neighbor_list[i] = []
               len[i] = 0
            end
            index_start = index_start + num_element
            index_end = min(index_end + num_element, size(y_cur)[1])
            i = i + 1
        end
    end

    Obj_new = Obj_cur + 2*ceil(demand_cur/C) - 2*ceil(demand_new/C)
    edges_changed::Array{Int64,2} = Array{Int64,2}(undef,sum_len,2)
    index::Int64 = 1
    change::Int64 = 0
    for i in 1:num_ball
        if index_fliped_set[i,2] == 1
            for j in 1:len[i]
                if w_cur_copy[index_fliped_set[i,1],neighbor_list[i][j]] != 0
                    if w_cur_copy[index_fliped_set[i,1],neighbor_list[i][j]] == 1
                        w_cur_copy[index_fliped_set[i,1],neighbor_list[i][j]] = -1
                        Obj_new = Obj_new - xij[index_fliped_set[i,1],neighbor_list[i][j]]
                        edges_changed[index,:] = [index_fliped_set[i,1],neighbor_list[i][j]]
                    else
                        w_cur_copy[index_fliped_set[i,1],neighbor_list[i][j]] = 1
                        Obj_new = Obj_new + xij[index_fliped_set[i,1],neighbor_list[i][j]]
                        edges_changed[index,:] = [index_fliped_set[i,1],neighbor_list[i][j]]
                    end
                else
                    if w_cur_copy[neighbor_list[i][j],index_fliped_set[i,1]] == 1
                        w_cur_copy[neighbor_list[i][j],index_fliped_set[i,1]] = -1
                        Obj_new = Obj_new - xij[neighbor_list[i][j],index_fliped_set[i,1]]
                        edges_changed[index,:] = [neighbor_list[i][j],index_fliped_set[i,1]]
                    else
                        w_cur_copy[neighbor_list[i][j],index_fliped_set[i,1]] = 1
                        Obj_new = Obj_new + xij[neighbor_list[i][j],index_fliped_set[i,1]]
                        edges_changed[index,:] = [neighbor_list[i][j],index_fliped_set[i,1]]
                    end
                end
                index = index + 1
            end
            change = 1
        end
    end
    ΔObj_value::Float64 = Obj_new - Obj_cur

    if change == 1
        if (ΔObj_value <= 0 || rand() < exp(-ΔObj_value/T_cur)) && mod == 1
            for i in 1:num_ball
                y_cur[index_fliped_set[i,1],1] = y_cur_copy[index_fliped_set[i,1],1]
            end
            for i in 1:sum_len
                if w_cur[edges_changed[i,1],edges_changed[i,2]] == w_cur_copy[edges_changed[i,1],edges_changed[i,2]]  
                    continue
                end
                w_cur[edges_changed[i,1],edges_changed[i,2]] = w_cur_copy[edges_changed[i,1],edges_changed[i,2]]
                ##neighbor을 랜덤으로 고를 때만 사용##
                if in_delta_S == 1 
                    if w_cur[edges_changed[i,1],edges_changed[i,2]] == 1
                        if edges_changed[i,1] != 1
                            if !haskey(neighbor_set, "$(edges_changed[i,1])")
                                neighbor_set["$(edges_changed[i,1])"] = 1
                            else
                                neighbor_set["$(edges_changed[i,1])"] = neighbor_set["$(edges_changed[i,1])"] + 1
                            end
                        end
                        if edges_changed[i,2] != 1
                            if !haskey(neighbor_set, "$(edges_changed[i,2])")
                                neighbor_set["$(edges_changed[i,2])"] =  1
                            else
                                neighbor_set["$(edges_changed[i,2])"] = neighbor_set["$(edges_changed[i,2])"] + 1
                            end
                        end
                    else
                        if edges_changed[i,1] != 1
                            if neighbor_set["$(edges_changed[i,1])"] == 1
                                delete!(neighbor_set, "$(edges_changed[i,1])")
                            else
                                neighbor_set["$(edges_changed[i,1])"] = neighbor_set["$(edges_changed[i,1])"] - 1
                            end
                        end
                        if edges_changed[i,2] != 1
                            if neighbor_set["$(edges_changed[i,2])"] == 1
                                delete!(neighbor_set, "$(edges_changed[i,2])")
                            else
                                neighbor_set["$(edges_changed[i,2])"] = neighbor_set["$(edges_changed[i,2])"] - 1
                            end
                        end
                    end
                end
                ##-----------------------------------------------------##
            end
            demand_cur = demand_new
            Obj_cur = Obj_new
        else
            for i in 1:num_ball
                y_cur_copy[index_fliped_set[i,1],1] = y_cur[index_fliped_set[i,1],1]
            end
            for i in 1:sum_len
                w_cur_copy[edges_changed[i,1],edges_changed[i,2]] = w_cur[edges_changed[i,1],edges_changed[i,2]]
            end
            demand_new = demand_cur
            Obj_new = Obj_cur
        end
    end

    if mod == 1
        return Vector{Any}([Obj_cur,demand_cur])
    else
        return ΔObj_value
    end

end

function neighbor_sol_n_filp(y_cur::Array{Int64,2}, y_cur_copy::Array{Int64,2}, w_cur::SparseMatrixCSC{Int64, Int64}, w_cur_copy::SparseMatrixCSC{Int64, Int64}, g::Dict{String, Vector{Int64}}, demand_cur::Int64, Obj_cur::Float64, xij::SparseMatrixCSC{Float64, Int64}, C::Int64, T_cur::Float64, mod::Int64, num_flip::Int64,neighbor_set::Dict{String, Int64},in_delta_S::Int64)
    Obj_new::Float64 = Obj_cur
    demand_new::Int64 = demand_cur
    index_fliped_set::Array{Int64, 2} = Array{Int64,2}(undef,num_flip,2)
    index_fliped::Int64 = 0
    index_candidate::Vector{Int64} = Vector{Int64}(undef,size(y_cur)[1])
    index_candidate_len = length(index_candidate)
    neighbor_list::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef,num_flip)
    len::Vector{Int64} = Vector{Int64}(undef,num_flip)
    sum_len::Int64 = 0
    i::Int64 = 1

    for i in 1:index_candidate_len
        index_candidate[i] = i
    end

    for i in 1:num_flip
        while index_fliped ==0
            index_fliped = index_candidate[rand(2:index_candidate_len)]
        end
        index_candidate[index_fliped] = 0
        index_fliped_set[i,1] = index_fliped
        if num_flip == 1
            index_fliped_set[i,2] = 1
        else
            index_fliped_set[i,2] = rand(0:1)
        end
        if index_fliped_set[i,2] == 1
            if y_cur_copy[index_fliped,1] == 1
                y_cur_copy[index_fliped,1] = 0 
                demand_new = demand_new - y_cur_copy[index_fliped,2]
            else
                y_cur_copy[index_fliped,1] = 1
                demand_new = demand_new + y_cur_copy[index_fliped,2]
            end
            neighbor_list[i] = g["$index_fliped"]
            len[i] = length(neighbor_list[i])
            sum_len = sum_len + len[i]
        else
           neighbor_list[i] = []
           len[i] = 0
        end
        index_fliped = 0
    end

    Obj_new = Obj_cur + 2*ceil(demand_cur/C) - 2*ceil(demand_new/C)
    edges_changed::Array{Int64,2} = Array{Int64,2}(undef,sum_len,2)
    index::Int64 = 1
    change::Int64 = 0
    for i in 1:num_flip
        if index_fliped_set[i,2] == 1
            for j in 1:len[i]
                if w_cur_copy[index_fliped_set[i,1],neighbor_list[i][j]] != 0
                    if w_cur_copy[index_fliped_set[i,1],neighbor_list[i][j]] == 1
                        w_cur_copy[index_fliped_set[i,1],neighbor_list[i][j]] = -1
                        Obj_new = Obj_new - xij[index_fliped_set[i,1],neighbor_list[i][j]]
                        edges_changed[index,:] = [index_fliped_set[i,1],neighbor_list[i][j]]
                    else
                        w_cur_copy[index_fliped_set[i,1],neighbor_list[i][j]] = 1
                        Obj_new = Obj_new + xij[index_fliped_set[i,1],neighbor_list[i][j]]
                        edges_changed[index,:] = [index_fliped_set[i,1],neighbor_list[i][j]]
                    end
                else
                    if w_cur_copy[neighbor_list[i][j],index_fliped_set[i,1]] == 1
                        w_cur_copy[neighbor_list[i][j],index_fliped_set[i,1]] = -1
                        Obj_new = Obj_new - xij[neighbor_list[i][j],index_fliped_set[i,1]]
                        edges_changed[index,:] = [neighbor_list[i][j],index_fliped_set[i,1]]
                    else
                        w_cur_copy[neighbor_list[i][j],index_fliped_set[i,1]] = 1
                        Obj_new = Obj_new + xij[neighbor_list[i][j],index_fliped_set[i,1]]
                        edges_changed[index,:] = [neighbor_list[i][j],index_fliped_set[i,1]]
                    end
                end
                index = index + 1
            end
            change = 1
        end
    end
    ΔObj_value::Float64 = Obj_new - Obj_cur

    if change == 1
        if (ΔObj_value <= 0 || rand() < exp(-ΔObj_value/T_cur)) && mod == 1
            for i in 1:num_flip
                y_cur[index_fliped_set[i,1],1] = y_cur_copy[index_fliped_set[i,1],1]
            end
            for i in 1:sum_len
                if w_cur[edges_changed[i,1],edges_changed[i,2]] == w_cur_copy[edges_changed[i,1],edges_changed[i,2]]  
                    continue
                end
                w_cur[edges_changed[i,1],edges_changed[i,2]] = w_cur_copy[edges_changed[i,1],edges_changed[i,2]]
                ##neighbor을 랜덤으로 고를 때만 사용##
                if in_delta_S == 1 
                    if w_cur[edges_changed[i,1],edges_changed[i,2]] == 1
                        if edges_changed[i,1] != 1
                            if !haskey(neighbor_set, "$(edges_changed[i,1])")
                                neighbor_set["$(edges_changed[i,1])"] = 1
                            else
                                neighbor_set["$(edges_changed[i,1])"] = neighbor_set["$(edges_changed[i,1])"] + 1
                            end
                        end
                        if edges_changed[i,2] != 1
                            if !haskey(neighbor_set, "$(edges_changed[i,2])")
                                neighbor_set["$(edges_changed[i,2])"] =  1
                            else
                                neighbor_set["$(edges_changed[i,2])"] = neighbor_set["$(edges_changed[i,2])"] + 1
                            end
                        end
                    else
                        if edges_changed[i,1] != 1
                            if neighbor_set["$(edges_changed[i,1])"] == 1
                                delete!(neighbor_set, "$(edges_changed[i,1])")
                            else
                                neighbor_set["$(edges_changed[i,1])"] = neighbor_set["$(edges_changed[i,1])"] - 1
                            end
                        end
                        if edges_changed[i,2] != 1
                            if neighbor_set["$(edges_changed[i,2])"] == 1
                                delete!(neighbor_set, "$(edges_changed[i,2])")
                            else
                                neighbor_set["$(edges_changed[i,2])"] = neighbor_set["$(edges_changed[i,2])"] - 1
                            end
                        end
                    end
                end
                ##-----------------------------------------------------##
            end
            demand_cur = demand_new
            Obj_cur = Obj_new
        else
            for i in 1:num_flip
                y_cur_copy[index_fliped_set[i,1],1] = y_cur[index_fliped_set[i,1],1]
            end
            for i in 1:sum_len
                w_cur_copy[edges_changed[i,1],edges_changed[i,2]] = w_cur[edges_changed[i,1],edges_changed[i,2]]
            end
            demand_new = demand_cur
            Obj_new = Obj_cur
        end
    end

    if mod == 1
        return Vector{Any}([Obj_cur,demand_cur])
    else
        return ΔObj_value
    end
end

function neighbor_sol_n_filp_in_delta_S(y_cur::Array{Int64,2}, y_cur_copy::Array{Int64,2}, w_cur::SparseMatrixCSC{Int64, Int64}, w_cur_copy::SparseMatrixCSC{Int64, Int64}, g::Dict{String, Vector{Int64}}, demand_cur::Int64, Obj_cur::Float64, xij::SparseMatrixCSC{Float64, Int64}, C::Int64, T_cur::Float64, mod::Int64, neighbor_set::Dict{String, Int64},num_flip::Int64)
    index_candidate::Vector{Int64} = Vector{Int64}(undef,length(neighbor_set))
    len_index_candidate::Int64 =length(index_candidate)
    index::Int64 = 1
    for (key,value) in neighbor_set
        index_candidate[index] = parse(Int64,key)
        index = index + 1
    end

    #=index_candidate::Vector{Int64} = []
    index::Int64 = 1
    for (key,value) in neighbor_set
        for i in 1:value
            push!(index_candidate,parse(Int64,key))
        end
    end
    len_index_candidate::Int64 =length(index_candidate)=#

    #=index_candidate::Vector{Int64} = []
    index::Int64 = 1
    for (key,value) in neighbor_set
        if value >= length(g[key]) - value
            push!(index_candidate,parse(Int64,key))
        end
    end
    len_index_candidate::Int64 =length(index_candidate)=#

    index_candidate2::Vector{Int64} = Vector{Int64}(undef,size(y_cur)[1])
    len_index_candidate2 = length(index_candidate2)
    for i in 1:len_index_candidate2
        index_candidate2[i] = i
    end

    Obj_new::Float64 = Obj_cur
    demand_new::Int64 = demand_cur
    index_fliped_set::Array{Int64, 2} = Array{Int64,2}(undef,num_flip,2)
    index_fliped::Int64 = 0
    neighbor_list::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef,num_flip)
    len::Vector{Int64} = Vector{Int64}(undef,num_flip)
    sum_len::Int64 = 0
    i::Int64 = 1
    k::Int64 = 0

    for i in 1:num_flip
        if len_index_candidate == 0
            while index_fliped == 0
                k = rand(2:len_index_candidate2)
                index_fliped = index_candidate2[k]
            end
            index_candidate2[k] = 0
        elseif len_index_candidate < num_flip
            if i <= len_index_candidate
                index_fliped = index_candidate[i]
            else
                index_fliped = 1
            end
        else
            while index_fliped == 0
                k = rand(1:len_index_candidate)
                index_fliped = index_candidate[k]
            end
            index_candidate[k] = 0
        end
        
        index_fliped_set[i,1] = index_fliped
        if num_flip == 1
            index_fliped_set[i,2] = 1
        elseif index_fliped == 1
            index_fliped_set[i,2] = 0
        else
            index_fliped_set[i,2] = rand(0:1)
        end
        if index_fliped_set[i,2] == 1
            if y_cur_copy[index_fliped,1] == 1
                y_cur_copy[index_fliped,1] = 0 
                demand_new = demand_new - y_cur_copy[index_fliped,2]
            else
                y_cur_copy[index_fliped,1] = 1
                demand_new = demand_new + y_cur_copy[index_fliped,2]
            end
            neighbor_list[i] = g["$index_fliped"]
            len[i] = length(neighbor_list[i])
            sum_len = sum_len + len[i]
        else
           neighbor_list[i] = []
           len[i] = 0
        end
        index_fliped = 0
    end

    Obj_new = Obj_cur + 2*ceil(demand_cur/C) - 2*ceil(demand_new/C)
    edges_changed::Array{Int64,2} = Array{Int64,2}(undef,sum_len,2)
    index = 1
    change::Int64 = 0
    for i in 1:num_flip
        if index_fliped_set[i,2] == 1
            for j in 1:len[i]
                if w_cur_copy[index_fliped_set[i,1],neighbor_list[i][j]] != 0
                    if w_cur_copy[index_fliped_set[i,1],neighbor_list[i][j]] == 1
                        w_cur_copy[index_fliped_set[i,1],neighbor_list[i][j]] = -1
                        Obj_new = Obj_new - xij[index_fliped_set[i,1],neighbor_list[i][j]]
                        edges_changed[index,:] = [index_fliped_set[i,1],neighbor_list[i][j]]
                    else
                        w_cur_copy[index_fliped_set[i,1],neighbor_list[i][j]] = 1
                        Obj_new = Obj_new + xij[index_fliped_set[i,1],neighbor_list[i][j]]
                        edges_changed[index,:] = [index_fliped_set[i,1],neighbor_list[i][j]]
                    end
                else
                    if w_cur_copy[neighbor_list[i][j],index_fliped_set[i,1]] == 1
                        w_cur_copy[neighbor_list[i][j],index_fliped_set[i,1]] = -1
                        Obj_new = Obj_new - xij[neighbor_list[i][j],index_fliped_set[i,1]]
                        edges_changed[index,:] = [neighbor_list[i][j],index_fliped_set[i,1]]
                    else
                        w_cur_copy[neighbor_list[i][j],index_fliped_set[i,1]] = 1
                        Obj_new = Obj_new + xij[neighbor_list[i][j],index_fliped_set[i,1]]
                        edges_changed[index,:] = [neighbor_list[i][j],index_fliped_set[i,1]]
                    end
                end
                index = index + 1
            end
            change = 1
        end
    end
    ΔObj_value::Float64 = Obj_new - Obj_cur

    if change == 1
        if (ΔObj_value <= 0 || rand() < exp(-ΔObj_value/T_cur)) && mod == 1
            for i in 1:num_flip
                y_cur[index_fliped_set[i,1],1] = y_cur_copy[index_fliped_set[i,1],1]
            end
            for i in 1:sum_len
                if w_cur[edges_changed[i,1],edges_changed[i,2]] == w_cur_copy[edges_changed[i,1],edges_changed[i,2]]  
                    continue
                end
                w_cur[edges_changed[i,1],edges_changed[i,2]] = w_cur_copy[edges_changed[i,1],edges_changed[i,2]]
                if w_cur[edges_changed[i,1],edges_changed[i,2]] == 1
                    if edges_changed[i,1] != 1
                        if !haskey(neighbor_set, "$(edges_changed[i,1])")
                            neighbor_set["$(edges_changed[i,1])"] = 1
                        else
                            neighbor_set["$(edges_changed[i,1])"] = neighbor_set["$(edges_changed[i,1])"] + 1
                        end
                    end
                    if edges_changed[i,2] != 1
                        if !haskey(neighbor_set, "$(edges_changed[i,2])")
                            neighbor_set["$(edges_changed[i,2])"] =  1
                        else
                            neighbor_set["$(edges_changed[i,2])"] = neighbor_set["$(edges_changed[i,2])"] + 1
                        end
                    end
                else
                    if edges_changed[i,1] != 1
                        if neighbor_set["$(edges_changed[i,1])"] == 1
                            delete!(neighbor_set, "$(edges_changed[i,1])")
                        else
                            neighbor_set["$(edges_changed[i,1])"] = neighbor_set["$(edges_changed[i,1])"] - 1
                        end
                    end
                    if edges_changed[i,2] != 1
                        if neighbor_set["$(edges_changed[i,2])"] == 1
                            delete!(neighbor_set, "$(edges_changed[i,2])")
                        else
                            neighbor_set["$(edges_changed[i,2])"] = neighbor_set["$(edges_changed[i,2])"] - 1
                        end
                    end
                end
            end
            demand_cur = demand_new
            Obj_cur = Obj_new
        else
            for i in 1:num_flip
                y_cur_copy[index_fliped_set[i,1],1] = y_cur[index_fliped_set[i,1],1]
            end
            for i in 1:sum_len
                w_cur_copy[edges_changed[i,1],edges_changed[i,2]] = w_cur[edges_changed[i,1],edges_changed[i,2]]
            end
            demand_new = demand_cur
            Obj_new = Obj_cur
        end
    end

    if mod == 1
        return Vector{Any}([Obj_cur,demand_cur])
    else
        return ΔObj_value
    end

end

function neighbor_sol_n_filp_in_delta_S_and_demand(y_cur::Array{Int64,2}, y_cur_copy::Array{Int64,2}, w_cur::SparseMatrixCSC{Int64, Int64}, w_cur_copy::SparseMatrixCSC{Int64, Int64}, g::Dict{String, Vector{Int64}}, demand_cur::Int64, Obj_cur::Float64, xij::SparseMatrixCSC{Float64, Int64}, C::Int64, T_cur::Float64, mod::Int64, neighbor_set::Dict{String, Int64},num_flip::Int64)
    index_candidate_in_S::Vector{Int64} = []
    index_candidate_not_in_S::Vector{Int64} = []
    index_candidate1::Vector{Int64} = Vector{Int64}(undef,size(y_cur)[1])
    len_index_candidate1 = length(index_candidate1)
    index::Int64 = 0
    num_remove::Int64 = rand(0:num_flip)
    num_add::Int64 = num_flip - num_remove
    Obj_new::Float64 = Obj_cur
    demand_new::Int64 = demand_cur
    k::Int64 = 0

    index_fliped_set::Array{Int64, 2} = Array{Int64,2}(undef,num_flip,2)
    index_fliped::Int64 = 0
    neighbor_list::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef,num_flip)
    len::Vector{Int64} = Vector{Int64}(undef,num_flip)
    sum_len::Int64 = 0
    demand_decreased::Int64 = 0
    for i in 1:len_index_candidate1 
        index_candidate1[i] = i
    end

    if num_remove != 0
        for (key,value) in neighbor_set
            index = parse(Int64,key)
            if y_cur[index,1] == 1
                if ceil((demand_cur - num_remove*y_cur[index,2])/C) == ceil((demand_cur)/C)
                    push!(index_candidate_in_S,index)
                end
            end
        end
        len_index_candidate_in_S::Int64 =length(index_candidate_in_S)

        index_candidate2::Vector{Int64} = [] ##그냥 S안에 있는 게아니라 S안에 있는것중 RHS 감소 안 시키는 거에 대해 뽑기로 설정해볼 수 있음
        for i in 2:size(y_cur)[1]
            if y_cur[i,1] == 1
                push!(index_candidate2,i)
            end
        end
        len_index_candidate2 = length(index_candidate2)

        for i in 1:num_remove
            if len_index_candidate_in_S == 0
                if len_index_candidate2 >= num_remove
                    while index_fliped == 0
                        k = rand(1:len_index_candidate2)
                        index_fliped = index_candidate2[k]
                    end
                    index_candidate2[k] = 0
                else
                    while index_fliped == 0
                        k = rand(2:len_index_candidate1)
                        index_fliped = index_candidate1[k]
                    end
                    index_candidate1[k] = 0
                end
            elseif len_index_candidate_in_S < num_remove
                if i <= len_index_candidate_in_S
                    index_fliped = index_candidate_in_S[i]
                else
                    index_fliped = 1
                end
            else
                while index_fliped == 0
                    k = rand(1:len_index_candidate_in_S)
                    index_fliped = index_candidate_in_S[k]
                end
                index_candidate_in_S[k] = 0
            end
            
            index_fliped_set[i,1] = index_fliped
            
            if index_fliped == 1
                index_fliped_set[i,2] = 0
            else
                index_fliped_set[i,2] = rand(0:1) ##무조건 빼기로 설정해볼 수 있음
            end

            if index_fliped_set[i,2] == 1
                if y_cur_copy[index_fliped,1] == 1
                    y_cur_copy[index_fliped,1] = 0 
                    demand_new = demand_new - y_cur_copy[index_fliped,2]
                    demand_decreased = demand_decreased + y_cur_copy[index_fliped,2]
                else
                    y_cur_copy[index_fliped,1] = 1
                    demand_new = demand_new + y_cur_copy[index_fliped,2]
                end
                neighbor_list[i] = g["$index_fliped"]
                len[i] = length(neighbor_list[i])
                sum_len = sum_len + len[i]
            else
            neighbor_list[i] = []
            len[i] = 0
            end
            index_fliped = 0
        end
    end

    index_fliped = 0
    if num_add != 0
        for (key,value) in neighbor_set
            index = parse(Int64,key)
            if y_cur[index,1] == 0
                if ceil(((demand_cur - demand_decreased) + num_add*y_cur[index,2])/C) > ceil((demand_cur - demand_decreased)/C)
                    push!(index_candidate_not_in_S,index)
                end
            end
        end
        len_index_candidate_not_in_S::Int64 =length(index_candidate_not_in_S)

        index_candidate3::Vector{Int64} = [] 
        for i in 2:size(y_cur)[1]
            if y_cur[i,1] == 0
                push!(index_candidate3,i)
            end
        end
        len_index_candidate3 = length(index_candidate3)

        for i in 1 + num_remove:num_add + num_remove
            if len_index_candidate_not_in_S < num_add
                if len_index_candidate3 >= num_add
                    while index_fliped == 0
                        k = rand(1:len_index_candidate3)
                        index_fliped = index_candidate3[k]
                    end
                    index_candidate3[k] = 0
                    index_fliped_set[i,1] = index_fliped
                    index_fliped_set[i,2] = rand(0:1)
                else
                    while index_fliped == 0
                        k = rand(2:len_index_candidate1)
                        index_fliped = index_candidate1[k]
                    end
                    index_candidate1[k] = 0
                    index_fliped_set[i,1] = index_fliped
                    index_fliped_set[i,2] = rand(0:1)
                end
            else
                while index_fliped == 0
                    k = rand(1:len_index_candidate_not_in_S)
                    index_fliped = index_candidate_not_in_S[k]
                end
                index_candidate_not_in_S[k] = 0
                index_fliped_set[i,1] = index_fliped
                index_fliped_set[i,2] = 1
            end

            if index_fliped_set[i,2] == 1
                if y_cur_copy[index_fliped,1] == 1
                    y_cur_copy[index_fliped,1] = 0 
                    demand_new = demand_new - y_cur_copy[index_fliped,2]
                else
                    y_cur_copy[index_fliped,1] = 1
                    demand_new = demand_new + y_cur_copy[index_fliped,2]
                end
                neighbor_list[i] = g["$index_fliped"]
                len[i] = length(neighbor_list[i])
                sum_len = sum_len + len[i]
            else
            neighbor_list[i] = []
            len[i] = 0
            end
            index_fliped = 0
        end
    end

    Obj_new = Obj_cur + 2*ceil(demand_cur/C) - 2*ceil(demand_new/C)

    edges_changed::Array{Int64,2} = Array{Int64,2}(undef,sum_len,2)
    index = 1
    change::Int64 = 0
    for i in 1:num_flip
        if index_fliped_set[i,2] == 1
            for j in 1:len[i]
                if w_cur_copy[index_fliped_set[i,1],neighbor_list[i][j]] != 0
                    if w_cur_copy[index_fliped_set[i,1],neighbor_list[i][j]] == 1
                        w_cur_copy[index_fliped_set[i,1],neighbor_list[i][j]] = -1
                        Obj_new = Obj_new - xij[index_fliped_set[i,1],neighbor_list[i][j]]
                        edges_changed[index,:] = [index_fliped_set[i,1],neighbor_list[i][j]]
                    else
                        w_cur_copy[index_fliped_set[i,1],neighbor_list[i][j]] = 1
                        Obj_new = Obj_new + xij[index_fliped_set[i,1],neighbor_list[i][j]]
                        edges_changed[index,:] = [index_fliped_set[i,1],neighbor_list[i][j]]
                    end
                else
                    if w_cur_copy[neighbor_list[i][j],index_fliped_set[i,1]] == 1
                        w_cur_copy[neighbor_list[i][j],index_fliped_set[i,1]] = -1
                        Obj_new = Obj_new - xij[neighbor_list[i][j],index_fliped_set[i,1]]
                        edges_changed[index,:] = [neighbor_list[i][j],index_fliped_set[i,1]]
                    else
                        w_cur_copy[neighbor_list[i][j],index_fliped_set[i,1]] = 1
                        Obj_new = Obj_new + xij[neighbor_list[i][j],index_fliped_set[i,1]]
                        edges_changed[index,:] = [neighbor_list[i][j],index_fliped_set[i,1]]
                    end
                end
                index = index + 1
            end
            change = 1
        end
    end
    ΔObj_value::Float64 = Obj_new - Obj_cur

    if change == 1
        if (ΔObj_value <= 0 || rand() < exp(-ΔObj_value/T_cur)) && mod == 1
            for i in 1:num_flip
                y_cur[index_fliped_set[i,1],1] = y_cur_copy[index_fliped_set[i,1],1]
            end
            for i in 1:sum_len
                if w_cur[edges_changed[i,1],edges_changed[i,2]] == w_cur_copy[edges_changed[i,1],edges_changed[i,2]]  
                    continue
                end
                w_cur[edges_changed[i,1],edges_changed[i,2]] = w_cur_copy[edges_changed[i,1],edges_changed[i,2]]
                if w_cur[edges_changed[i,1],edges_changed[i,2]] == 1
                    if edges_changed[i,1] != 1
                        if !haskey(neighbor_set, "$(edges_changed[i,1])")
                            neighbor_set["$(edges_changed[i,1])"] = 1
                        else
                            neighbor_set["$(edges_changed[i,1])"] = neighbor_set["$(edges_changed[i,1])"] + 1
                        end
                    end
                    if edges_changed[i,2] != 1
                        if !haskey(neighbor_set, "$(edges_changed[i,2])")
                            neighbor_set["$(edges_changed[i,2])"] =  1
                        else
                            neighbor_set["$(edges_changed[i,2])"] = neighbor_set["$(edges_changed[i,2])"] + 1
                        end
                    end
                else
                    if edges_changed[i,1] != 1
                        if neighbor_set["$(edges_changed[i,1])"] == 1
                            delete!(neighbor_set, "$(edges_changed[i,1])")
                        else
                            neighbor_set["$(edges_changed[i,1])"] = neighbor_set["$(edges_changed[i,1])"] - 1
                        end
                    end
                    if edges_changed[i,2] != 1
                        if neighbor_set["$(edges_changed[i,2])"] == 1
                            delete!(neighbor_set, "$(edges_changed[i,2])")
                        else
                            neighbor_set["$(edges_changed[i,2])"] = neighbor_set["$(edges_changed[i,2])"] - 1
                        end
                    end
                end
            end
            demand_cur = demand_new
            Obj_cur = Obj_new
        else
            for i in 1:num_flip
                y_cur_copy[index_fliped_set[i,1],1] = y_cur[index_fliped_set[i,1],1]
            end
            for i in 1:sum_len
                w_cur_copy[edges_changed[i,1],edges_changed[i,2]] = w_cur[edges_changed[i,1],edges_changed[i,2]]
            end
            demand_new = demand_cur
            Obj_new = Obj_cur
        end
    end

    if mod == 1
        return Vector{Any}([Obj_cur,demand_cur])
    else
        return ΔObj_value
    end
end

function initial_temperature(y_cur::Array{Int64,2}, y_cur_copy::Array{Int64,2}, w_cur::SparseMatrixCSC{Int64, Int64}, w_cur_copy::SparseMatrixCSC{Int64, Int64}, g::Dict{String, Vector{Int64}}, demand_cur::Int64, Obj_cur::Float64, xij::SparseMatrixCSC{Float64, Int64}, C::Int64, neighbor_set::Dict{String, Int64})
    Δf::Float64 = 0
    m1::Int64 = 0
    m2::Int64 = 0
    T0::Float64 = 0.0
    ΔObj_value::Float64 = 0.0
    num_hamming_ball::Int64 = 2
    num_flip::Int64 = rand(2:3)
    ver::Int64 = rand(1:3)
    for i in 1:1000
        ver = rand(1:3)
        num_flip = rand(2:3)

        ver == 3
        num_flip = 3
        if ver == 1
            ΔObj_value = neighbor_sol_hamming_ball(y_cur,y_cur_copy,w_cur,w_cur_copy,g,demand_cur,Obj_cur,xij,C,T0,0,num_hamming_ball,neighbor_set,1)
        elseif ver == 2
            ΔObj_value= neighbor_sol_n_filp(y_cur,y_cur_copy,w_cur,w_cur_copy,g,demand_cur,Obj_cur,xij,C,T0,0,num_flip,neighbor_set,1)
        else
            ΔObj_value= neighbor_sol_n_filp_in_delta_S(y_cur,y_cur_copy,w_cur,w_cur_copy,g,demand_cur,Obj_cur,xij,C,T0,0,neighbor_set,num_flip)
        end
        #ΔObj_value = neighbor_sol_n_filp_in_delta_S_and_demand(y_cur,y_cur_copy,w_cur,w_cur_copy,g,demand_cur,Obj_cur,xij,C,T0,0,neighbor_set,num_flip)
        if ΔObj_value > 0
            Δf = Δf + ΔObj_value
            m2 = m2 + 1
        else
            m1 = m1 + 1
        end
    end

    Δf_avg::Float64 = Δf/m2
    
    if (m2*accept_rate - m1*(1-accept_rate)) <= 0
        T0 = 100.0
    else
        T0 = Δf_avg/log(m2/(m2*accept_rate - m1*(1-accept_rate)))
    end

    if isnan(T0)
        T0 = 100.0
    end

    return T0
end

function simulated_annealing(y0::Array{Int64,2}, w0::SparseMatrixCSC{Int64, Int64}, g::Dict{String, Vector{Int64}}, demand0::Int64, Obj0::Float64, xij::SparseMatrixCSC{Float64, Int64}, C::Int64, T0::Float64, neighbor_set::Dict{String, Int64})
    #deepcopy 대신 내부 리스트 하나하나 copy하는 것으로 수정하면 더 빨라질 듯
    y_opt::Array{Int64,2} = deepcopy(y0)
    w_opt::SparseMatrixCSC{Int64, Int64} = deepcopy(w0)
    Obj_opt::Float64 = Obj0
    y_cur::Array{Int64,2} = y0
    y_cur_copy::Array{Int64,2} = deepcopy(y_cur) 
    w_cur::SparseMatrixCSC{Int64, Int64} = w0
    w_cur_copy::SparseMatrixCSC{Int64, Int64} = deepcopy(w_cur)
    Obj_cur::Float64 = Obj0
    demand_cur::Int64 = demand0
    T::Float64 = T0
    Obj_cur2::Float64 = Obj0
    demand_cur2::Int64 = demand0
    result::Int64 = 0
    num_hamming_ball::Int64 = 2
    num_flip::Int64 = 2
    ver::Int64 = rand(1:3)
    improve::Int64 = rand(0:1)
    Obj_cur3::Float64 = Obj0
    demand_cur3::Int64 = demand0
    y_cur2::Array{Int64,2} = y0
    y_cur_copy2::Array{Int64,2} = y0
    w_cur2::SparseMatrixCSC{Int64, Int64} = w0
    w_cur_copy2::SparseMatrixCSC{Int64, Int64} = w0
    neighbor_set2::Dict{String, Int64} = deepcopy(neighbor_set)
    println(T0)
    while(T >= 0.001)
        for i in 1:L#ceil(size(y_cur)[1]/num_hamming_ball)^num_hamming_ball #2*length(neighbor_set)*(length(neighbor_set)-1)
            num_flip = rand(2:3)
            ver = rand(1:3)
            improve = rand(0:1)
            
            num_flip = 3
            ver = 3
            if ver == 1
                Obj_cur, demand_cur= neighbor_sol_hamming_ball(y_cur,y_cur_copy,w_cur,w_cur_copy,g,demand_cur,Obj_cur,xij,C,T,1,num_hamming_ball,neighbor_set,1)
            elseif ver == 2
                Obj_cur, demand_cur= neighbor_sol_n_filp(y_cur,y_cur_copy,w_cur,w_cur_copy,g,demand_cur,Obj_cur,xij,C,T,1,num_flip,neighbor_set,1)
            else
                Obj_cur, demand_cur= neighbor_sol_n_filp_in_delta_S(y_cur,y_cur_copy,w_cur,w_cur_copy,g,demand_cur,Obj_cur,xij,C,T,1,neighbor_set,num_flip)
            end
            #Obj_cur, demand_cur = neighbor_sol_n_filp_in_delta_S_and_demand(y_cur,y_cur_copy,w_cur,w_cur_copy,g,demand_cur,Obj_cur,xij,C,T,1,neighbor_set,num_flip)

            #=
            y_cur2 = deepcopy(y_cur)
            y_cur_copy2 = deepcopy(y_cur)
            w_cur2 = deepcopy(w_cur)
            w_cur_copy2 = deepcopy(w_cur)
            neighbor_set2 = deepcopy(neighbor_set)
            Obj_cur3, demand_cur3 = neighbor_sol_n_filp_in_delta_S(y_cur2,y_cur_copy2,w_cur2,w_cur_copy2,g,demand_cur,Obj_cur,xij,C,T,1,neighbor_set2,num_flip)
            if Obj_cur3 < Obj_cur
                Obj_cur = Obj_cur3
                demand_cur = demand_cur3
                y_cur = deepcopy(y_cur2)
                y_cur_copy = deepcopy(y_cur_copy2)
                w_cur = deepcopy(w_cur2)
                w_cur_copy = deepcopy(w_cur_copy2)
                neighbor_set = deepcopy(neighbor_set2)
            end
            =#

            Obj_cur2, demand_cur2,result = neighbor_sol_demand(y_cur,y_cur_copy,w_cur,w_cur_copy,g,demand_cur,Obj_cur,xij,C,T,neighbor_set,1)
            if result != 0
                Obj_cur = Obj_cur2
                demand_cur = demand_cur2
            end

            improve = 0
            if improve == 1
                y_cur2 = deepcopy(y_cur)
                y_cur_copy2 = deepcopy(y_cur)
                w_cur2 = deepcopy(w_cur)
                w_cur_copy2 = deepcopy(w_cur)
                neighbor_set2 = deepcopy(neighbor_set)

                Obj_cur2, demand_cur2,result = neighbor_sol_demand(y_cur2,y_cur_copy2,w_cur2,w_cur_copy2,g,demand_cur,Obj_cur,xij,C,T,neighbor_set2,1)
                if result != 0
                    if Obj_cur2 < Obj_cur
                        Obj_cur = Obj_cur2
                        demand_cur = demand_cur2
                        y_cur = deepcopy(y_cur2)
                        y_cur_copy = deepcopy(y_cur_copy2)
                        w_cur = deepcopy(w_cur2)
                        w_cur_copy = deepcopy(w_cur_copy2)
                        neighbor_set = deepcopy(neighbor_set2)
                    end
                end
            end

            if Obj_cur < Obj_opt
                y_opt = deepcopy(y_cur)
                w_opt = deepcopy(w_cur)  
                Obj_opt = Obj_cur
            end
        end
        T = α*T
    end

    return (y_opt, Obj_opt)

end


α::Float64  = 0.98
accept_rate::Float64  = 0.9
#data = CSV.read("test.csv",DataFrame)
data = CSV.read("rci_separation_cvrp200.csv",DataFrame)
L::Int64  = 2*data[1,"customers"]*(data[1,"customers"]-1)
#L = data[1,"customers"]

#=#test
input = read_Kth_data(data, 2)
g = support_graph(input)
sol_initial = initial_sol_CVRPSep(input[4],input[5],input[7], input[3],input[1],input[6])
if sol_initial == 0 
    sol_initial = initial_sol_max_subset(input[1],input[4],input[5], input[6], input[3])
end
println(sol_initial[1])
sol_initial_copy_y = deepcopy(sol_initial[1])
sol_initial_copy_w = deepcopy(sol_initial[2])
T0 = initial_temperature(sol_initial[1],sol_initial_copy_y,sol_initial[2],sol_initial_copy_w ,g,sol_initial[3], sol_initial[4], input[6],input[3],sol_initial[5])
sol_opt = simulated_annealing(sol_initial[1], sol_initial[2], g, sol_initial[3], sol_initial[4], input[6],input[3], T0, sol_initial[5])
println(sol_opt[2], sol_opt[1])=#

@time begin
    len = size(data)[1]
    sol_obj::Vector{Float64} = Vector{Float64}(undef,len)
    sol_set::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef,len)
    for i in 1:len
        sol_set[i] = []
    end
    for i in 1:len
        input = read_Kth_data(data,i)  
        g = support_graph(input)
        #sol_initial = initial_sol_CVRPSep(input[4],input[5],input[7], input[3],input[1],input[6])
        sol_initial = initial_sol_CVRPSep_best(input[4],input[5],input[7], input[3],input[1],input[6])
        if sol_initial == 0 
            sol_initial = initial_sol_random(input[1],input[4],input[5], input[6], input[3])
        end
        #sol_initial = initial_sol_max_subset(input[1],input[4],input[5], input[6], input[3])
        #sol_initial = initial_sol_random(input[1],input[4],input[5], input[6], input[3])
        #sol_initial = initial_sol_one_subset(input[1],input[4],input[5], input[6], input[3], g)
        sol_initial_copy_y = deepcopy(sol_initial[1])
        sol_initial_copy_w = deepcopy(sol_initial[2])
        T0 = initial_temperature(sol_initial[1],sol_initial_copy_y,sol_initial[2],sol_initial_copy_w ,g,sol_initial[3], sol_initial[4], input[6],input[3],sol_initial[5])
        sol_opt = simulated_annealing(sol_initial[1], sol_initial[2], g, sol_initial[3], sol_initial[4], input[6],input[3], T0, sol_initial[5])
        sol_obj[i] = -sol_opt[2]
        for j in 1:size(sol_opt[1])[1]
            if sol_opt[1][j][1] == 1
                push!(sol_set[i],j)
            end
        end
        println(i,": ",sol_obj[i],sol_set[i])
    end
end

data_comparative = CSV.read("rci_exact_violation_max_200.csv",DataFrame)
violation_exact::Vector{Float64} = data_comparative[!,"violation_optimal"]

performance::Vector{Float64} = Vector{Float64}(undef,length(violation_exact))
success::Vector{Int64} = Vector{Int64}(undef,length(violation_exact))


for i in 1:length(performance)
    if sol_obj[i] >= 0.0001
        success[i] = 1
    else
        success[i] = 0
    end
    #performance[i] = sol_obj[i]
    if violation_exact[i] !=0
        performance[i] = min(100,round((violation_exact[i]-sol_obj[i])/violation_exact[i]*100.0,digits = 2))
    else
        performance[i] = 0
    end
end
df = DataFrame(
    optimal_gap = getindex.(performance, 1), 
    success = getindex.(success, 1),
    objective_value = getindex.(sol_obj)
)
CSV.write("output_ver2_2/performance_200_ver40.csv", df)