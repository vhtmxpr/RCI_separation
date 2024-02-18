using CSV, DataFrames, JSON
α::Float64  = 0.95
accept_rate::Float64  = 0.8
data = CSV.read("rci_separation_cvrp50.csv",DataFrame)
L::Int64  = 2*data[1,"customers"]*(data[1,"customers"]-1)

function read_Kth_data(data, k)

    customers::Int64 = data[k,"customers"]
    
    K::Int64 = data[k,"vehicles"]
    
    C::Int64 = data[k,"capacity"]

    edges_raw::Vector{Vector{Int64}} = Vector{Vector{Int64}}(JSON.parse(data[k,"edges"]))
   
    xij_raw::Vector{Float64} = Vector{Float64}(JSON.parse(data[k,"x_bar"]))
    
    demands::Vector{Int64} = Vector{Int64}(JSON.parse(data[k,"demand"]))
    
    edges::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef,Int64(length(edges_raw)/2))
    
    xij::Dict{String, Float64} = Dict{String, Float64}()
    
    index::Int64 = 1
    
    for i in 1:2:length(edges_raw)  # get undirected support graph s.t Xij >0
        if xij_raw[i] > 0 # Xij = 0인 edge는 가져올 필요가 없음
            edges[index] = edges_raw[i] 
            xij["($(edges_raw[i][1]),$(edges_raw[i][2]))"] =  xij_raw[i]
            index = index + 1
        end
    end

    input = [customers, K, C, demands, edges, xij]
    
    return input
end

function support_graph(input)
    
    g::Dict{String, Vector{Int64}} = Dict{String, Vector{Int64}}()
    
    for i in 1:input[1]+1
        g["$i"] = []
    end
    
    for i in 1:length(input[5])
        push!(g["$(input[5][i][1])"],input[5][i][2])
        push!(g["$(input[5][i][2])"],input[5][i][1])
    end

    return g
end

function initial_sol_random(customer::Int64, demand::Vector{Int64}, edges::Vector{Vector{Int64}}, xij::Dict{String,Float64}, C::Int64)
#function initial_sol_random(input::Vector{Any}) //일단은 위의 방식이 약간 더 빨라 보이나 거의 비슷
    #customer::Int64 = input[1]
    #demand::Vector{Int64} = input[4]
    #edges::Vector{Vector{Int64}} = input[5]
    #xij::Dict{String,Float64} = input[6]

    sol_initial_y::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef,customer + 1)
    sum_demand::Int64 = 0

    sol_initial_y[1] = [0,0]
    num::Int64 = 0

    for i in 2:customer + 1 
        num = rand(0:1)
        sol_initial_y[i] = [num, demand[i]]
        if num == 1
            sum_demand = sum_demand + demand[i]
        end
    end

    sol_initial_w::Dict{String, Int64} = Dict{String, Int64}()
    value::Int64 = 0
    Obj_value::Float64 = 0.0

    for i in 1:length(edges)
        if sol_initial_y[edges[i][1]][1] == 1 && sol_initial_y[edges[i][2]][1] == 0
            value = 1
        elseif sol_initial_y[edges[i][1]][1] == 0 && sol_initial_y[edges[i][2]][1] == 1
            value = 1
        else
            value = 0
        end
        sol_initial_w["($(edges[i][1]),$(edges[i][2]))"] = value
        if value == 1
            Obj_value = Obj_value + xij["($(edges[i][1]),$(edges[i][2]))"]*value
        end
    end

    Obj_value = Obj_value - 2*ceil(sum_demand/C)

    return [sol_initial_y,sol_initial_w, sum_demand, Obj_value]
end

function initial_sol_max_subset(customer::Int64, demand::Vector{Int64}, edges::Vector{Vector{Int64}}, xij::Dict{String,Float64}, C::Int64)
    sol_initial_y::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef,customer + 1)
    sum_demand::Int64 = sum(demand)
    
    sol_initial_y[1] = [0,0]
    for i in 2:customer + 1 
        sol_initial_y[i] = [1, demand[i]]
    end

    sol_initial_w::Dict{String, Int64} = Dict{String, Int64}()
    value::Int64 = 0
    Obj_value::Float64 = 0.0

    for i in 1:length(edges)
        if sol_initial_y[edges[i][1]][1] == 1 && sol_initial_y[edges[i][2]][1] == 0
            value = 1
        elseif sol_initial_y[edges[i][1]][1] == 0 && sol_initial_y[edges[i][2]][1] == 1
            value = 1
        else
            value = 0
        end
        sol_initial_w["($(edges[i][1]),$(edges[i][2]))"] = value
        if value == 1
            Obj_value = Obj_value + xij["($(edges[i][1]),$(edges[i][2]))"]*value
        end
    end

    Obj_value = Obj_value - 2*ceil(sum_demand/C)

    return [sol_initial_y,sol_initial_w, sum_demand, Obj_value]
    
end

function initial_sol_one_subset(customer::Int64, demand::Vector{Int64}, edges::Vector{Vector{Int64}}, xij::Dict{String,Float64}, C::Int64, g::Dict{String, Vector{Int64}})
    sol_initial_y::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef,customer + 1)

    sol_initial_y[1] = [0,0]
    for i in 2:customer + 1 
        sol_initial_y[i] = [0, demand[i]]
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
    sol_initial_y[index_max_neighbor][1] = 1
    sum_demand::Int64 = demand[index_max_neighbor]
    
    sol_initial_w::Dict{String, Int64} = Dict{String, Int64}()
    neighbor_set::Dict{String, Vector{Int64}} =Dict{String, Vector{Int64}}()

    value::Int64 = 0
    Obj_value::Float64 = 0.0

    for i in 1:length(edges)
        if sol_initial_y[edges[i][1]][1] == 1 && sol_initial_y[edges[i][2]][1] == 0
            value = 1
        elseif sol_initial_y[edges[i][1]][1] == 0 && sol_initial_y[edges[i][2]][1] == 1
            value = 1
        else
            value = 0
        end
        sol_initial_w["($(edges[i][1]),$(edges[i][2]))"] = value
        if value == 1
            if edges[i][1] != 1
                if !haskey(neighbor_set, "$(edges[i][1])")
                    neighbor_set["$(edges[i][1])"] = [edges[i][1],1]
                else
                    neighbor_set["$(edges[i][1])"] = [neighbor_set["$(edges[i][1])"][1],neighbor_set["$(edges[i][1])"][2] + 1]
                end
            end
            if edges[i][2] != 1
                if !haskey(neighbor_set, "$(edges[i][2])") 
                    neighbor_set["$(edges[i][2])"] = [edges[i][2],1]
                else
                    neighbor_set["$(edges[i][2])"] = [neighbor_set["$(edges[i][2])"][1],neighbor_set["$(edges[i][2])"][2] + 1]
                end
            end 
            Obj_value = Obj_value + xij["($(edges[i][1]),$(edges[i][2]))"]*value
        end
    end

    Obj_value = Obj_value - 2*ceil(sum_demand/C)

    return [sol_initial_y,sol_initial_w, sum_demand, Obj_value,neighbor_set]
    
end

function stop(T::Float64)

    stop::Int64 = 0
    
    if(T < 0.01)
        stop = 1
    end

    return stop
end

function neighbor_sol(y_cur::Vector{Vector{Int64}}, y_cur_copy::Vector{Vector{Int64}}, w_cur::Dict{String, Int64}, w_cur_copy::Dict{String, Int64}, g::Dict{String, Vector{Int64}}, demand_cur::Int64, Obj_cur::Float64, xij::Dict{String, Float64}, C::Int64, T_cur::Float64, mod::Int64)

    index_fliped::Int64 = rand(2:length(y_cur_copy))
    Obj_new::Float64 = Obj_cur
    demand_new::Int64 = demand_cur
    
    if y_cur_copy[index_fliped][1] == 1
        y_cur_copy[index_fliped][1] = 0
        demand_new = demand_new - y_cur_copy[index_fliped][2]
    else
        y_cur_copy[index_fliped][1] = 1
        demand_new = demand_new + y_cur_copy[index_fliped][2]
    end
    
    Obj_new = Obj_cur + 2*ceil(demand_cur/C) - 2*ceil(demand_new/C)
    
    neighbor_list::Vector{Int64} = g["$index_fliped"]
    len::Int64 = length(neighbor_list)
    edges_changed::Vector{Vector{Int64}}  = Vector{Vector{Int64}}(undef,len)
    
    for i in 1:len
        if haskey(w_cur_copy,"($(index_fliped),$(neighbor_list[i]))") == true
            if w_cur_copy["($(index_fliped),$(neighbor_list[i]))"] == 1
                w_cur_copy["($(index_fliped),$(neighbor_list[i]))"] = 0
                Obj_new = Obj_new - xij["($(index_fliped),$(neighbor_list[i]))"]
                edges_changed[i] = [index_fliped,neighbor_list[i]]
            else
                w_cur_copy["($(index_fliped),$(neighbor_list[i]))"] = 1
                Obj_new = Obj_new + xij["($(index_fliped),$(neighbor_list[i]))"]
                edges_changed[i] = [index_fliped,neighbor_list[i]]
            end
        else
            if w_cur_copy["($(neighbor_list[i]),$(index_fliped))"] == 1
                w_cur_copy["($(neighbor_list[i]),$(index_fliped))"] = 0
                Obj_new = Obj_new - xij["($(neighbor_list[i]),$(index_fliped))"]
                edges_changed[i] = [neighbor_list[i],index_fliped]
                
            else
                w_cur_copy["($(neighbor_list[i]),$(index_fliped))"] = 1
                Obj_new = Obj_new + xij["($(neighbor_list[i]),$(index_fliped))"]
                edges_changed[i] = [neighbor_list[i],index_fliped]
            end
        end
    end
    
    ΔObj_value::Float64 = Obj_new - Obj_cur

    if mod == 1
        if ΔObj_value <= 0 || rand() < exp(-ΔObj_value/T_cur)
            y_cur[index_fliped][1] = y_cur_copy[index_fliped][1]
            for i in 1:len
                w_cur["($(edges_changed[i][1]),$(edges_changed[i][2]))"] = w_cur_copy["($(edges_changed[i][1]),$(edges_changed[i][2]))"]
            end
            demand_cur = demand_new
            Obj_cur = Obj_new
        else
            y_cur_copy[index_fliped][1] = y_cur[index_fliped][1]
            for i in 1:len
                w_cur_copy["($(edges_changed[i][1]),$(edges_changed[i][2]))"] = w_cur["($(edges_changed[i][1]),$(edges_changed[i][2]))"]
            end
            demand_new = demand_cur
            Obj_new = Obj_cur
        end
    else
        y_cur_copy[index_fliped][1] = y_cur[index_fliped][1]
        for i in 1:len
            w_cur_copy["($(edges_changed[i][1]),$(edges_changed[i][2]))"] = w_cur["($(edges_changed[i][1]),$(edges_changed[i][2]))"]
        end
        demand_new = demand_cur
        Obj_new = Obj_cur
    end
    
    
    if mod == 1
        return [Obj_cur, demand_cur]
    else
        return ΔObj_value
    end
end

function neighbor_sol_2flip(y_cur::Vector{Vector{Int64}}, y_cur_copy::Vector{Vector{Int64}}, w_cur::Dict{String, Int64}, w_cur_copy::Dict{String, Int64}, g::Dict{String, Vector{Int64}}, demand_cur::Int64, Obj_cur::Float64, xij::Dict{String, Float64}, C::Int64, T_cur::Float64, mod::Int64)

    index_fliped1::Int64 = rand(2:length(y_cur_copy))
    index_fliped2::Int64 = index_fliped1
    while index_fliped1 == index_fliped2
        index_fliped2 = rand(2:length(y_cur_copy))
    end

    Obj_new::Float64 = Obj_cur
    demand_new::Int64 = demand_cur
    rand1::Int64 = rand(0:1)
    rand2::Int64 = rand(0:1)

    if rand1 == 1
        if y_cur_copy[index_fliped1][1] == 1
            y_cur_copy[index_fliped1][1] = 0
            demand_new = demand_new - y_cur_copy[index_fliped1][2]
        else
            y_cur_copy[index_fliped1][1] = 1
            demand_new = demand_new + y_cur_copy[index_fliped1][2]
        end
    end
    if rand2 == 1
        if y_cur_copy[index_fliped2][1] == 1
            y_cur_copy[index_fliped2][1] = 0
            demand_new = demand_new - y_cur_copy[index_fliped2][2]
        else
            y_cur_copy[index_fliped2][1] = 1
            demand_new = demand_new + y_cur_copy[index_fliped2][2]
        end
    end
    Obj_new = Obj_cur + 2*ceil(demand_cur/C) - 2*ceil(demand_new/C)
    
    neighbor_list::Vector{Int64} = []
    neighbor_list2::Vector{Int64} = []

    if rand1 == 1
        neighbor_list = g["$index_fliped1"]
    end

    if rand2 == 1
        neighbor_list2 = g["$index_fliped2"]
    end

    len1::Int64 = length(neighbor_list)
    len2::Int64 = length(neighbor_list2)
    len::Int64 = len1 + len2
    edges_changed::Vector{Vector{Int64}}  = Vector{Vector{Int64}}(undef,len)
    
    if rand1 == 1
        for i in 1:len1
            if haskey(w_cur_copy,"($(index_fliped1),$(neighbor_list[i]))") == true
                if w_cur_copy["($(index_fliped1),$(neighbor_list[i]))"] == 1
                    w_cur_copy["($(index_fliped1),$(neighbor_list[i]))"] = 0
                    Obj_new = Obj_new - xij["($(index_fliped1),$(neighbor_list[i]))"]
                    edges_changed[i] = [index_fliped1,neighbor_list[i]]
                else
                    w_cur_copy["($(index_fliped1),$(neighbor_list[i]))"] = 1
                    Obj_new = Obj_new + xij["($(index_fliped1),$(neighbor_list[i]))"]
                    edges_changed[i] = [index_fliped1,neighbor_list[i]]
                end
            else
                if w_cur_copy["($(neighbor_list[i]),$(index_fliped1))"] == 1
                    w_cur_copy["($(neighbor_list[i]),$(index_fliped1))"] = 0
                    Obj_new = Obj_new - xij["($(neighbor_list[i]),$(index_fliped1))"]
                    edges_changed[i] = [neighbor_list[i],index_fliped1]
                    
                else
                    w_cur_copy["($(neighbor_list[i]),$(index_fliped1))"] = 1
                    Obj_new = Obj_new + xij["($(neighbor_list[i]),$(index_fliped1))"]
                    edges_changed[i] = [neighbor_list[i],index_fliped1]
                end
            end
        end
    end

    if rand2 == 1
        for i in 1:len2
            if haskey(w_cur_copy,"($(index_fliped2),$(neighbor_list2[i]))") == true
                if w_cur_copy["($(index_fliped2),$(neighbor_list2[i]))"] == 1
                    w_cur_copy["($(index_fliped2),$(neighbor_list2[i]))"] = 0
                    Obj_new = Obj_new - xij["($(index_fliped2),$(neighbor_list2[i]))"]
                    edges_changed[i+len1] = [index_fliped2,neighbor_list2[i]]
                else
                    w_cur_copy["($(index_fliped2),$(neighbor_list2[i]))"] = 1
                    Obj_new = Obj_new + xij["($(index_fliped2),$(neighbor_list2[i]))"]
                    edges_changed[i+len1] = [index_fliped2,neighbor_list2[i]]
                end
            else
                if w_cur_copy["($(neighbor_list2[i]),$(index_fliped2))"] == 1
                    w_cur_copy["($(neighbor_list2[i]),$(index_fliped2))"] = 0
                    Obj_new = Obj_new - xij["($(neighbor_list2[i]),$(index_fliped2))"]
                    edges_changed[i+len1] = [neighbor_list2[i],index_fliped2]
                    
                else
                    w_cur_copy["($(neighbor_list2[i]),$(index_fliped2))"] = 1
                    Obj_new = Obj_new + xij["($(neighbor_list2[i]),$(index_fliped2))"]
                    edges_changed[i+len1] = [neighbor_list2[i],index_fliped2]
                end
            end
        end
    end
    
    ΔObj_value::Float64 = Obj_new - Obj_cur

    if rand1 == 1||rand2 == 1
        if mod == 1
            if ΔObj_value <= 0 || rand() < exp(-ΔObj_value/T_cur)
                y_cur[index_fliped1][1] = y_cur_copy[index_fliped1][1]
                y_cur[index_fliped2][1] = y_cur_copy[index_fliped2][1]
                for i in 1:len
                    w_cur["($(edges_changed[i][1]),$(edges_changed[i][2]))"] = w_cur_copy["($(edges_changed[i][1]),$(edges_changed[i][2]))"]
                end
                demand_cur = demand_new
                Obj_cur = Obj_new
            else
                y_cur_copy[index_fliped1][1] = y_cur[index_fliped1][1]
                y_cur_copy[index_fliped2][1] = y_cur[index_fliped2][1]
                for i in 1:len
                    w_cur_copy["($(edges_changed[i][1]),$(edges_changed[i][2]))"] = w_cur["($(edges_changed[i][1]),$(edges_changed[i][2]))"]
                end
                demand_new = demand_cur
                Obj_new = Obj_cur
            end
        else
            y_cur_copy[index_fliped1][1] = y_cur[index_fliped1][1]
            y_cur_copy[index_fliped2][1] = y_cur[index_fliped2][1]
            for i in 1:len
                w_cur_copy["($(edges_changed[i][1]),$(edges_changed[i][2]))"] = w_cur["($(edges_changed[i][1]),$(edges_changed[i][2]))"]
            end
            demand_new = demand_cur
            Obj_new = Obj_cur
        end
    end

    if mod == 1
        return [Obj_cur, demand_cur]
    else
        return ΔObj_value
    end
end

function neighbor_sol_2flip_in_delat_S(y_cur::Vector{Vector{Int64}}, y_cur_copy::Vector{Vector{Int64}}, w_cur::Dict{String, Int64}, w_cur_copy::Dict{String, Int64}, g::Dict{String, Vector{Int64}}, demand_cur::Int64, Obj_cur::Float64, xij::Dict{String, Float64}, C::Int64, T_cur::Float64, mod::Int64, neighbor_set::Dict{String,Vector{Int64}})    

    Obj_new::Float64 = Obj_cur
    demand_new::Int64 = demand_cur
    index_filped_set::Vector{Int64} = Vector{Int64}(undef,length(neighbor_set))
    index::Int64 = 1
    for (key,value) in neighbor_set
        index_filped_set[index] = value[1]
        index = index + 1
    end

    index_fliped1::Int64 = 0
    index_fliped2::Int64 = 0
    rand1::Int64 = 1
    rand2::Int64 = 1
    
    while true
        if length(index_filped_set) == 0
            index_fliped1 = rand(2:length(y_cur))
            index_fliped2 = index_fliped1
            while index_fliped1 == index_fliped2
                index_fliped2 = rand(2:length(y_cur))
            end
            rand1 = rand(0:1)
            rand2 = rand(0:1)
        elseif length(index_filped_set) == 1
            index_fliped1 = index_filped_set[1]
            index_fliped2 = index_fliped1
            rand1 = 1
            rand2 = 0
        else
            index_fliped1 = index_filped_set[rand(1:length(index_filped_set))]
            index_fliped2 = index_fliped1
            while index_fliped1 == index_fliped2
                index_fliped2 = index_filped_set[rand(1:length(index_filped_set))]
            end
            rand1 = rand(0:1)
            rand2 = rand(0:1)
        end

        if rand1 == 1
            if y_cur_copy[index_fliped1][1] == 1
                y_cur_copy[index_fliped1][1] = 0
                demand_new = demand_new - y_cur_copy[index_fliped1][2]
            else
                y_cur_copy[index_fliped1][1] = 1
                demand_new = demand_new + y_cur_copy[index_fliped1][2]
            end
        end
        if rand2 == 1
            if y_cur_copy[index_fliped2][1] == 1
                y_cur_copy[index_fliped2][1] = 0
                demand_new = demand_new - y_cur_copy[index_fliped2][2]
            else
                y_cur_copy[index_fliped2][1] = 1
                demand_new = demand_new + y_cur_copy[index_fliped2][2]
            end
        end
        if demand_new == 0
            if rand1 == 1
                if y_cur_copy[index_fliped1][1] == 1
                    y_cur_copy[index_fliped1][1] = 0
                    demand_new = demand_new - y_cur_copy[index_fliped1][2]
                else
                    y_cur_copy[index_fliped1][1] = 1
                    demand_new = demand_new + y_cur_copy[index_fliped1][2]
                end
            end
            if rand2 == 1
                if y_cur_copy[index_fliped2][1] == 1
                    y_cur_copy[index_fliped2][1] = 0
                    demand_new = demand_new - y_cur_copy[index_fliped2][2]
                else
                    y_cur_copy[index_fliped2][1] = 1
                    demand_new = demand_new + y_cur_copy[index_fliped2][2]
                end
            end
        else
            break
        end
    end
    Obj_new = Obj_cur + 2*ceil(demand_cur/C) - 2*ceil(demand_new/C)

    neighbor_list::Vector{Int64} = []
    neighbor_list2::Vector{Int64} = []

    if rand1 == 1
        neighbor_list = g["$index_fliped1"]
    end

    if rand2 == 1
        neighbor_list2 = g["$index_fliped2"]
    end

    len1::Int64 = length(neighbor_list)
    len2::Int64 = length(neighbor_list2)
    len::Int64 = len1 + len2
    edges_changed::Vector{Vector{Int64}}  = Vector{Vector{Int64}}(undef,len)
    
    if rand1 == 1
        for i in 1:len1
            if haskey(w_cur_copy,"($(index_fliped1),$(neighbor_list[i]))") == true
                if w_cur_copy["($(index_fliped1),$(neighbor_list[i]))"] == 1
                    w_cur_copy["($(index_fliped1),$(neighbor_list[i]))"] = 0
                    Obj_new = Obj_new - xij["($(index_fliped1),$(neighbor_list[i]))"]
                    edges_changed[i] = [index_fliped1,neighbor_list[i]]

                else
                    w_cur_copy["($(index_fliped1),$(neighbor_list[i]))"] = 1
                    Obj_new = Obj_new + xij["($(index_fliped1),$(neighbor_list[i]))"]
                    edges_changed[i] = [index_fliped1,neighbor_list[i]]
                end
            else
                if w_cur_copy["($(neighbor_list[i]),$(index_fliped1))"] == 1
                    w_cur_copy["($(neighbor_list[i]),$(index_fliped1))"] = 0
                    Obj_new = Obj_new - xij["($(neighbor_list[i]),$(index_fliped1))"]
                    edges_changed[i] = [neighbor_list[i],index_fliped1]
                    
                else
                    w_cur_copy["($(neighbor_list[i]),$(index_fliped1))"] = 1
                    Obj_new = Obj_new + xij["($(neighbor_list[i]),$(index_fliped1))"]
                    edges_changed[i] = [neighbor_list[i],index_fliped1]
                end
            end
        end
    end

    if rand2 == 1
        for i in 1:len2
            if haskey(w_cur_copy,"($(index_fliped2),$(neighbor_list2[i]))") == true
                if w_cur_copy["($(index_fliped2),$(neighbor_list2[i]))"] == 1
                    w_cur_copy["($(index_fliped2),$(neighbor_list2[i]))"] = 0
                    Obj_new = Obj_new - xij["($(index_fliped2),$(neighbor_list2[i]))"]
                    edges_changed[i+len1] = [index_fliped2,neighbor_list2[i]]
                else
                    w_cur_copy["($(index_fliped2),$(neighbor_list2[i]))"] = 1
                    Obj_new = Obj_new + xij["($(index_fliped2),$(neighbor_list2[i]))"]
                    edges_changed[i+len1] = [index_fliped2,neighbor_list2[i]]
                end
            else
                if w_cur_copy["($(neighbor_list2[i]),$(index_fliped2))"] == 1
                    w_cur_copy["($(neighbor_list2[i]),$(index_fliped2))"] = 0
                    Obj_new = Obj_new - xij["($(neighbor_list2[i]),$(index_fliped2))"]
                    edges_changed[i+len1] = [neighbor_list2[i],index_fliped2]
                    
                else
                    w_cur_copy["($(neighbor_list2[i]),$(index_fliped2))"] = 1
                    Obj_new = Obj_new + xij["($(neighbor_list2[i]),$(index_fliped2))"]
                    edges_changed[i+len1] = [neighbor_list2[i],index_fliped2]
                end
            end
        end
    end
    ΔObj_value::Float64 = Obj_new - Obj_cur

    if rand1 == 1||rand2 == 1
        if mod == 1
            if ΔObj_value <= 0 || rand() < exp(-ΔObj_value/T_cur)
                y_cur[index_fliped1][1] = y_cur_copy[index_fliped1][1]
                y_cur[index_fliped2][1] = y_cur_copy[index_fliped2][1]
                for i in 1:len
                    #2 flip에서 하나는 set에 있던 노드 i가 나가고, 하나는 없던 노드 j가 들어오면 wij와 wij_copy 값은 같음
                    #이 경우에 edges_changed에는 (i,j)가 두번 저장되게 되는데 이 경우 neighbor_set[i][2] 값을 2번 변화시키는데, 변화시키면 안됨(wij값이 변하지 않으므로)
                    if w_cur["($(edges_changed[i][1]),$(edges_changed[i][2]))"] == w_cur_copy["($(edges_changed[i][1]),$(edges_changed[i][2]))"]  
                        continue
                    end
                    w_cur["($(edges_changed[i][1]),$(edges_changed[i][2]))"] = w_cur_copy["($(edges_changed[i][1]),$(edges_changed[i][2]))"]
                    if w_cur["($(edges_changed[i][1]),$(edges_changed[i][2]))"] == 1
                        if edges_changed[i][1] != 1
                            if !haskey(neighbor_set, "$(edges_changed[i][1])")
                                neighbor_set["$(edges_changed[i][1])"] = [edges_changed[i][1],1]
                            else
                                neighbor_set["$(edges_changed[i][1])"][2] = neighbor_set["$(edges_changed[i][1])"][2]+1
                            end
                        end
                        if edges_changed[i][2] != 1
                            if !haskey(neighbor_set, "$(edges_changed[i][2])")
                                neighbor_set["$(edges_changed[i][2])"] = [edges_changed[i][2],1]
                            else
                                neighbor_set["$(edges_changed[i][2])"][2] = neighbor_set["$(edges_changed[i][2])"][2]+1
                            end
                        end
                    else
                        if edges_changed[i][1] != 1
                            if neighbor_set["$(edges_changed[i][1])"][2] == 1
                                delete!(neighbor_set, "$(edges_changed[i][1])")
                            else
                                neighbor_set["$(edges_changed[i][1])"][2] = neighbor_set["$(edges_changed[i][1])"][2]-1
                            end
                        end
                        if edges_changed[i][2] != 1
                            if neighbor_set["$(edges_changed[i][2])"][2] == 1
                                delete!(neighbor_set, "$(edges_changed[i][2])")
                            else
                                neighbor_set["$(edges_changed[i][2])"][2] = neighbor_set["$(edges_changed[i][2])"][2]-1
                            end
                        end
                    end
                end
                demand_cur = demand_new
                Obj_cur = Obj_new
                #println("a")
            else
                y_cur_copy[index_fliped1][1] = y_cur[index_fliped1][1]
                y_cur_copy[index_fliped2][1] = y_cur[index_fliped2][1]
                for i in 1:len
                    w_cur_copy["($(edges_changed[i][1]),$(edges_changed[i][2]))"] = w_cur["($(edges_changed[i][1]),$(edges_changed[i][2]))"]
                end
                demand_new = demand_cur
                Obj_new = Obj_cur
                #println("b")
            end
        else
            y_cur_copy[index_fliped1][1] = y_cur[index_fliped1][1]
            y_cur_copy[index_fliped2][1] = y_cur[index_fliped2][1]
            for i in 1:len
                w_cur_copy["($(edges_changed[i][1]),$(edges_changed[i][2]))"] = w_cur["($(edges_changed[i][1]),$(edges_changed[i][2]))"]
            end
            demand_new = demand_cur
            Obj_new = Obj_cur
        end
    end
    
    if mod == 1
        return [Obj_cur, demand_cur]
    else
        return ΔObj_value
    end
end

function initial_temperature(y_cur::Vector{Vector{Int64}}, y_cur_copy::Vector{Vector{Int64}}, w_cur::Dict{String, Int64}, w_cur_copy::Dict{String, Int64}, g::Dict{String, Vector{Int64}}, demand_cur::Int64, Obj_cur::Float64, xij::Dict{String, Float64}, C::Int64,neighbor_set::Dict{String,Vector{Int64}})

    Δf::Float64 = 0
    m1::Int64 = 0
    m2::Int64 = 0
    T0::Float64 = 0.0
    ΔObj_value::Float64 = 0.0

    for i in 1:1000
        ΔObj_value = neighbor_sol_2flip_in_delat_S(y_cur,y_cur_copy,w_cur,w_cur_copy,g,demand_cur,Obj_cur,xij,C,T0,0,neighbor_set)

        if ΔObj_value > 0
            Δf = Δf + ΔObj_value
            m2 = m2 + 1
        else
            m1 = m1 + 1
        end
    end

    Δf_avg::Float64 = Δf/m2
    T0 = Δf_avg/log(m2/(m2*accept_rate - m1*(1-accept_rate)))

    if isnan(T0)
        T0 = 10.0
    end

    return T0
end

function simulated_annealing(y0::Vector{Vector{Int64}}, w0::Dict{String, Int64}, g::Dict{String, Vector{Int64}}, demand0::Int64, Obj0::Float64, xij::Dict{String, Float64}, C::Int64, T0::Float64,neighbor_set::Dict{String,Vector{Int64}})

    #deepcopy 대신 내부 리스트 하나하나 copy하는 것으로 수정하면 더 빨라질 듯
    y_opt::Vector{Vector{Int64}} = deepcopy(y0)
    w_opt::Dict{String, Int64} = deepcopy(w0)
    Obj_opt::Float64 = Obj0
    y_cur::Vector{Vector{Int64}} = y0
    y_cur_copy::Vector{Vector{Int64}} = deepcopy(y_cur) 
    w_cur::Dict{String, Int64} = w0
    w_cur_copy::Dict{String, Int64} = deepcopy(w_cur)
    Obj_cur::Float64 = Obj0
    demand_cur::Int64 = demand0
    T::Float64 = T0


    while(T >= 0.01)
        for i in 1:2*length(neighbor_set)*(length(neighbor_set)-1)
            #println(neighbor_set)
            Obj_cur, demand_cur= neighbor_sol_2flip_in_delat_S(y_cur,y_cur_copy,w_cur,w_cur_copy,g,demand_cur,Obj_cur,xij,C,T,1,neighbor_set)
            if Obj_cur < Obj_opt
                y_opt = deepcopy(y_cur)
                w_opt = deepcopy(w_cur)  
                Obj_opt = Obj_cur
            end
        end
        T = α*T
    end

    return (y_opt, Obj_opt)
    #return (y_cur, Obj_cur)


end

#test
#input = read_Kth_data(data, 51)  
#g = support_graph(input)
#sol_initial = initial_sol_one_subset(input[1],input[4],input[5], input[6], input[3], g)
#sol_initial_copy_y = deepcopy(sol_initial[1])
#sol_initial_copy_w = deepcopy(sol_initial[2])
#T0 = initial_temperature(sol_initial[1],sol_initial_copy_y,sol_initial[2],sol_initial_copy_w ,g,sol_initial[3], sol_initial[4], input[6],input[3],sol_initial[5])
#sol_opt = simulated_annealing(sol_initial[1], sol_initial[2], g, sol_initial[3], sol_initial[4], input[6],input[3], T0,sol_initial[5])
#println(sol_opt[2], sol_opt[1])

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
        #sol_initial = initial_sol_max_subset(input[1],input[4],input[5], input[6], input[3])
        #sol_initial = initial_sol_random(input[1],input[4],input[5], input[6], input[3])
        sol_initial = initial_sol_one_subset(input[1],input[4],input[5], input[6], input[3], g)
        sol_initial_copy_y = deepcopy(sol_initial[1])
        sol_initial_copy_w = deepcopy(sol_initial[2])
        T0 = initial_temperature(sol_initial[1],sol_initial_copy_y,sol_initial[2],sol_initial_copy_w ,g,sol_initial[3], sol_initial[4], input[6],input[3],sol_initial[5])
        sol_opt = simulated_annealing(sol_initial[1], sol_initial[2], g, sol_initial[3], sol_initial[4], input[6],input[3], T0,sol_initial[5])
        sol_obj[i] = -sol_opt[2]
        for j in 1:length(sol_opt[1])
            if sol_opt[1][j][1] == 1
                push!(sol_set[i],j)
            end
        end
    end
end

data_comparative = CSV.read("rci_exact_violation_max_50.csv",DataFrame)
violation_exact::Vector{Float64} = data_comparative[!,"violation_optimal"]

performance::Vector{Float64} = Vector{Float64}(undef,length(violation_exact))
success::Vector{Int64} = Vector{Int64}(undef,length(violation_exact))

for i in 1:length(performance)
    if sol_obj[i] > 0
        success[i] = 1
    else
        success[i] = 0
    end
    performance[i] = round((violation_exact[i]-sol_obj[i])/violation_exact[i]*100.0,digits = 2)
end

df = DataFrame(
    optimal_gap = getindex.(performance, 1), 
    success = getindex.(success, 1)
)
CSV.write("performance_50_ver6.csv", df)