using Plots
using DataFrames
using JuMP


function plot_all_columns_df(df, title, filename) 
    index_values = 1:nrow(df)
    a = plot(title=title)
    for col in names(df)
        plot!(index_values, df[:, col], label=col, legend=:outertopright)
    end
    display(a)
    savefig(a, "C:/Users/vegardvk/PycharmProjects/Bernstein/results/$filename")
end

function dense_array_to_df(weights)
    cols = length(weights[1, :])
    array = collect(weights)
    column_names = ["$c" for c in 1:cols]
    df = DataFrame(array, column_names)
    return df
end

function print_model_info(model::Model)
    # Get variables and constraints

    variables = all_variables(model)
    num_binary = 0
    num_integer = 0
    num_continuous = 0
    for var in variables
        if is_binary(var)
            num_binary += 1
        elseif is_integer(var)
            num_integer += 1
        else
            num_continuous += 1
        end
    end
    println("Number of binary variables: ", num_binary)
    println("Number of integer variables: ", num_integer)
    println("Number of continuous variables: ", num_continuous)
    # num_binary = num_constraints(model, VariableRef, MOI.ZeroOne)

    constraints = all_constraints(model, include_variable_in_set_constraints=true)
    println("Number of constraints: ", length(constraints))
    print("\n\n")
end



function get_load(areas, timesteps)
    file_path = "consumption.xlsx"
    xf = XLSX.readxlsx(file_path)
    sh = xf["Sheet1"]
    demand = sh[2:end, 3]
    multiplier = 200
    demand = [multiplier * demand[i] for i in 1:(length(demand))]

    df = DataFrame()
    for a in 1:areas
        demand_array = demand[(1 + timesteps * (a-1)):(timesteps+timesteps*(a-1))]
        df[!, "$a"] = demand_array
    end
    return df
end


function add_power_plant_data!(data_dict)
    df = CSV.read("gen.csv", DataFrame)

    gen_ub_array = df[:, "PMax MW"]
    gen_ub_dict = Dict(i => v for (i, v) in enumerate(gen_ub_array))
    gen_lb_array = df[:, "PMin MW"]
    gen_lb_dict = Dict(i => v for (i, v) in enumerate(gen_lb_array))
    fuel_price_array = df[:, "Fuel Price \$/MMBTU"]
    fuel_price_dict = Dict(i => v for (i, v) in enumerate(fuel_price_array))

    n_thermal_plants = length(gen_ub_array)
    data_dict["n_thermal_plants"] =  n_thermal_plants
    # data_dict["P_t"] = 1:n_thermal_plants
    data_dict["gen_ub"] = gen_ub_dict
    data_dict["gen_lb"] = gen_lb_dict
    data_dict["fuel_price"] = fuel_price_dict

    p_max = df[:, "PMax MW"]
    p_min = df[:, "PMin MW"]
    fuel = df[:, "Fuel"]
    min_down_time = df[:, "Min Down Time Hr"]
    min_up_time = df[:, "Min Up Time Hr"]
    fuel_price = df[:, "Fuel Price \$/MMBTU"]
    a = 3
    return data_dict
end

function get_module_data!(data_dict)
    file_path = "C:/Users/vegardvk/vscodeProjects/bernstein/Input/moduldata.xlsx"
    xf = XLSX.readxlsx(file_path)
    sh = xf["Sheet1"]
    data = Dict()
    data["modnr"] = sh[2:end, 1][:, 1]
    data["modnavn"] = sh[2:end, 2][:, 1]
    data["vassdrag"] = sh[2:end, 3][:, 1]
    data["kap_mag"] = sh[2:end, 5][:, 1] .* 1000
    data["kap_forb"] = sh[2:end, 8][:, 1] .* 1000
    data["kap_gen"] =  sh[2:end, 9][:, 1] .* 1000
    data["enekv"] = sh[2:end, 10][:, 1]
    data["topo_gen"] = sh[2:end, 15][:, 1]
    data["topo_forb"] = sh[2:end, 16][:, 1]
    data["topo_flom"] = sh[2:end, 17][:, 1]
    data["tilsig_reg"] = sh[2:end, 27][:, 1]
    data["tilsig_reg_serie"] = sh[2:end, 28][:, 1]
    data["tilsig_ureg"] = sh[2:end, 29][:, 1]
    data["tilsig_ureg_serie"] = sh[2:end, 30][:, 1]
    
    full_df = DataFrame(data)

    df = full_df[full_df.vassdrag .== "NIDELVA_H.DETD", :]
    
    n_hydro_plants = length(df[!, "modnr"])

    # n_thermal_plants = data_dict["n_thermal_plants"]
    # P_h = (n_thermal_plants+1):(n_thermal_plants+n_hydro_plants)
    # P = vcat(data_dict["P_t"], P_h)

    df[!, "P_h"] = df[!, "modnr"]
    data_dict["P"] = vcat(1:20, df[!, "P_h"])
    data_dict["P_h"] = df[!, "modnr"]

    key_list = ["enekv", "kap_spill", "kap_forb", "kap_mag", "kap_gen", "tilsig_reg", "tilsig_ureg", "topo_flom", "topo_forb", "topo_gen", "starting_reservoir", "end_wv"]
    for key in key_list
        data_dict[key] = Dict()
    end
    for i in 1:nrow(df)
        data_dict["gen_ub"][df[i, "P_h"]] = df[i, "kap_gen"]
        data_dict["gen_lb"][df[i, "P_h"]] = 0
        data_dict["enekv"][df[i, "P_h"]] = df[i, "enekv"]
        data_dict["kap_gen"][df[i, "P_h"]] = df[i, "kap_gen"]
        data_dict["kap_forb"][df[i, "P_h"]] = df[i, "kap_forb"]
        data_dict["kap_mag"][df[i, "P_h"]] = df[i, "kap_mag"]
        data_dict["kap_spill"][df[i, "P_h"]] = 100000
        data_dict["tilsig_reg"][df[i, "P_h"]] = df[i, "tilsig_reg_serie"]
        data_dict["tilsig_ureg"][df[i, "P_h"]] = df[i, "tilsig_ureg_serie"]
        data_dict["topo_flom"][df[i, "P_h"]] = df[i, "topo_flom"]
        data_dict["topo_forb"][df[i, "P_h"]] = df[i, "topo_forb"]
        data_dict["topo_gen"][df[i, "P_h"]] = df[i, "topo_gen"]
    end
    find_connected_plants!(data_dict, df)
    add_dummy_hydro_parameters!(data_dict, df)
end

function find_connected_plants!(data_dict, df)
    I_disch = Dict(key => [] for key in data_dict["P_h"])
    I_bypass = Dict(key => [] for key in data_dict["P_h"])
    I_spill = Dict(key => [] for key in data_dict["P_h"])

    for p in data_dict["P_h"]
        row_index = findfirst(df.P_h .== p)
        push!(get!(I_disch, df.topo_gen[row_index], []), p)
        push!(get!(I_spill, df.topo_flom[row_index], []), p)
        push!(get!(I_bypass, df.topo_forb[row_index], []), p)
    end
    data_dict["I_disch"] = I_disch
    data_dict["I_spill"] = I_spill
    data_dict["I_bypass"] = I_bypass
end

function add_dummy_hydro_parameters!(data_dict, df)
    for p in data_dict["P_h"]
        row_index = findfirst(df.P_h .== p)
        if p == 49927
            data_dict["starting_reservoir"][p] = df[row_index, "kap_mag"]*0.99
        else
            data_dict["starting_reservoir"][p] = df[row_index, "kap_mag"]*0.99
        end
        data_dict["fuel_price"][p] = 0
    end

    power_plant_list = data_dict["P_h"]
    end_reservoir = 49900
    energy_price = 3
    set_wv!(data_dict["I_disch"], data_dict["fuel_price"], end_reservoir, data_dict["enekv"], energy_price)

end


function get_inflow(timesteps)
    data = Dict()
    add_power_plant_data!(data)
    get_module_data!(data)
    plants = data["P_h"]
    df = DataFrame()
    for p in plants
        if data["kap_mag"][p] > 0
            if p == 49927
                df[!, "$p"] = [500 for t in 1:timesteps]
            else
                df[!, "$p"] = [500 for t in 1:timesteps]
            end
        else
            df[!, "$p"] = [0 for t in 1:timesteps]
        end
    end
    # display(df)
    return df
end

function set_wv!(I_disch, fuel_price, p, enekv, energy_price)
    for p2 in I_disch[p]
        fuel_price[p2] = energy_price * enekv[p2] + fuel_price[p]
        # println("Fuel price in $p2 = $energy_price * $(enekv[p2]) + $(fuel_price[p]) = $(fuel_price[p2])")
        if length(I_disch[p2]) > 0
            set_wv!(I_disch, fuel_price, p2, enekv, energy_price)
        end
    end
end