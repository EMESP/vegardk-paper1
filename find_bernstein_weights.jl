using Combinatorics
using JuMP
using Plots
using Ipopt
using DataFrames
using CSV
using AxisArrays
using XLSX

include("C:/Users/vegardvk/vscodeProjects/bernstein/helper_functions.jl")

function define_problem(parameters, cont_constraints=true)
    B = parameters["B"]
    T = parameters["T"]
    S = parameters["S"]
    bernstein_degree = parameters["bernstein_degree"]
    
    model = Model()
    set_optimizer(model, Ipopt.Optimizer)

    @variable(model, weights[b in B, t in T] ≥ 0)
    @variable(model, deviation[s in S, t in T])

    @constraint(model, deviation_summation[s in S, t in T], deviation[s, t] == sum(parameters["bernstein_curves"][b+1, s+1] * weights[b, t] for b in B) - parameters["target_values"][s+1, t])
    if cont_constraints
        @constraint(model, continuity_constraint1[t in T[1:end-1]], weights[bernstein_degree, t] == weights[0, t+1])
        @constraint(model, continuity_constraint2[t in T[1:end-1]], weights[bernstein_degree, t] - weights[bernstein_degree-1, t] == weights[1, t+1] - weights[0, t+1])
    end
    
    @objective(model, Min, sum(deviation[s, t]^2 for s in S for t in T))
    optimize!(model)

    weights = value.(weights)
    df = dense_array_to_df(weights)
    
    # df = get_converted_df_separated_timesteps(weights, S[end])
    converted_list = get_converted_list(df, S[end])
    # display(plot(converted_list))

    # plot_all_columns_df(df)

    # print(model)
    return df

end


function define_parameters(target_values, bernstein_degree, sampling_points)
    time_steps = length(target_values)
    T = 1:time_steps
    B = 0:bernstein_degree
    S = 0:sampling_points
    bernstein_curves = zeros(bernstein_degree+1, sampling_points+1)
    target_values_expanded = zeros(sampling_points+1, time_steps)
    for s in S
        for b in B
            bernstein_curves[b+1, s+1] = get_bernstein_val(bernstein_degree, b, s/sampling_points)
        end
        for t in T
            target_values_expanded[s+1, t] = target_values[t]
        end
    end

    data = Dict(
        "bernstein_degree" => bernstein_degree,
        "bernstein_curves" => bernstein_curves,
        "target_values" => target_values_expanded,
        "B" => B,
        "S" => S,
        "T" => T
    )
    return data
end


function get_bernstein_val(B, b, s)
    n_choose_i = binomial(B, b)
    val = n_choose_i * (s^b)*((1-s)^(B-b))
    return val
end


function get_converted_list(weights_df, sampling_points)
    bernstein_degree = length(weights_df[:,1])-1
    time_steps = length(weights_df[1, :])
    converted_list = zeros(time_steps * sampling_points+1)

    for b in 0:bernstein_degree
        converted_list[1] += get_bernstein_val(bernstein_degree, b, 0) * weights_df[b+1, 1]        
        for t in 1:time_steps
            for s in 1:sampling_points
                # println((t-1)*sampling_points+s+1)
                converted_list[(t-1)*sampling_points+s+1] += get_bernstein_val(bernstein_degree, b, s/sampling_points) * weights_df[b+1, t]
            end
        end
    end
    a = 3
    return converted_list
end


function get_converted_df_separated_timesteps(weights, sampling_points)
    bernstein_degree = length(weights[:,1])-1
    time_steps = length(weights[0, :])
    df = DataFrame()
    for t in 1:time_steps
        # array = zeros(time_steps * sampling_points+1)
        array = fill(NaN, time_steps * sampling_points+1)
        for s in 0:sampling_points
            array[(t-1) * sampling_points+s+1] = 0
            for b in 0:bernstein_degree
                array[(t-1) * sampling_points+s+1] += get_bernstein_val(bernstein_degree, b, s/sampling_points) * weights[b, t]
            end
        end
        df[!, "$t"] = array
        # empty_df.t = array
        # hcat!(empty_df, new_df)
    end
    return df
end


function find_and_write_demand_weights(bernstein_degree)
    load_df = DataFrame(XLSX.readtable("output/load_data.xlsx", "Sheet1", infer_eltypes=true))
    weights_df = DataFrame()
    A = unique(load_df.area)
    for a in A
        area_load = load_df[load_df.area .== a, :Forbruk]
        parameters = define_parameters(area_load, bernstein_degree, 100)
        df = define_problem(parameters)
        df.b = 0:bernstein_degree
        df_long = stack(df, Not(:b), variable_name="timestep", value_name="Forbruk")
        df_long.timestep = parse.(Int, df_long.timestep)
        df_long.area .= a
        weights_df = vcat(weights_df, df_long)
    end
    weights_df.Forbruk = round.(weights_df.Forbruk, digits=2)
    XLSX.writetable("output/load_weights.xlsx", weights_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end

function find_and_write_shedding_weights(bernstein_degree, cont_constraints)
    ts_df = DataFrame(XLSX.readtable("discrete_results/results.xlsx", "area_results", infer_eltypes=true))
    weights_df = DataFrame()
    A = unique(ts_df.area)
    for a in A
        area_ts_df = ts_df[ts_df.area .== a, :load_shedding]
        parameters = define_parameters(area_ts_df, bernstein_degree, 100)
        df = define_problem(parameters, cont_constraints)
        df.b = 0:bernstein_degree
        df_long = stack(df, Not(:b), variable_name="timestep", value_name="load_shedding")
        df_long.timestep = parse.(Int, df_long.timestep)
        df_long.area .= a
        weights_df = vcat(weights_df, df_long)
    end
    weights_df.load_shedding = round.(weights_df.load_shedding, digits=2)
    XLSX.writetable("output/shedding_weights.xlsx", weights_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end

function find_and_write_dumping_weights(bernstein_degree, cont_constraints)
    ts_df = DataFrame(XLSX.readtable("discrete_results/results.xlsx", "area_results", infer_eltypes=true))
    weights_df = DataFrame()
    A = unique(ts_df.area)
    for a in A
        area_ts_df = ts_df[ts_df.area .== a, :power_dumping]
        parameters = define_parameters(area_ts_df, bernstein_degree, 100)
        df = define_problem(parameters, cont_constraints)
        df.b = 0:bernstein_degree
        df_long = stack(df, Not(:b), variable_name="timestep", value_name="power_dumping")
        df_long.timestep = parse.(Int, df_long.timestep)
        df_long.area .= a
        weights_df = vcat(weights_df, df_long)
    end
    weights_df.power_dumping = round.(weights_df.power_dumping, digits=2)
    XLSX.writetable("output/dumping_weights.xlsx", weights_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end

function find_and_write_production_weights(bernstein_degree, cont_constraints)
    prod_df = DataFrame(XLSX.readtable("discrete_results/results.xlsx", "production", infer_eltypes=true))
    weights_df = DataFrame()
    P = unique(prod_df.plant_id)
    for p in P
        powerplant_prod = prod_df[prod_df.plant_id .== p, :production]
        parameters = define_parameters(powerplant_prod, bernstein_degree, 100)
        df = define_problem(parameters, cont_constraints)
        df.b = 0:bernstein_degree
        df_long = stack(df, Not(:b), variable_name="timestep", value_name="production")
        df_long.timestep=parse.(Int, df_long.timestep)
        df_long.plant_id .= p
        weights_df = vcat(weights_df, df_long)
    end
    weights_df.production = round.(weights_df.production, digits=2)
    XLSX.writetable("output/production_weights.xlsx", weights_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end

function find_and_write_wind_weights(bernstein_degree)
    wind_df = DataFrame(XLSX.readtable("output/wind_ts_data.xlsx", "Sheet1", infer_eltypes=true))
    weights_df = DataFrame()
    P = unique(wind_df.plant_id)
    for p in P
        plant_wind_ts = wind_df[wind_df.plant_id .== p, :wind_power]
        parameters = define_parameters(plant_wind_ts, bernstein_degree, 100)
        df = define_problem(parameters)
        df.b = 0:bernstein_degree
        df_long = stack(df, Not(:b), variable_name="timestep", value_name="wind_power")
        df_long.timestep = parse.(Int, df_long.timestep)
        df_long.plant_id .= p
        weights_df = vcat(weights_df, df_long)
    end
    weights_df.wind_power = round.(weights_df.wind_power, digits=2)
    XLSX.writetable("output/wind_ts_weights.xlsx", weights_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end


function find_and_write_inflow_weights(bernstein_degree)
    inflow_df = DataFrame(XLSX.readtable("output/inflow_data.xlsx", "Sheet1", infer_eltypes=true))
    weights_df = DataFrame()
    P = unique(inflow_df.plant_id)
    for p in P
        plant_inflow_ts = inflow_df[inflow_df.plant_id .== p, :inflow]
        parameters = define_parameters(plant_inflow_ts, bernstein_degree, 100)
        df = define_problem(parameters)
        df.b = 0:bernstein_degree
        df_long = stack(df, Not(:b), variable_name="timestep", value_name="inflow")
        df_long.timestep = parse.(Int, df_long.timestep)
        df_long.plant_id .= p
        weights_df = vcat(weights_df, df_long)
    end
    weights_df.inflow = round.(weights_df.inflow, digits=2)
    XLSX.writetable("output/inflow_weights.xlsx", weights_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end


nB = 3
# find_and_write_demand_weights(nB)
# find_and_write_wind_weights(nB)
# find_and_write_inflow_weights(nB)


# find_and_write_capacity_weights()
# find_bernstein_weights()
