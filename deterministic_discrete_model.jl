using DataFrames
using JuMP
using XLSX
using CSV
# using Ipopt
using GLPK
using CPLEX
using MathOptInterface

# const MOI = MathOptInterface

include("C:/Users/vegardvk/vscodeProjects/bernstein/helper_functions.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/write_discrete_results_to_file.jl")

function define_parameters()
    P = 1:20
    T = 1:24
    A = 1:3
    L = 1:3

    n_areas = 3

    L_in = Dict(
        1 => [],
        2 => [1],
        3 => [2, 3]
    )
    L_out = Dict(
        1 => [1, 2],
        2 => [3],
        3 => []
    )
    C = Dict(1 => 1, 2 => 2, 3=> 3, 4=> 1000)
    # C_shedding = 10000
    # C_dumping = 100
    # C_startup = 100

    cap_line = Dict(
        1 => 50,
        2 => 50,
        3 => 50
    )
    load_df = get_load(n_areas, T[end])
    inflow_df = get_inflow(T[end])
    parameter_dict = Dict(
        "L_in" => L_in,
        "L_out" => L_out,
        "cap_line" => cap_line,
        "load" => load_df,
        "inflow" => inflow_df,
        "P" => P,
        "L" => L,
        "T" => T,
        "A" => A,
    )

    add_power_plant_data!(parameter_dict)
    get_module_data!(parameter_dict)
    return parameter_dict
end


function define_model()
    model = Model()
    p_dict = define_parameters()
    P = p_dict["P"]
    L = p_dict["L"]
    T = p_dict["T"]
    A = p_dict["A"]
    P_h = p_dict["P_h"]# Set of hydropower plants, subset of p
    P_t = p_dict["P_t"]
    P_w = p_dict["P_w"]
    I_bypass = p_dict["I_bypass"] # Dictionary where I_bypass[p] gives set of powerplants that bypass into p
    I_disch = p_dict["I_disch"]
    I_spill = p_dict["I_spill"]

    println(P)
    # set_optimizer(model, Ipopt.Optimizer)
    # set_optimizer(model, GLPK.Optimizer)
    set_optimizer(model, CPLEX.Optimizer)


    @variable(model, production[p in P, t in T] >= 0)

    @variable(model, load_shedding[a in A, t in T] >= 0)
    @variable(model, power_dumping[a in A, t in T] >= 0)

    @variable(model, transmission[l in L, t in T])

    # @variable(model, 0 ≤ x1[p in P, t in T] ≤ 1)
    # @variable(model, 0 ≤ x2[p in P, t in T] ≤ 1)
    # @variable(model, startup[p in P, t in T] ≥ 0)
    # @constraint(model, production_linking[p in P, t in T], production[p, t] == x1[p, t] * p_dict["gen_lb"][p] + x2[p, t] * (p_dict["gen_ub"][p] - p_dict["gen_lb"][p]))
    # @constraint(model, start_before_prod[p in P, t in T], x1[p, t] ≥ x2[p, t])
    # @constraint(model, startup_count[p in P, t in T[2:end]], x1[p, t] - x1[p, t-1]  == startup[p, t])

    @variable(model, status[p in P, t in T], Bin)
    @variable(model, startup[p in P, t in T], Bin)
    @constraint(model, gen_ub[p in P, t in T], production[p, t] ≤ p_dict["gen_ub"][p]*status[p, t])
    @constraint(model, wind_prod[p in P_w, t in T], production[p, t] == p_dict["wind_ts"][t, "$p"])
    @constraint(model, gen_lb[p in P, t in T], production[p, t] ≥ p_dict["gen_lb"][p]*status[p, t])
    @constraint(model, gen_on_off[p in P, t in T[2:end]], startup[p, t] ≥ status[p, t] - status[p, t-1])

    @variable(model, flow_disch[p in P_h, t in T] ≥ 0) # Antar at alle moduler bare har ett discharge-segment
    @variable(model, flow_bypass[p in P_h, t in T] ≥ 0)
    @variable(model, flow_spill[p in P_h, t in T] ≥ 0)
    @variable(model, total_flow_in[p in P_h, t in T] ≥ 0)
    @variable(model, total_flow_out[p in P_h, t in T] ≥ 0)
    @variable(model, volume[p in P_h, t in 0:T[end]] ≥ 0)

    @constraint(model, controlled_inflow[p in P_h, t in T], total_flow_in[p, t] == sum(flow_disch[i, t] for i in I_disch[p]) 
                                                                                + sum(flow_bypass[i, t] for i in I_bypass[p]) 
                                                                                + sum(flow_spill[i, t] for i in I_spill[p])
                                                                                + p_dict["inflow"][t, "$p"])


    @constraint(model, controlled_outflow[p in P_h, t in T], total_flow_out[p, t] == flow_disch[p, t] + flow_bypass[p, t] + flow_spill[p, t])
    @constraint(model, starting_reservoir[p in P_h], volume[p, 0] == p_dict["starting_reservoir"][p])
    @constraint(model, reservoir_balance[p in P_h, t in 1:(T[end])], volume[p, t] - volume[p, t-1] == total_flow_in[p, t] - total_flow_out[p, t])

    @constraint(model, hydro_production[p in P_h, t in T], production[p, t] == p_dict["enekv"][p] * flow_disch[p, t])
    @constraint(model, vol_ub[p in P_h, t in T], volume[p, t] ≤ p_dict["kap_mag"][p])
    @constraint(model, bypass_ub[p in P_h, t in T], flow_bypass[p, t] ≤ p_dict["kap_forb"][p])
    @constraint(model, prod_ub[p in P_h, t in T], flow_disch[p, t] ≤ p_dict["kap_gen"][p])
    @constraint(model, spill_ub[p in P_h, t in T], flow_spill[p, t] ≤ p_dict["kap_spill"][p])

    @constraint(model, transmission_ub[l in L, t in T], transmission[l, t] ≤ p_dict["cap_line"][l])
    @constraint(model, transmission_lb[l in L, t in T], transmission[l, t] ≥ -p_dict["cap_line"][l])

    @constraint(model, energy_balance[a in A, t in T], sum(production[p, t] for p in p_dict["P_a"][a]) 
                + load_shedding[a, t] - power_dumping[a, t] 
                + sum(transmission[l, t] for l in p_dict["L_in"][a]) - sum(transmission[l, t] for l in p_dict["L_out"][a]) == p_dict["load"][t, a])
    println(p_dict["P_t"])
    @objective(model, Min, sum(production[p, t] * p_dict["fuel_price"][p] for p in P_t for t in T)
                + sum((-volume[p, T[end]] + volume[p, 0]) * p_dict["fuel_price"][p] for p in P_h) 
                + sum(load_shedding[a, t] * p_dict["C_shedding"] for a in A for t in T) 
                + sum(power_dumping[a, t] * p_dict["C_dumping"] for a in A for t in T)
                + sum(startup[p, t] * p_dict["C_startup"] for p in P_t for t in T))
    # println(model)
    print_model_info(model)
    optimize!(model)

    # println(objective_value(model))

    return model, p_dict
end

function print_simple_results(model)
    p_dict = define_parameters()
    A = p_dict["A"]
    L = p_dict["L"]
    P = p_dict["P"]
    T = p_dict["T"]
    
    P = p_dict["P"]
    P_h = p_dict["P_h"]
    P_t = p_dict["P_t"]

    load_shedding = value.(model[:load_shedding])
    # areas = length(load_shedding[:, 1])
    println("Load shedding")
    for a in A
        println("\tArea $a: ", sum(load_shedding[a, :]))
    end
    println("Power dumping:")
    power_dumping = value.(model[:power_dumping])
    for a in A
        println("\tArea $a: ", sum(power_dumping[a, :]))
    end
    
    transmission = value.(model[:transmission])
    # lines = length(transmission[:, 1])
    println("Transmission:")
    for l in L
        println("\tLine $l: ")
        println("\t\tPositive direction: ", sum([t for t in transmission[l, :] if t>0]))
        println("\t\tNegative direction: ", sum([t for t in transmission[l, :] if t<0]))
    end

    production = value.(model[:production])
    # power_plants = length(production[:, 1])
    println("Production")
    for p in p_dict["P"]
        println("\t Power plant $p: ", sum(production[p, :]))
    end
    # a = 3

    load_shedding = value.(model[:load_shedding])
    power_dumping = value.(model[:power_dumping])
    startup = value.(model[:startup])
    volume = value.(model[:volume])

    prod_costs      = sum([production[p, t] * p_dict["fuel_price"][p] for p in P_t for t in T])
    startup_costs   = sum([startup[p, t] * p_dict["C_startup"] for p in P_t for t in T])
    shed_costs      = sum(load_shedding[a, t] * p_dict["C_shedding"] for a in A for t in T)
    dumping_costs   = sum(power_dumping[a, t] * p_dict["C_dumping"] for a in A for t in T)
    volume_costs    = sum((-volume[p, T[end]] + volume[p, 0]) * p_dict["fuel_price"][p] for p in P_h) 

    objective = objective_value(model)

    println("Objective: $objective")
    println("\t Production costs: $prod_costs")
    println("\t Volume costs: $volume_costs")
    println("\t Startup costs: $startup_costs")
    println("\t shed_costs: $shed_costs")
    println("\t dump costs: $dumping_costs")

end

function plot_energy_balance(model)
    p_dict = define_parameters()

    areas = last(p_dict["A"])
    lines = last(p_dict["L"])
    power_plants = last(p_dict["P"])
    time_steps = last(p_dict["T"])

    for a in 1:areas
        df = DataFrame()
        transmission = value.(model[:transmission])

        imported_power = zeros(time_steps)
        exported_power = zeros(time_steps)

        for l in p_dict["L_in"][a]
            exported_power -= [min(0, val) for val in collect(transmission[l, :])]
            imported_power += [max(0, val) for val in collect(transmission[l, :])]
        end

        for l in p_dict["L_out"][a]
            exported_power += [max(0, val) for val in collect(transmission[l, :])]
            imported_power -= [min(0, val) for val in collect(transmission[l, :])]
        end
        # df[!, "Imported power"] = imported_power
        # df[!, "Exported power"] = exported_power
        df[!, "Net position"] = exported_power - imported_power
        df[!, "Demand"] = p_dict["load"][1:time_steps, a]

        production = zeros(time_steps)
        prod_df = DataFrame()
        for p in p_dict["P_a"][a]
            prod_df[!, "Power plant $p"] = collect(value.(model[:production])[p, :])
            production += collect(value.(model[:production])[p, :])
        end
        df[!, "Production"] = production
        # df = hcat(df, prod_df)
        # load_df = p_dict["load"])
        df[!, "Shedding"] = collect(value.(model[:load_shedding]))[a, :]
        df[!, "Dumping"] = collect(value.(model[:power_dumping]))[a, :]
        # display(df)
        plot_all_columns_df(df, "Discrete model, area $a", "discrete_area$a.png")
    end
end

function plot_hydro_balance(model)
    p_dict = define_parameters()
    T = p_dict["T"]

    volume = value.(model[:volume])
    total_flow_in = value.(model[:total_flow_in])
    total_flow_out = value.(model[:total_flow_out])

    for p in p_dict["P_h"]
        df = DataFrame()

        
        inflow = collect(total_flow_in[p, :])
        pushfirst!(inflow, 0)
        outflow = collect(total_flow_out[p, :])
        pushfirst!(outflow, 0)

        df[!, "volume"] = collect(volume[p, :])
        df[!, "total_flow_in"] = inflow
        df[!, "total_flow_out"] = outflow
        plot_all_columns_df(df, "Discrete model, hydro balance, plant $p", "discrete_hydro_balance$p.png")

    end
end

function plot_unit_commitment(model)
    p_dict = define_parameters()
    T = p_dict["T"]
    thermal_uc = value.(model[:status])

    for a in p_dict["A"]
        df = DataFrame()
        peak = zeros(size(thermal_uc, 2))
        for p in p_dict["P_a"][a]
            # if !(p in p_dict["P_t"]) 
            #     continue
            # end
            uc = collect(thermal_uc[p, :])
            if sum(uc) != 0
                df[!, "$p"] = uc + peak
                peak .+= uc
            end
        end
        if ncol(df) > 0
            plot_all_columns_df(df, "Discrete model, unit commitment, area $a", "discrete_uc_area$a.png", true)
        end
    end
end

# p_max = get_power_plant_data()
# println(p_max)
# p_dict = define_parameters()

model, data = define_model()
write_output_to_file(model, data)
# print_simple_results(model)
# plot_hydro_balance(model)
# plot_energy_balance(model)
# plot_unit_commitment(model)
