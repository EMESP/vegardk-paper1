using DataFrames
using XLSX
using CSV
using Interpolations
using Statistics


include("C:/Users/vegardvk/vscodeProjects/bernstein/get_hydro_data.jl")



function process_plant_data(steps_per_hour)
    input_df = CSV.read("input/gen.csv", DataFrame)
    input_df[!, :gen_index] = 1:nrow(input_df)

    wind_df = filter(row -> row.Category == "Wind", input_df)
    
    P_t = 1:10
    P_w = wind_df[:, "gen_index"]
    P_a = Dict(
        1 => P_t,
        2 => P_w,
        # 3 => [i for i in 21:30],
    )
    output_df = DataFrame(plant_id=Int[], area=Int[], gen_ub=Float64[], gen_lb=Float64[], fuel_type=String[], fuel_price=Float64[])
    for a in 1:2
        for p in P_a[a]
            p_max = input_df[input_df.gen_index .== p, "PMax MW"][1]
            p_min = input_df[input_df.gen_index .== p, "PMin MW"][1]
            fuel_type = input_df[input_df.gen_index .== p, "Category"][1]
            fuel_price = input_df[input_df.gen_index .== p, "Fuel Price \$/MMBTU"][1]
            if fuel_type != "Wind"
                fuel_type = "Thermal"
            end
            row = (p, a, p_max, p_min, fuel_type, fuel_price)
            push!(output_df, row)
        end
    end

    hydro_df = process_hydro_data(steps_per_hour)
    for hydro_row in eachrow(hydro_df)
        output_row = (hydro_row["plant_id"], 3, hydro_row["kap_gen_mw"], 0, "Hydro", hydro_row["fuel_price"])
        push!(output_df, output_row)
    end
    XLSX.writetable("output/plant_data.xlsx", output_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
end

function process_wind_ts_data(steps_per_hour)
    input_df = CSV.read("input/gen.csv", DataFrame)
    input_df[!, :gen_index] = 1:nrow(input_df)
    wind_df = filter(row -> row.Category == "Wind", input_df)

    wind_ts_df = get_wind_ts()
    name_mapping = Dict(Symbol(wind_df[!, "GEN UID"][i]) .=> Symbol(wind_df[!, "gen_index"][i]) for i in 1:nrow(wind_df))
    rename!(wind_ts_df, name_mapping)
    wind_ts_df.timestep = 1:nrow(wind_ts_df)
    df_long = stack(wind_ts_df, Not(:timestep), variable_name="plant_id", value_name="wind_power")
    df_long.plant_id = parse.(Int, df_long.plant_id)
    P = unique(df_long.plant_id)
    extended_df = DataFrame()
    if steps_per_hour != 1
        for p in P        
            extended_plant_df = DataFrame()
            extended_plant_df.wind_power = extend_ts(df_long[df_long.plant_id .== p, :wind_power], steps_per_hour)
            extended_plant_df.timestep = 1:(24*steps_per_hour)
            extended_plant_df.plant_id .= p
            append!(extended_df, extended_plant_df)
        end
        XLSX.writetable("output/wind_ts_data.xlsx", extended_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
    else
        XLSX.writetable("output/wind_ts_data.xlsx", df_long, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
    end
end

function get_wind_ts(year=2020, month=1, day=1)
    file_path = "input/DAY_AHEAD_wind.csv"
    df = CSV.read(file_path, DataFrame)
    filtered_df = subset(df, :Year => ByRow(==(year)), :Month => ByRow(==(month)), :Day => ByRow(==(day)))
    select!(filtered_df, Not([:Year, :Month, :Day, :Period]))
    multiplier = 0.3
    filtered_df .= filtered_df .*multiplier
    return filtered_df 
end

function process_hydro_data(steps_per_hour)
    file_path = "C:/Users/vegardvk/vscodeProjects/bernstein/Input/moduldata.xlsx"
    # xf = XLSX.readxlsx(file_path)
    # sh = xf["Sheet1"]

    df = DataFrame(XLSX.readtable(file_path, "Sheet1", infer_eltypes=true))
    df = df[df.vassdrag .== "NIDELVA_H.DETD", :]
    df = select(df, :modnr, :modnavn, :kap_mag_mm3, :kap_gen_m3s, :kap_forb_m3s, :kap_gen_mw, 
                :enekv, :topo_gen, :topo_forb, :topo_flom, :tilsig_reg_mm3, :tilsig_ureg_mm3)
    df[!, "kap_spill"] .= 100000
    rename!(df, :modnr => "plant_id", :kap_mag_mm3 => "kap_mag", :kap_forb_m3s => "kap_forb")
    df.kap_mag .*= 1000
    df.kap_forb .*= 1000

    # P = df[:, "plant_id"]
    I_disch, I_spill, I_bypass = find_connected_plants(df)

    end_reservoir = 49900
    energy_price = 10
    df[!, "fuel_price"] .= 0.0
    set_wv!(I_disch, df, end_reservoir, energy_price)
    add_start_reservoir!(df)
    write_inflow_table(df, 24 * steps_per_hour)
    XLSX.writetable("output/hydro_data.xlsx", df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
    return df
end


function process_load_data(steps_per_hour)
    file_path = "input/consumption.xlsx"
    df = DataFrame(XLSX.readtable(file_path, "Sheet1", infer_eltypes=true))
    multiplier = 200
    df.Forbruk .*= multiplier
    df = select(df, :Forbruk)
    df = df[1:24*3, :]
    df.timestep = repeat(1:24, outer=3)
    df.area = repeat(1:3, inner=24)
    df = select(df, :area, :timestep, :Forbruk)

    if steps_per_hour != 1
        A = unique(df.area)
        extended_df = DataFrame()
        for a in A        
            extended_area_df = DataFrame()
            extended_area_df.Forbruk = extend_ts(df[df.area .== a, :Forbruk], steps_per_hour)
            extended_area_df.timestep = 1:(24*steps_per_hour)
            extended_area_df.area .= a
            append!(extended_df, extended_area_df)
        end
        XLSX.writetable("output/load_data.xlsx", extended_df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
    else
        XLSX.writetable("output/load_data.xlsx", df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
    end
end

function extend_ts(ts, n_steps)
    new_ts = ts_interpolation(n_steps, ts)

    for i in 1:length(ts)
        original_avg = ts[i]
        segment = new_ts[((i-1)*n_steps + 1):i*n_steps]
        segment_avg = mean(segment)
        
        scaling_factor = original_avg/segment_avg
        new_ts[((i-1)*n_steps + 1):i*n_steps] .= segment .* scaling_factor
    end
    return new_ts

end

function ts_interpolation(n_steps, input_ts)
    steplength = 1/(n_steps * 2)
    output_vals = []

    for j in 1:(length(input_ts)-1)
        timestep = 0:1
        val = [input_ts[j], input_ts[j+1]]
        itp = LinearInterpolation(timestep, val, extrapolation_bc=Line())
        
        if j == 1
            for i in [steplength + k*2*steplength for k in (-n_steps/2):(n_steps-1)]
                push!(output_vals, max(itp(i), 0))
            end
        elseif j == (length(input_ts)-1)
            for i in [steplength + k*2*steplength for k in 0:((n_steps-1) + (n_steps/2))]
                push!(output_vals, max(itp(i), 0))
            end
        else
            for i in [steplength + k*2*steplength for k in 0:(n_steps-1)]
                push!(output_vals, max(itp(i), 0))
            end
        end
    end
    return output_vals
end

# process_hydro_data(4)
# process_wind_ts_data(4)
# ts_interpolation(4)
