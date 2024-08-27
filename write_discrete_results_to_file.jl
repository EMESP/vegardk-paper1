using DataFrames
using XLSX

function write_output_to_file(model, data)
    write_plant_proudction(model, data)
    write_plant_data(model, data)
end


function write_plant_proudction(model, data)
    production = value.(model[:production])
    df = DataFrame(id=Int[], plant_id=Int[], timestep=Int[], production=Float64[])
    id_counter = 1
    for p in data["P"]
        for t in data["T"]
            row = (id_counter, p, t, production[p, t])
            push!(df, row)
            id_counter += 1
        end
    end
    XLSX.writetable("output/plant_production.xlsx", df, overwrite=true, sheetname="production", anchor_cell="A1")
end

# function write_plant_data(model, data)
#     df = DataFrame(plant_id=Int[], area=Int[], gen_ub=Float64[], gen_lb=Float64[])
#     for a in data["A"]
#         for p in data["P_a"][a]
#             row = (p, a, data["gen_ub"][p], data["gen_lb"][p])
#             push!(df, row)
#         end
#     end
#     XLSX.writetable("output/plant_data.xlsx", df, overwrite=true, sheetname="Sheet1", anchor_cell="A1")
# end