include("C:/Users/vegardvk/vscodeProjects/bernstein/process_input_data.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/continuous_model.jl")
# include("C:/Users/vegardvk/vscodeProjects/bernstein/continuous_model_testing.jl")
# include("C:/Users/vegardvk/vscodeProjects/bernstein/continuous_model_exogenous_prod.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/plot_results.jl")
include("C:/Users/vegardvk/vscodeProjects/bernstein/find_bernstein_weights.jl")

steps_per_hour = 4
process_wind_ts_data(steps_per_hour)
process_load_data(steps_per_hour)
process_plant_data(steps_per_hour)

bernstein_degree = 3
find_and_write_inflow_weights(bernstein_degree)
find_and_write_wind_weights(bernstein_degree)
find_and_write_demand_weights(bernstein_degree)

find_and_write_production_weights(bernstein_degree, false)
find_and_write_shedding_weights(bernstein_degree, false)
find_and_write_dumping_weights(bernstein_degree, false)


model = define_and_solve_model()
write_results(model)
calculate_objective_components_continuous()