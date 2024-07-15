using JuMP

# Create a model
model = Model()

@variable(model, x >= 0, Int)
@constraint(model, 2x <= 1)
# num_constraints(model; count_variable_in_set_constraints = true)

num_constraints(model, VariableRef, MOI.ZeroOne)
# num_constraints(model; count_variable_in_set_constraints = false)