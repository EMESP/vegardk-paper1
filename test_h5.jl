using HDF5



# Recursively list all groups and datasets in the file
function explore_h5(file, path="/")
    # Get the contents at the current path
    for name in keys(file[path])
        println(name)
        current_path = joinpath(path, name)
        println(current_path)
        
        # If it's a group, recursively explore it
        if isa(file[current_path], HDF5.Group)
            explore_h5(file, current_path)
        end
    end
end

path = "C:/Users/vegardvk/vscodeProjects/res100-dataset/datasets/4area"
filename = "TidsserieData.h5"
filepath = join([path, filename], "/")
println(filepath)
file = h5open(filepath, "r")
explore_h5(file)
# h5open(join([folder, filename], "/"), "r") do file
#     data = read(file, "dataset_name")
#     println(data)
# end