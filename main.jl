##### loading packages #####
using JuMP, Gurobi, XLSX, Plots

recourse_type = 0 # 0 for continuous recourse only and 1 for mixed-integer recourse


# any parameter changes with respect to the process network should be changed here
include("compressor_train_data.jl")
#######################################

include("utilities.jl")

# breakpoints can be changed in the params_recourse_type.jl file
# Currenlty, one breakpoint is placed for each w_t at 0.2
if recourse_type == 1
    include("params_mixed_integer_recourse.jl")
else
    include("params_continuous_recourse.jl")
end

# The robust counterpart finds:
#   1. First stage variables (amount of IL to be provided, nominal production decisions)
#   2. Coefficients in the decision rules
include("robust_counterpart.jl")
# This completes the problem
#######################################################################################

# To reproduce plots in the paper, the following files are used:
# In Fig. 7, the scenario used is generated using
load_reduction_requests = 3 # total load reduction requests have be to >= load_reduction_requests
include("worst_case.jl")

# All the variables are now in the generalized model format
# We now convert the variables into the uncertain scheduling problem format,
# summarize the results in a excel sheet, and make plots
include("solution_and_plots.jl")