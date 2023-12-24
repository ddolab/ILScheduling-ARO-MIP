The main.jl file can be run directly to obtain results for mixed-integer recourse with \bar{\zeta} = 11

For continuous recourse only:
         set recourse_type = 0 in main.jl

In the compressor_train_data.jl file, the following parameters can be changed to reproduce results from the paper:
         uncertainty set parameters
         storage tank sizes and position of tank 2
         number of operating modes
         \beta (the weight on TC^{recourse}) and \bar{\zeta}
         the day of operation, that is, electricity and interruptible load prices
         


In the params_recourse_type.jl file:
         breakpoint settings can be changed
         no other parameters should need tuning in this file

To run a different process network:
         change the compressor_train_data.jl file 
         params_recourse_type.jl, robust_counterpart.jl, and worst_case.jl need not be modified
         solution_and_plots.jl has to be modified naturally