## Code for scheduling in a power-intensive manufacturing plant providing interruptible load

Some notes:
1. The main.jl file can be run directly to obtain results for mixed-integer recourse with $\bar{\zeta} = 11$

2. For continuous recourse only set _recourse_type = 0_ in main.jl

3. In the compressor_train_data.jl file, the following parameters can be changed to reproduce results from the paper:
   * uncertainty set parameters
   * storage tank sizes and position of tank 2
   * number of operating modes
   * $\beta$ (the weight on $\hat{C}$) and $\bar{\zeta}$
   * the day of operation, that is, electricity and interruptible load prices
         
4. In the params_recourse_type.jl file:
   * breakpoint settings can be changed
   * no other parameters should need tuning in this file

5. To run a different process network:
   * change the compressor_train_data.jl file 
   * params_recourse_type.jl, robust_counterpart.jl, and worst_case.jl need not be modified
   * solution_and_plots.jl has to be modified naturally
