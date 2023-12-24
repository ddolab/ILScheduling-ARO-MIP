##### All the parameters in the uncertain scheduling problem are listed here #####
##### that is, parameters for Problem (11) #####

###### Case study specific parameters ###############
## Compressor parameters are specified here ##
# Input flowrate of CO2 into the compressor train
input_flowrate = 63 # x 1000 m3/h or equivalently, 125 tons/h at 1 bar 

# utilization is the ratio of input flowrate and maximum flowrate in the compressor
utilization = 0.9 

# turndown is (max flowrate - min flowrate)/max flowrate in the compressor
turndown = 0.2

# using the turndown and utilization values, limits on compressor capacity are calculated as:
max_flowrate = input_flowrate/utilization
min_flowrate = max_flowrate*(1-turndown)

startup_max_flowrate = 5 # x 1000 m3/h
startup_min_flowrate = 5 # x 1000 m3/h

shutdown_max_flowrate = 5 # x 1000 m3/h
shutdown_min_flowrate = 5 # x 1000 m3/h

###### Parameter specific to our uncertain scheduling problem ########
# Scheduling horizon
T = 1:24
tfin = T[end]

## Parameters for the plant model
# sets of processes and products
I = 1:5 # for a five compressor train
J = [:C1, :C2, :C3, :C4, :C5, :C6] # labels for each CO2 stream in the compressor train

# density of CO2 in each stream
density = Dict(j => 1.87 for j in J) # kg/m3
density[:C2] = 4.5
density[:C3] = 10.77
density[:C4] = 18.58
density[:C5] = 60

# sizes of storage tanks in the compressor train
# note that the volume of the tanks is fixed even when tank 2 is moved along the compressor train
tank1_size = 110 # x 1000 m3 or equivalently, 210 tons at 1 bar
tank2_size = 40 # x 1000 m3 or equivalently, 430 tons at 5.8 bar when tank 2 is between compressors 2 and 3
tank2_position = J[3] # can take values from 2 to 5, '3' corresponds to tank 2 between compressors 2 and 3

# sets of processes that consume a product (Ibar) and those that produce a product (Ihat)
Ihat = Dict()
Ibar = Dict()
Ihat[:C1] = []
Ibar[:C1] = [:1]
Ihat[:C2] = [:1]
Ibar[:C2] = [:2]
Ihat[:C3] = [:2]
Ibar[:C3] = [:3]
Ihat[:C4] = [:3]
Ibar[:C4] = [:4]
Ihat[:C5] = [:4]
Ibar[:C5] = [:5]
Ihat[:C6] = [:5]
Ibar[:C6] = []

# inventory limits and initial conditions
IVmin = Dict((j, t) => 0.0 for j in J for t in T)
IVmax = Dict((j, t) => 0.0 for j in J for t in T)
for t in T
    IVmax[:C1, t] = tank1_size
    IVmax[tank2_position, t] = tank2_size*density[tank2_position]/density[:C1]
    # note that, we convert all volumes to 1 bar pressure volumes from here on
    # necessary adjustments in power consumption calculations are made automatically
    IVmax[:C6, t] = 5000
end
IVini = Dict(j => 0.0 for j in J)
IVini[:C1] = tank1_size

# we can vary the number of operating modes in each compressor here
number_modes = 3 # other inputs are 2 modes (on and off) and 4 modes (on, off, startup, and shutdown) 

# labels for the different operating modes in a process
M = Dict()
for i in I
    if number_modes == 2
        M[i] = [:off1, :on]
    elseif number_modes == 3
        M[i] = [:off1, :startup, :on]
    elseif number_modes == 4
        M[i] = [:off1, :startup, :on, :shutdown]
    else println("The number of modes in a compressor can only take the values 2,3, or 4. \n 
        Please change the input for 'number_modes'")
            break
    end
end

# If convex region surrogate (CRS) models are used to formulate the processes, then we may have subregions in each operating mode which are stored in the set R
# In the compressor train case study, we do not have multiple subregions for each operating mode.
R = Dict((i, m) => [1] for i in I for m in M[i])

# The set F again comes from CRS model
# It's used to list the hyperplanes necessary to define the process model in each subregion 'r' of the opreating mode 'm'
F = Dict()
for i in I
    F[i, :off1, 1] = []
    F[i, :on, 1] = 1:4 
    if number_modes >= 3
        F[i, :startup, 1] = 1:4
    end
    if number_modes == 4
        F[i, :shutdown, 1] = 1:4
    end
end

# Here all the coefficients of the hyperplanes that define the process model in subregion 'r' of operating mode 'm' are written
# General form of hyperplane: \sum_j n[i,m,r,f,j] P[i,m,r,f,j] \geq y_[i,m,r,f] \sum_j n[i,m,r,f,j] a[i,m,r,f,j] 
a = Dict((i, m, r, f, j) => 0.0 for i in I for m in M[i] for r in R[i,m] for f in F[i, m, r] for j in J)
n = Dict((i, m, r, f, j) => 0.0 for i in I for m in M[i] for r in R[i,m] for f in F[i, m, r] for j in J)

for j in J
    for i in Ibar[j]
        a[i, :on, 1, 1, j] = min_flowrate
        n[i, :on, 1, 1, j] = 1
        a[i, :on, 1, 2, j] = max_flowrate
        n[i, :on, 1, 2, j] = -1
        n[i, :on, 1, 3, j] = 1
        n[i, :on, 1, 4, j] = -1
    end
    for i in Ihat[j]
        n[i, :on, 1, 3, j] = -1
        n[i, :on, 1, 4, j] = 1
    end
end

if number_modes >= 3
    for j in J
        for i in Ibar[j]
            a[i, :startup, 1, 1, j] = startup_min_flowrate
            n[i, :startup, 1, 1, j] = 1
            a[i, :startup, 1, 2, j] = startup_max_flowrate
            n[i, :startup, 1, 2, j] = -1
            n[i, :startup, 1, 3, j] = 1
            n[i, :startup, 1, 4, j] = -1
        end
        for i in Ihat[j]
            n[i, :startup, 1, 3, j] = -1
            n[i, :startup, 1, 4, j] = 1
        end
    end
end

if number_modes == 4
    for j in J
        for i in Ibar[j]
            a[i, :shutdown, 1, 1, j] = shutdown_min_flowrate
            n[i, :shutdown, 1, 1, j] = 1
            a[i, :shutdown, 1, 2, j] = shutdown_max_flowrate
            n[i, :shutdown, 1, 2, j] = -1
            n[i, :shutdown, 1, 3, j] = 1
            n[i, :shutdown, 1, 4, j] = -1
        end
        for i in Ihat[j]
            n[i, :shutdown, 1, 3, j] = -1
            n[i, :shutdown, 1, 4, j] = 1
        end
    end
end

# Also required parameters in the CRS model
# Maximum production rate possible in subregion 'r' of operating mode 'm' in process 'i' for material 'j'
PDmax = Dict((i, m, r, j) => 0.0 for i in I for m in M[i] for r in R[i,m] for j in J)
for j in J
    for i in Ihat[j]
        PDmax[i, :on, 1, j] = max_flowrate
    end
    for i in Ibar[j]
        PDmax[i, :on, 1, j] = max_flowrate
    end
end
if number_modes >= 3
    for j in J
        for i in Ihat[j]
            PDmax[i, :startup, 1, j] = startup_max_flowrate
        end
        for i in Ibar[j]
            PDmax[i, :startup, 1, j] = startup_max_flowrate
        end
    end
end
if number_modes == 4
    for j in J
        for i in Ihat[j]
            PDmax[i, :shutdown, 1, j] = shutdown_max_flowrate
        end
        for i in Ibar[j]
            PDmax[i, :shutdown, 1, j] = shutdown_max_flowrate
        end
    end
end

# electricity consumption coefficients for each process wrt each product
# EC = \sum_i \sum_m \sum_r (δ[i,m,r]*y[i,m,r] + \sum_j γ[i,m,r,j] PD[i,m,r,j])
δ = Dict((i, m, r) => 0.0 for i in I for m in M[i] for r in R[i,m])
γ = Dict((i, m, r, j) => 0.0 for i in I for m in M[i] for r in R[i,m] for j in J)
for j in J
    for i in Ibar[j]
        # this parameter is obtained from a isentropic single stage compressor relation
        γ[i, :on, 1, j] = 1.2 
    end
end
if number_modes >= 3
    for i in I
        δ[i, :startup, 1] = 5
    end
end
if number_modes == 4
    for i in I
        δ[i, :shutdown, 1] = 5
    end
end

# possible mode transitions and sequences
TR = Dict()
for i in I
    if number_modes == 2
        TR[i, 1] = [:off1, :on]
        TR[i, 2] = [:on, :off1]
    end
    if number_modes == 3
        TR[i, 1] = [:off1, :startup]
        TR[i, 2] = [:startup, :on]
        TR[i, 3] = [:on, :off1]
    end
    if number_modes == 4
        TR[i, 1] = [:off1, :startup]
        TR[i, 2] = [:startup, :on]
        TR[i, 3] = [:on, :shutdown]
        TR[i, 4] = [:shutdown, :off1]
    end
end
SQ = Dict() # no sequences in the compressor train case study

# mode switching information required for transition constraints
θ = Dict()
for i in I
    if number_modes == 2
        θ[i, :off1, :on] = 1
        θ[i, :on, :off1] = 1
    end
    if number_modes == 3
        θ[i, :off1, :startup] = 1
        θ[i, :startup, :on] = 2
        θ[i, :on, :off1] = 1
    end
    if number_modes == 4
        θ[i, :off1, :startup] = 1
        θ[i, :startup, :on] = 2
        θ[i, :on, :shutdown] = 1
        θ[i, :shutdown, :off1] = 2
    end
end
θ_ = Dict()
for i in I
    if number_modes == 3
        θ_[i, :off1, :startup, :on] = 2
    end
    if number_modes == 4
        θ_[i, :off1, :startup, :on] = 2
        θ_[i, :on, :shutdown, :off1] = 2
    end
end
θmax = max(
    maximum(θ[i, m, m_] for ((i, j), (m, m_)) in TR),
    maximum(θ_[i, m, m_, m__] for ((i, j), (m, m_, m__)) in SQ; init = 0),
)

# # no ramping constraints are assumed for the compressor train case study
# # ramping up and down parameters for each process/compressor
# ru = Dict((i, m, j) => PDmax[i, m, 1, j] for i in I for m in M[i] for j in J)
# rd = Dict((i, m, j) => 0.0 for i in I for m in M[i] for j in J)


# limits on amoutn of interruptible load that can be provided
ILmin = 0 * ones(T[end])
ILmax = 500 * ones(T[end])

# Inital conditions 
yini = Dict()
for i in I
    yini[i, :off1] = 0
    yini[i, :startup] = 0
    yini[i, :on] = 1
    yini[i, :shutdown] = 0
end
zini = Dict()
for ((i, j), (m, m_)) in TR
    for t = -θmax:0
        zini[i, m, m_, t] = 0
    end
end

# β can be adjusted to control the impact pf recourse decisions on the objective function
# β = 0 ⟹ Objective is nominal TC
# β = 1 ⟹ Objective is worst-case TC
β = 1

# zeta_ describes the number of previous time periods to be considered in the decision rule
zeta_ = 11

# parameters for the uncertainty set
# budget parameter
Γ = Dict(t => t for t in T)
for t = 1:8
    for k = 0:5
        Γ[t+8*k] = k + 1
    end
end
# there can be at most ν load reduction requests in μ+1 time periods
μ = 3-1
ν = 2

# Demand
Demand = Dict((j,t) => 0.0 for j in J for t in T)
for t in T
    Demand[:C1, t] = -input_flowrate # CO2 flowing into the first compressor in each time period t
end
Demand[:C6, T[end]] = T[end]*input_flowrate # all the CO2 has tobe compressed by the end of scheduling horizon


# importing electricity and interruptible load prices
xf = XLSX.readxlsx("./Prices_from_ERCOT/Prices.xlsx")
operating_day = xf["high suitability"]["D25"]
alphaEC = xf["high suitability"]["A2:A25"]
alphaIL = xf["high suitability"]["B2:B25"] # $/MWh

println("Operating on ", operating_day)