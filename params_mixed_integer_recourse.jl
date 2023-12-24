# The data given in the compressor_train_data.jl is converted to the general matrix form ###
#### The general form can be written as:                                                ####
#### minimize 
#### such that A_1 x_1 + D_1 y_1 ≤ b_1              (deterministic constraints)         ####
####           A_t(ξ^[t]) x_1 + D_t(ξ^[t]) y_1       (uncertain constraints)            ####
####                + \sum_{t'=2}^t [Atilde_tt' x_t'(ξ^[t']) + Dtilde_tt' y_t'(ξ^[t'])] ####
####                           ≤ b^t(ξ^[t])    \forall t \in stages, ξ^[t] \in Ξ^[t]    ####

ST = tfin + 1 # number of stages

############################ UNCERTAINTY SET ###############################################
####     Ξ^[t] = {ξ^[t]: W_t ξ^[t] ≤ u_t, ξ_11 = 1}                                     ####
# Define the number of uncertain parameters in each stage
K = Dict(st => 1 for st = 1:ST) # we only have one uncertain parameter in each stage
Num_K = Dict(st => sum(K[stt] for stt = 1:st) for st = 1:ST) # The number of uncertain parameters up to stage st

# Number of constraints in original uncertainty set in each stage
# Mu = [3t-1 for t in 1:ST]
Mu = [2t-1 for t in 1:ST]
Mu[1] = 1
Mu = trunc.(Int, Mu)

# uncertainty set parameters
W = Dict((st, stt, i) => zeros(Mu[st]) for st = 2:ST for stt = 1:st for i = 1:K[stt])
u = Dict(st => zeros(Mu[st]) for st = 2:ST)

for st in 2:ST
    for tt in 2:st
        ## w_tt ≤ 1
        W[st, 1, 1][tt-1] = -1
        W[st, tt, 1][tt-1] = 1
        ## w_tt ≥ 0
        W[st, 1, 1][tt+st-2] = -1
        # ## no. of load reductions in consecutive μ+1 time periods is less than ν
        # for ttt in max(tt-μ,2):tt
        #     W[st, ttt, 1][2*(st-1)+tt-1] = 1
        # end
        # W[st, 1, 1][2*(st-1)+tt-1] = -ν
    end
    ## cummulative load reduction 
    # W[st, 1, 1][3*st-1] = -Γ[st-1]
    # for tt in 2:st
    #     W[st, tt, 1][3*st-1] = 1
    # end
    W[st, 1, 1][2*st-1] = -Γ[st-1]
    for tt in 2:st
        W[st, tt, 1][2*st-1] = 1
    end
end

## These details are relevant to the Lifted Uncertainty Sets
# Define xi_max and xi_min to facilitate the application of general reformulation
xi_min = Dict((st, i) => 0.0 for st = 1:ST for i = 1:K[st])
xi_max = Dict((st, i) => 1.0 for st = 1:ST for i = 1:K[st])

# Define r[st, i] whose value - 1 is the number of breakpoints for xi[st, i]
r = Dict((st, i) => 2 for st = 1:ST for i = 1:K[st]) # the dimensions of ξ_bar
r[1, 1] = 1
g = Dict((st, i) => max(1, r[st,i]-1) for st = 1:ST for i = 1:K[st]) # the dimensions of ξ_hat

# breakpoints 
# 0.2 was obtained using the power consumption values of the compressor
breakpoints = Dict((t, i) => [0, 0.2, 1] for t = 2:ST for i = 1:K[t])
breakpoints[1, 1] = [1]

# vertices in the lifted marginal support of uncertain parameters
v_bar = vertices_generation(breakpoints, ST, K, r)
v_hat = G_generation(ST, K, r)


########################### VARIABLE INFO ##################################################
# Number of integer variables
Q = zeros(ST)
Q[1] = tfin + (θmax) * (sum(size(M[i])[1]^2 for i in I)) +
    tfin * (sum(sum(length(R[i, m]) for m in M[i]) for i in I)) # x_t, {z_imm't for t = -theta_max+1:0}, yimrt_nominal
for t = 2:ST
    Q[t] = (sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + (sum(size(M[i])[1]^2 for i in I)) #  y_imrt, z_imm't> 0
end

# Number of continuous variables
P = zeros(ST)
P[1] = 1 + tfin + sum(size(M[i])[1] for i in I) +
    tfin * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)# epigraph variable, IL_t, y_im0, PDhat
for t = 2:ST
    P[t] = size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + sum(size(M[i])[1] for i in I) # PDtilde, y_imt
end


########################## CONSTRAINTS INFO ##########################################################
zetaset = zeta_sets(1:tfin+1, zeta_)

# the number of deterministic constraints in each decision stage (basic number)
N = zeros(ST)
N[1] = 2 * tfin + # IL limits
    2 * sum(size(M[i])[1] for i in I) + # initial conditions for y_imrt
    2 * length(TR) * (θmax - 1) + # initial conditions for z_imm't
    tfin * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) # PDtilde ≥ 0

for t = 2:ST-1
    N[t] = 2 * sum(size(M[i])[1] for i in I) + # ∑_r y_imrt = 1, can operate in only one subregion
        2 * size(I)[1] + # ∑_i yimt = 1, can operate in only one operating mode
        2 * sum(size(M[i])[1] for i in I) + # transition constraint
        length(TR) + # minimum stay time in one mode
        2 * length(SQ) + # transition sequences
        2 * size(J)[1] + # mass balance
        1 + # load reduction constraint 
        (
            sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) + # hyperplane for operating in subregion r of mode m
            2 * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) # max production rate in subregion r of mode m for material j
        ) +
        2 * sum(size(M[i])[1] for i in I) + # y_mt >= 0, y_m t <= 1
        2 * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + #  y_mrt >= 0, y_mrt <= 1
        2 * sum(size(M[i])[1]^2 for i in I) # z_mm't >= 0, z_mm't <= 1            
end
N[ST] = N[ST-1] + 1 # +1 for objective function

P = trunc.(Int, P)
Q = trunc.(Int, Q)
N = trunc.(Int, N)

# deterministic model parameters
A = Dict((st, sttt, i) => zeros(N[st], P[1]) for st = 1:ST for sttt = 1:st for i = 1:K[sttt])
D = Dict((st, sttt, i) => zeros(N[st], Q[1]) for st = 1:ST for sttt = 1:st for i = 1:K[sttt])
b = Dict((st, sttt, i) => zeros(N[st]) for st = 1:ST for sttt = 1:st for i = 1:K[sttt])
A_tilde = Dict((st, stt) => zeros(N[st], P[stt]) for st = 2:ST for stt = 2:st)
D_tilde = Dict((st, stt) => zeros(N[st], Q[stt]) for st = 2:ST for stt = 2:st)


# constraint matrices
A[ST, 1, 1][N[ST], 1] = -1 # epigraph variable
A[ST, 1, 1][N[ST], 2:ST] = -alphaIL[1:ST-1] # IL_t
A[ST, 1, 1][N[ST], ST + sum(size(M[i])[1] for i in I)+1:1 + tfin + sum(size(M[i])[1] for i in I) +
    tfin * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)] = 
    [alphaEC[t]*γ[i,m,r,j] for t in 1:tfin for j in J for i in I for m in M[i] for r in R[i,m]] # PDhat
D[ST, 1, 1][N[ST], tfin + (θmax) * (sum(size(M[i])[1]^2 for i in I)) + 1: tfin + (θmax) * (sum(size(M[i])[1]^2 for i in I)) +
    tfin * (sum(sum(length(R[i, m]) for m in M[i]) for i in I))] = [alphaEC[t]*δ[i,m,r] for t in 1:tfin for i in I for m in M[i] for r in R[i,m]] # y_imrt
for t in 2:ST
    A_tilde[ST, t][N[ST], 1:size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)] = 
        [β*alphaEC[t-1]*γ[i,m,r,j] for j in J for i in I for m in M[i] for r in R[i,m]] # PDtilde
    D_tilde[ST, t][N[ST], 1:(sum(sum(length(R[i, m]) for m in M[i]) for i in I))] = 
        [β*alphaEC[t-1]*δ[i,m,r] for i in I for m in M[i] for r in R[i,m]] #ytilde_imrt
end

# 1e-1 
for t = 2:ST
    local count_cons = 1
    local count_var_y_ = 1
    for i in I
        for m in M[i]
            A_tilde[t,t][count_cons, size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count_cons] = -1
            for r in R[i,m]
                D[t,1,1][count_cons, tfin +
                (θmax - 1) * (sum(size(M[i])[1]^2 for i in I)) + (t-2)*(sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + count_var_y_] = 1
                D_tilde[t,t][count_cons, count_var_y_] = 1
                count_var_y_ += 1
            end
            count_cons += 1
        end
    end
end

# 1e-2
for t = 2:ST
    local count_cons = 1
    local count_var_y_ = 1
    for i in I
        for m in M[i]
            A_tilde[t,t][sum(size(M[i])[1] for i in I) + count_cons, size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count_cons] = 1
            for r in R[i,m]
                D[t,1,1][sum(size(M[i])[1] for i in I) + count_cons, tfin +
                (θmax - 1) * (sum(size(M[i])[1]^2 for i in I)) + (t-2)*(sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + count_var_y_] = -1
                D_tilde[t,t][sum(size(M[i])[1] for i in I) + count_cons, count_var_y_] = -1
                count_var_y_ += 1
            end
            count_cons += 1
        end
    end
end

# 1f-1
for t in 2:ST
    local count_var = 1
    local count_i = 1
    for i in I
        for m in M[i]
            A_tilde[t,t][2 * sum(size(M[i])[1] for i in I) + count_i, size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count_var] = 1
            count_var += 1
        end
        b[t,1,1][2 * sum(size(M[i])[1] for i in I) + count_i] = 1
        count_i += 1
    end
end

# 1f-2
for t in 2:ST
    local count_var = 1
    local count_i = 1
    for i in I
        for m in M[i]
            A_tilde[t,t][2 * sum(size(M[i])[1] for i in I) + size(I)[1] + count_i, size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count_var] = -1
            count_var += 1
        end
        b[t,1,1][2 * sum(size(M[i])[1] for i in I) + size(I)[1] + count_i] = -1
        count_i += 1
    end
end

M_trans = Dict()
for i in I
    local trans_var = 1
    for m in M[i]
        M_trans[i,m] = trans_var
        trans_var += 1
    end
end

# 2-1
# for t = 1, we need a different constraint as it involves zmm',0 and y_m,0
count_z = 0
count_z2 = 0
for i in I
    for m in M[i]
        TR_ = []
        TRhat = []
        for ((x,y),(m_,m__)) in TR
            if x == i
                if m_ == m
                    TRhat = append!(TRhat, [m__])
                end
                if m__ == m
                    TR_ = append!(TR_, [m_])
                end
            end
        end
        for m_ in TR_
            D_tilde[2,2][2 * sum(size(M[i])[1] for i in I) +
            2 * size(I)[1] + count_z + M_trans[i,m], (sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + count_z2 + (M_trans[i,m_]-1)*length(M[i]) + M_trans[i,m]] = 1
        end
        for m_ in TRhat
            D_tilde[2,2][2 * sum(size(M[i])[1] for i in I) +
            2 * size(I)[1] + count_z + M_trans[i,m], (sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + count_z2 + (M_trans[i,m]-1)*length(M[i]) + M_trans[i,m_]] = -1
        end
        A_tilde[2,2][2 * sum(size(M[i])[1] for i in I) +
        2 * size(I)[1] + count_z + M_trans[i,m], size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + 
        count_z + M_trans[i,m]] = -1
        A[2,1,1][2 * sum(size(M[i])[1] for i in I) +
        2 * size(I)[1] + count_z + M_trans[i,m], 1 +
        tfin + count_z + M_trans[i,m]] = 1
    end
    global count_z = count_z + size(M[i])[1]
    global count_z2 = count_z2 + size(M[i])[1]^2
end

# multistages
for t in 3:ST
    local count_z = 0
    local count_z2 = 0
    for i in I
        for m in M[i]
            TR_ = []
            TRhat = []
            for ((x,y),(m_,m__)) in TR
                if x == i
                    if m_ == m
                        TRhat = append!(TRhat, [m__])
                    end
                    if m__ == m
                        TR_ = append!(TR_, [m_])
                    end
                end
            end
            for m_ in TR_
                D_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
                2 * size(I)[1] + count_z + M_trans[i,m], (sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + count_z2 + (M_trans[i,m_]-1)*length(M[i]) + M_trans[i,m]] = 1
            end
            for m_ in TRhat
                D_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
                2 * size(I)[1] + count_z + M_trans[i,m], (sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + count_z2 + (M_trans[i,m]-1)*length(M[i]) + M_trans[i,m_]] = -1
            end
            A_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
            2 * size(I)[1] + count_z + M_trans[i,m], size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + 
            count_z + M_trans[i,m]] = -1
            A_tilde[t,t-1][2 * sum(size(M[i])[1] for i in I) +
            2 * size(I)[1] + count_z + M_trans[i,m], size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + 
            count_z + M_trans[i,m]] = 1
        end
        count_z = count_z + size(M[i])[1]
        count_z2 = count_z2 + size(M[i])[1]^2
    end
end

# 2-2
# for t = 1, we need a different constraint as it involves zmm',0 and y_m,0
count_z = 0
count_z2 = 0
for i in I
    for m in M[i]
        TR_ = []
        TRhat = []
        for ((x,y),(m_,m__)) in TR
            if x == i
                if m_ == m
                    TRhat = append!(TRhat, [m__])
                end
                if m__ == m
                    TR_ = append!(TR_, [m_])
                end
            end
        end
        for m_ in TR_
            D_tilde[2,2][2 * sum(size(M[i])[1] for i in I) +
            2 * size(I)[1] + sum(size(M[i])[1] for i in I) + count_z + M_trans[i,m], (sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + count_z2 + (M_trans[i,m_]-1)*length(M[i]) + M_trans[i,m]] = -1
        end
        for m_ in TRhat
            D_tilde[2,2][2 * sum(size(M[i])[1] for i in I) +
            2 * size(I)[1] + sum(size(M[i])[1] for i in I) + count_z + M_trans[i,m], (sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + count_z2 + (M_trans[i,m]-1)*length(M[i]) + M_trans[i,m_]] = 1
        end
        A_tilde[2,2][2 * sum(size(M[i])[1] for i in I) +
        2 * size(I)[1] + sum(size(M[i])[1] for i in I) + count_z + M_trans[i,m], size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + 
        count_z + M_trans[i,m]] = 1
        A[2,1,1][2 * sum(size(M[i])[1] for i in I) +
        2 * size(I)[1] + sum(size(M[i])[1] for i in I) + count_z + M_trans[i,m], 1 +
        tfin + count_z + M_trans[i,m]] = -1
    end
    global count_z = count_z + size(M[i])[1]
    global count_z2 = count_z2 + size(M[i])[1]^2
end

# multistages
for t in 3:ST
    local count_z = 0
    local count_z2 = 0
    for i in I
        for m in M[i]
            TR_ = []
            TRhat = []
            for ((x,y),(m_,m__)) in TR
                if x == i
                    if m_ == m
                        TRhat = append!(TRhat, [m__])
                    end
                    if m__ == m
                        TR_ = append!(TR_, [m_])
                    end
                end
            end
            for m_ in TR_
                D_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
                2 * size(I)[1] + sum(size(M[i])[1] for i in I) + count_z + M_trans[i,m], (sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + count_z2 + (M_trans[i,m_]-1)*length(M[i]) + M_trans[i,m]] = -1
            end
            for m_ in TRhat
                D_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
                2 * size(I)[1] + sum(size(M[i])[1] for i in I) + count_z + M_trans[i,m], (sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + count_z2 + (M_trans[i,m]-1)*length(M[i]) + M_trans[i,m_]] = 1
            end
            A_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
            2 * size(I)[1] + sum(size(M[i])[1] for i in I) + count_z + M_trans[i,m], size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + 
            count_z + M_trans[i,m]] = 1
            A_tilde[t,t-1][2 * sum(size(M[i])[1] for i in I) +
            2 * size(I)[1] + sum(size(M[i])[1] for i in I) + count_z + M_trans[i,m], size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + 
            count_z + M_trans[i,m]] = -1
        end
        count_z = count_z + size(M[i])[1]
        count_z2 = count_z2 + size(M[i])[1]^2
    end
end


# 3
var_y = Dict()
var_z = Dict()
count_temp = 0
count_temp2 = 0
for i in I
    var_y[i] = count_temp
    var_z[i] = count_temp2
    global count_temp += size(M[i])[1]
    global count_temp2 += size(M[i])[1]^2
end

for t in 2:ST
    count_tr = 1
    for ((i,j),(m,m_)) in TR
        A_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
        2 * size(I)[1] +
        2 * sum(size(M[i])[1] for i in I) + count_tr, size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + 
        var_y[i] + M_trans[i,m_]] = -1
        for k in 0:(θ[i,m,m_]-1)
            if t-1-k <= 0
                D[t,1,1][2 * sum(size(M[i])[1] for i in I) +
                2 * size(I)[1] +
                2 * sum(size(M[i])[1] for i in I) + count_tr, tfin + (t-2-k+θmax-1)*(sum(size(M[i])[1]^2 for i in I)) + var_z[i] + (M_trans[i,m]-1)*size(M[i])[1]+M_trans[i,m_]] = 1
            elseif t-1-k > 0
                D_tilde[t,t-k][2 * sum(size(M[i])[1] for i in I) +
                2 * size(I)[1] +
                2 * sum(size(M[i])[1] for i in I) + count_tr, (sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + var_z[i] + (M_trans[i,m]-1)*size(M[i])[1]+M_trans[i,m_]] = 1
            end
        end
        count_tr += 1
    end
end


# 4-1
for t in 2:ST
    local count_sq = 1
    for ((i,j), (m,m_,m__)) in SQ
        if t-1-θ_[i,m,m_,m__] <= 0
            D[t,1,1][2 * sum(size(M[i])[1] for i in I) +
            2 * size(I)[1] +
            2 * sum(size(M[i])[1] for i in I) +
            length(TR) + count_sq, tfin + (t-2-θ_[i,m,m_,m__] + θmax)*(sum(size(M[i])[1]^2 for i in I)) + var_z[i] + (M_trans[i,m]-1)*size(M[i])[1]+M_trans[i,m_]] = 1
        else
            D_tilde[t,t-θ_[i,m,m_,m__]][2 * sum(size(M[i])[1] for i in I) +
            2 * size(I)[1] +
            2 * sum(size(M[i])[1] for i in I) +
            length(TR) + count_sq, (sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + var_z[i] + (M_trans[i,m]-1)*size(M[i])[1]+M_trans[i,m_]] = 1
        end
        D_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
        2 * size(I)[1] +
        2 * sum(size(M[i])[1] for i in I) +
        length(TR) + count_sq, (sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + var_z[i] + (M_trans[i,m_]-1)*size(M[i])[1] + M_trans[i,m__]] = -1
        count_sq += 1
    end
end


# 4-2
for t in 2:ST
    local count_sq = 1
    for ((i,j), (m,m_,m__)) in SQ
        if t-1-θ_[i,m,m_,m__] <= 0
            D[t,1,1][2 * sum(size(M[i])[1] for i in I) +
            2 * size(I)[1] +
            2 * sum(size(M[i])[1] for i in I) +
            length(TR) + length(SQ) + count_sq, tfin + (t-2-θ_[i,m,m_,m__] + θmax)*(sum(size(M[i])[1]^2 for i in I)) + var_z[i] + (M_trans[i,m]-1)*size(M[i])[1]+M_trans[i,m_]] = -1
        else
            D_tilde[t,t-θ_[i,m,m_,m__]][2 * sum(size(M[i])[1] for i in I) +
            2 * size(I)[1] +
            2 * sum(size(M[i])[1] for i in I) +
            length(TR) + length(SQ) + count_sq, (sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + var_z[i] + (M_trans[i,m]-1)*size(M[i])[1]+M_trans[i,m_]] = -1
        end
        D_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
        2 * size(I)[1] +
        2 * sum(size(M[i])[1] for i in I) +
        length(TR) + length(SQ) + count_sq, (sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + var_z[i] + (M_trans[i,m_]-1)*size(M[i])[1] + M_trans[i,m__]] = 1
        count_sq += 1
    end
end


# 5-1
for t in 2:ST
    for k in 1:t-1
        local count_j = 1
        for j in J
            for i in Ihat[j]
                local count_mr = 1
                for m in M[i]
                    for r in R[i,m]
                        A[t,1,1][2 * sum(size(M[i])[1] for i in I) +
                        2 * size(I)[1] +
                        2 * sum(size(M[i])[1] for i in I) +
                        length(TR) +
                        2 * length(SQ) + count_j, 1 +
                        tfin +
                        sum(size(M[i])[1] for i in I) + (k-1)*size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + var_y[i] + count_mr] = -1
                        A_tilde[t,k+1][2 * sum(size(M[i])[1] for i in I) + 2 * size(I)[1] +
                        2 * sum(size(M[i])[1] for i in I) + length(TR) +
                        2 * length(SQ) + count_j, (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + var_y[i] + count_mr] = -1
                        count_mr += 1
                    end
                end
            end
            for i in Ibar[j]
                local count_mr = 1
                for m in M[i]
                    for r in R[i,m]
                        A[t,1,1][2 * sum(size(M[i])[1] for i in I) +
                        2 * size(I)[1] +
                        2 * sum(size(M[i])[1] for i in I) +
                        length(TR) +
                        2 * length(SQ) + count_j, 1 +
                        tfin +
                        sum(size(M[i])[1] for i in I) + (k-1)*size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + var_y[i] + count_mr] = 1
                        A_tilde[t,k+1][2 * sum(size(M[i])[1] for i in I) + 2 * size(I)[1] +
                        2 * sum(size(M[i])[1] for i in I) + length(TR) +
                        2 * length(SQ) + count_j, (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + var_y[i] + count_mr] = 1
                        count_mr += 1
                    end
                end
            end
            count_j += 1
        end
    end
    local count_j = 1
    for j in J
        b[t,1,1][2 * sum(size(M[i])[1] for i in I) +
        2 * size(I)[1] +
        2 * sum(size(M[i])[1] for i in I) +
        length(TR) +
        2 * length(SQ) + count_j] = IVini[j] - sum(Demand[j,k] for k in 1:(t-1)) - IVmin[j,t-1]
        # A[t, 1, 1][2 * sum(size(M[i])[1] for i in I) +
        # 2 * size(I)[1] +
        # 2 * sum(size(M[i])[1] for i in I) +
        # length(TR) +
        # 2 * length(SQ) + count_j,1 +
        # tfin +
        # sum(size(M[i])[1] for i in I) +
        # tfin * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count_j] = -1
        # b[t,1,1][2 * sum(size(M[i])[1] for i in I) +
        # 2 * size(I)[1] +
        # 2 * sum(size(M[i])[1] for i in I) +
        # length(TR) +
        # 2 * length(SQ) + count_j] = -sum(Demand[j,k] for k in 1:(t-1)) - IVmin[j,t-1]
        count_j += 1
    end
end

# 5-2
for t in 2:ST
    for k in 1:t-1
        local count_j = 1
        for j in J
            for i in Ihat[j]
                local count_mr = 1
                for m in M[i]
                    for r in R[i,m]
                        A[t,1,1][2 * sum(size(M[i])[1] for i in I) +
                        2 * size(I)[1] +
                        2 * sum(size(M[i])[1] for i in I) +
                        length(TR) +
                        2 * length(SQ) + size(J)[1] + count_j, 1 +
                        tfin +
                        sum(size(M[i])[1] for i in I) + (k-1)*size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + var_y[i] + count_mr] = 1
                        A_tilde[t,k+1][2 * sum(size(M[i])[1] for i in I) + 2 * size(I)[1] +
                        2 * sum(size(M[i])[1] for i in I) + length(TR) +
                        2 * length(SQ) + size(J)[1] + count_j, (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + var_y[i] + count_mr] = 1
                        count_mr += 1
                    end
                end
            end
            for i in Ibar[j]
                local count_mr = 1
                for m in M[i]
                    for r in R[i,m]
                        A[t,1,1][2 * sum(size(M[i])[1] for i in I) +
                        2 * size(I)[1] +
                        2 * sum(size(M[i])[1] for i in I) +
                        length(TR) +
                        2 * length(SQ) + size(J)[1] + count_j, 1 +
                        tfin +
                        sum(size(M[i])[1] for i in I) + (k-1)*size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + var_y[i] + count_mr] = -1
                        A_tilde[t,k+1][2 * sum(size(M[i])[1] for i in I) + 2 * size(I)[1] +
                        2 * sum(size(M[i])[1] for i in I) + length(TR) +
                        2 * length(SQ) + size(J)[1] + count_j, (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + var_y[i] + count_mr] = -1
                        count_mr += 1
                    end
                end
            end
            count_j += 1
        end
    end
    local count_j = 1
    for j in J
        b[t,1,1][2 * sum(size(M[i])[1] for i in I) +
        2 * size(I)[1] +
        2 * sum(size(M[i])[1] for i in I) +
        length(TR) +
        2 * length(SQ) + size(J)[1] + count_j] = -IVini[j] + sum(Demand[j,k] for k in 1:(t-1)) + IVmax[j,t-1]
        # A[t,1,1][2 * sum(size(M[i])[1] for i in I) + 2 * size(I)[1] +
        # 2 * sum(size(M[i])[1] for i in I) + length(TR) +
        # 2 * length(SQ) + size(J)[1] + count_j, 1 +
        # tfin +
        # sum(size(M[i])[1] for i in I) +
        # tfin * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count_j] = 1
        # b[t,1,1][2 * sum(size(M[i])[1] for i in I) + 2 * size(I)[1] +
        # 2 * sum(size(M[i])[1] for i in I) + length(TR) +
        # 2 * length(SQ) + size(J)[1] + count_j] = sum(Demand[j,k] for k in 1:(t-1)) + IVmax[j,t-1]
        count_j += 1
    end
end


# 7
Omega = maximum(maximum(maximum(dECmax(i, m, r, J, F, a, n, γ) for r in R[i,m]) for m in [:on]) for i in I)
for t in 2:ST
    A[t,t,1][2 * sum(size(M[i])[1] for i in I) +
    2 * size(I)[1] +
    2 * sum(size(M[i])[1] for i in I) +
    length(TR) +
    2 * length(SQ) +
    2 * size(J)[1] + 1, 1+t-1] = 1
    b[t,1,1][2 * sum(size(M[i])[1] for i in I) +
    2 * size(I)[1] +
    2 * sum(size(M[i])[1] for i in I) +
    length(TR) +
    2 * length(SQ) +
    2 * size(J)[1] + 1] = Omega
    D[t,1,1][2 * sum(size(M[i])[1] for i in I) +
    2 * size(I)[1] +
    2 * sum(size(M[i])[1] for i in I) +
    length(TR) +
    2 * length(SQ) +
    2 * size(J)[1] + 1, t-1] = Omega
    local count_var = 1
    for j in J
        for i in I
            for m in M[i]
                for r in R[i,m]
                    A_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
                    2 * size(I)[1] +
                    2 * sum(size(M[i])[1] for i in I) +
                    length(TR) +
                    2 * length(SQ) +
                    2 * size(J)[1] + 1, count_var] = γ[i,m,r,j]
                    count_var += 1
                end
            end
        end
    end
    local count_imr = 1
    for i in I
        for m in M[i]
            for r in R[i,m]
                D_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
                2 * size(I)[1] +
                2 * sum(size(M[i])[1] for i in I) +
                length(TR) +
                2 * length(SQ) +
                2 * size(J)[1] + 1, count_imr] = δ[i,m,r]
                count_imr += 1
            end
        end
    end
end

# 8
A[1, 1, 1][1:tfin, 2:ST] = [(t==tt) ? -1 : 0 for t in 1:tfin, tt in 1:tfin]
D[1, 1, 1][1:tfin, 1:tfin] = [(t==tt) ? ILmin[t] : 0 for t in 1:tfin, tt in 1:tfin]

A[1, 1, 1][tfin + 1:2*tfin, 2:ST] = [(t==tt) ? 1 : 0 for t in 1:tfin, tt in 1:tfin]
D[1, 1, 1][tfin + 1:2*tfin, 1:tfin] = [(t==tt) ? -ILmax[t] : 0 for t in 1:tfin, tt in 1:tfin]

# 9b-1
count_cons = 1
for i in I
    for m in M[i]
        A[1,1,1][2*tfin + count_cons, 1+ tfin + count_cons] = 1
        b[1,1,1][2*tfin + count_cons] = yini[i,m]
        global count_cons += 1
    end
end

# 9b-2
count_cons = 1
for i in I
    for m in M[i]
        A[1,1,1][2*tfin + sum(size(M[i])[1] for i in I) + count_cons, 1+ tfin + count_cons] = -1
        b[1,1,1][2*tfin + sum(size(M[i])[1] for i in I) + count_cons] = -yini[i,m]
        global count_cons += 1
    end
end

# 9c-1
count_cons = 1
for t in -θmax+2:0
    for ((i,j),(m,m_)) in TR
        D[1,1,1][2 * tfin +
        2 * sum(size(M[i])[1] for i in I) + count_cons, tfin + (t-2+θmax)*(sum(size(M[i])[1]^2 for i in I)) + var_z[i] + (M_trans[i,m]-1)*size(M[i])[1]+M_trans[i,m_]]= 1
        b[1,1,1][2 * tfin +
        2 * sum(size(M[i])[1] for i in I) + count_cons] = zini[i,m,m_,t]
        global count_cons += 1
    end
end

# 9c-2
count_cons = 1
count_cons = 1
for t in -θmax+2:0
    for ((i,j),(m,m_)) in TR
        D[1,1,1][2 * tfin +
        2 * sum(size(M[i])[1] for i in I) + length(TR) * (θmax - 1) + count_cons, tfin + (t-2+θmax)*(sum(size(M[i])[1]^2 for i in I)) + var_z[i] + (M_trans[i,m]-1)*size(M[i])[1]+M_trans[i,m_]]= -1
        b[1,1,1][2 * tfin +
        2 * sum(size(M[i])[1] for i in I) + length(TR) * (θmax - 1) + count_cons] = -zini[i,m,m_,t]
        global count_cons += 1
    end
end

# 16a
##### for modes: :MPL
for t in 2:ST
    local count_cons = 1
    local count_var = 1
    for i in I
        for m in M[i]
            for r in R[i,m]
                for f in F[i,m,r]
                    D[t,1,1][2 * sum(size(M[i])[1] for i in I) +
                    2 * size(I)[1] +
                    2 * sum(size(M[i])[1] for i in I) +
                    length(TR) +
                    2 * length(SQ) +
                    2 * size(J)[1] +
                    1 + count_cons, tfin +
                    (θmax - 1) * (sum(size(M[i])[1]^2 for i in I)) + (t-2)*(sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + count_var] = sum(a[i, m, r, f, j] * n[i, m, r, f, j] for j in J)
                    D_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
                    2 * size(I)[1] +
                    2 * sum(size(M[i])[1] for i in I) +
                    length(TR) +
                    2 * length(SQ) +
                    2 * size(J)[1] +
                    1 + count_cons, count_var] = sum(a[i, m, r, f, j] * n[i, m, r, f, j] for j in J)
                    count_cons += 1
                end
                count_var += 1
            end
        end
    end
end

for t in 2:ST
    local count_j = 1
    for j in J
        local count_cons = 1
        local count_var = 1
        for i in I
            for m in M[i]
                for r in R[i,m]
                    for f in F[i,m,r]
                        A[t,1,1][2 * sum(size(M[i])[1] for i in I) +
                        2 * size(I)[1] +
                        2 * sum(size(M[i])[1] for i in I) +
                        length(TR) +
                        2 * length(SQ) +
                        2 * size(J)[1] +
                        1 + count_cons, 1 +
                        tfin +
                        sum(size(M[i])[1] for i in I) + (t-2)*size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count_var] = -n[i,m,r,f,j]
                        A_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
                        2 * size(I)[1] +
                        2 * sum(size(M[i])[1] for i in I) +
                        length(TR) +
                        2 * length(SQ) +
                        2 * size(J)[1] +
                        1 + count_cons, (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count_var] = -n[i,m,r,f,j]
                        count_cons += 1
                    end
                    count_var += 1
                end
            end
        end
        count_j += 1
    end
end


# 16b-1
for t in 2:ST
    local count_j = 1
    for j in J
        local count_var = 1
        for i in I
            for m in M[i]
                for r in R[i,m]
                    A[t,1,1][2 * sum(size(M[i])[1] for i in I) +
                    2 * size(I)[1] +
                    2 * sum(size(M[i])[1] for i in I) +
                    length(TR) +
                    2 * length(SQ) +
                    2 * size(J)[1] +
                    1 + sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) + (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count_var, 1 +
                    tfin +
                    sum(size(M[i])[1] for i in I) + (t-2)*size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count_var] = 1
                    A_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
                    2 * size(I)[1] +
                    2 * sum(size(M[i])[1] for i in I) +
                    length(TR) +
                    2 * length(SQ) +
                    2 * size(J)[1] +
                    1 + sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) + (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count_var, (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count_var] = 1
                    D[t,1,1][2 * sum(size(M[i])[1] for i in I) +
                    2 * size(I)[1] +
                    2 * sum(size(M[i])[1] for i in I) +
                    length(TR) +
                    2 * length(SQ) +
                    2 * size(J)[1] +
                    1 + sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) + (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count_var, tfin +
                    (θmax - 1) * (sum(size(M[i])[1]^2 for i in I)) + (t-2)*(sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + count_var] = -PDmax[i,m,r,j]
                    D_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
                    2 * size(I)[1] +
                    2 * sum(size(M[i])[1] for i in I) +
                    length(TR) +
                    2 * length(SQ) +
                    2 * size(J)[1] +
                    1 + sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) + (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count_var, count_var] = -PDmax[i,m,r,j]
                    count_var += 1
                end
            end
        end
        count_j += 1
    end
end


# 16b-2
for t in 2:ST
    local count_j = 1
    for j in J
        local count_var = 1
        for i in I
            for m in M[i]
                for r in R[i,m]
                    A[t,1,1][2 * sum(size(M[i])[1] for i in I) +
                    2 * size(I)[1] +
                    2 * sum(size(M[i])[1] for i in I) +
                    length(TR) +
                    2 * length(SQ) +
                    2 * size(J)[1] +
                    1 + sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) + size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count_var, 1 +
                    tfin +
                    sum(size(M[i])[1] for i in I) + (t-2)*size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count_var] = -1
                    A_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
                    2 * size(I)[1] +
                    2 * sum(size(M[i])[1] for i in I) +
                    length(TR) +
                    2 * length(SQ) +
                    2 * size(J)[1] +
                    1 + sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) + size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count_var, (count_j-1)*sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count_var] = -1
                    count_var += 1
                end
            end
        end
        count_j += 1
    end
end


# 17f - PD_hat >= 0
count = 1
for t in 1:tfin
    for j in J
        for i in I
            for m in M[i]
                for r in R[i,m]
                    A[1,1,1][2 * tfin +
                    2 * sum(size(M[i])[1] for i in I) +
                    2 * length(TR) * (θmax - 1) + count, 1 +
                    tfin +
                    sum(size(M[i])[1] for i in I) + count] = -1
                    global count += 1
                end
            end
        end
    end
end


# y_mt >= 0
for t in 2:ST
    A_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
    2 * size(I)[1] +
    2 * sum(size(M[i])[1] for i in I) +
    length(TR) +
    2 * length(SQ) +
    2 * size(J)[1] +
    1 +
    (sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) +
        2 * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)
    )+1: 2 * sum(size(M[i])[1] for i in I) +
    2 * size(I)[1] +
    2 * sum(size(M[i])[1] for i in I) +
    length(TR) +
    2 * length(SQ) +
    2 * size(J)[1] +
    1 +
    (sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) +
        2 * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)
    )+ sum(size(M[i])[1] for i in I), 
    size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + 1:size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)  + sum(size(M[i])[1] for i in I)] = [(k==kk) ? -1 : 0 for k in 1:sum(size(M[i])[1] for i in I), kk in 1:sum(size(M[i])[1] for i in I)]
end


# y_mt <= 1
for t in 2:ST
    local count = 1
    for i in I
        for m in M[i]
            A_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
            2 * size(I)[1] +
            2 * sum(size(M[i])[1] for i in I) +
            length(TR) +
            2 * length(SQ) +
            2 * size(J)[1] +
            1 +
            (
                sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) +
                2 * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)
            ) + sum(size(M[i])[1] for i in I) + count, size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) +
            count] = 1
            b[t,1,1][2 * sum(size(M[i])[1] for i in I) +
            2 * size(I)[1] +
            2 * sum(size(M[i])[1] for i in I) +
            length(TR) +
            2 * length(SQ) +
            2 * size(J)[1] +
            1 +
            (
                sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) +
                2 * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)
            ) + sum(size(M[i])[1] for i in I) + count] = 1
            count += 1
        end
    end
end

# need to add constraints on limits for all recourse integer variables
# y_mrt >= 0
for t in 2:ST
    local count = 1
    for i in I
        for m in M[i]
            for r in R[i,m]
                D[t,1,1][2 * sum(size(M[i])[1] for i in I) +
                2 * size(I)[1] +
                2 * sum(size(M[i])[1] for i in I) +
                length(TR) +
                2 * length(SQ) +
                2 * size(J)[1] +
                1 +
                (
                    sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) +
                    2 * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)
                ) +
                2 * sum(size(M[i])[1] for i in I) + count, tfin +
                (θmax - 1) * (sum(size(M[i])[1]^2 for i in I)) + (t-2)*(sum(sum(length(R[i, m]) for m in M[i]) for i in I))+count] = -1
                D_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
                2 * size(I)[1] +
                2 * sum(size(M[i])[1] for i in I) +
                length(TR) +
                2 * length(SQ) +
                2 * size(J)[1] +
                1 +
                (
                    sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) +
                    2 * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)
                ) +
                2 * sum(size(M[i])[1] for i in I) + count, count] = -1
                count += 1
            end
        end
    end
end

# y_mrt <= 1
for t in 2:ST
    local count = 1
    for i in I
        for m in M[i]
            for r in R[i,m]
                D[t,1,1][2 * sum(size(M[i])[1] for i in I) +
                2 * size(I)[1] +
                2 * sum(size(M[i])[1] for i in I) +
                length(TR) +
                2 * length(SQ) +
                2 * size(J)[1] +
                1 +
                (
                    sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) +
                    2 * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)
                ) +
                2 * sum(size(M[i])[1] for i in I) + sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count, tfin +
                (θmax - 1) * (sum(size(M[i])[1]^2 for i in I)) + (t-2)*(sum(sum(length(R[i, m]) for m in M[i]) for i in I))+count] = 1
                D_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
                2 * size(I)[1] +
                2 * sum(size(M[i])[1] for i in I) +
                length(TR) +
                2 * length(SQ) +
                2 * size(J)[1] +
                1 +
                (
                    sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) +
                    2 * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)
                ) +
                2 * sum(size(M[i])[1] for i in I) + sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count, count] = 1
                b[t,1,1][2 * sum(size(M[i])[1] for i in I) +
                2 * size(I)[1] +
                2 * sum(size(M[i])[1] for i in I) +
                length(TR) +
                2 * length(SQ) +
                2 * size(J)[1] +
                1 +
                (
                    sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) +
                    2 * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)
                ) +
                2 * sum(size(M[i])[1] for i in I) + sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count] = 1
                count += 1
            end
        end
    end
end

# z_mm't >= 0
for t in 2:ST
    local count = 1
    for i in I
        for m_ in M[i]
            for m__ in M[i]
                D_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
                2 * size(I)[1] +
                2 * sum(size(M[i])[1] for i in I) +
                length(TR) +
                2 * length(SQ) +
                2 * size(J)[1] +
                1 +
                (
                    sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) +
                    2 * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)
                ) +
                2 * sum(size(M[i])[1] for i in I) +
                2 * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count, (sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + count] = -1
                count += 1
            end
        end
    end
end

# z_mm't <= 1
for t in 2:ST
    local count = 1
    for i in I
        for m_ in M[i]
            for m__ in M[i]
                D_tilde[t,t][2 * sum(size(M[i])[1] for i in I) +
                2 * size(I)[1] +
                2 * sum(size(M[i])[1] for i in I) +
                length(TR) +
                2 * length(SQ) +
                2 * size(J)[1] +
                1 +
                (
                    sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) +
                    2 * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)
                ) +
                2 * sum(size(M[i])[1] for i in I) +
                2 * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + sum(size(M[i])[1]^2 for i in I) + count, (sum(sum(length(R[i, m]) for m in M[i]) for i in I)) + count] = 1
                b[t,1,1][2 * sum(size(M[i])[1] for i in I) +
                2 * size(I)[1] +
                2 * sum(size(M[i])[1] for i in I) +
                length(TR) +
                2 * length(SQ) +
                2 * size(J)[1] +
                1 +
                (
                    sum(sum(sum(length(F[i, m, r]) for r in R[i, m]) for m in M[i]) for i in I) +
                    2 * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)
                ) +
                2 * sum(size(M[i])[1] for i in I) +
                2 * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + sum(size(M[i])[1]^2 for i in I) + count] = 1
                count += 1
            end
        end
    end
end

penalty = Dict((st) => 0.01*ones(1, size(J)[1] * sum(sum(length(R[i,m]) for m in M[i]) for i in I)) for st = 2:ST)