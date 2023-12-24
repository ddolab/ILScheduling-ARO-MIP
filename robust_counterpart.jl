model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "MIPGap", 0.005)

# first stage variables
@variable(model, x1[1:P[1]])
@variable(model, y1[1:Q[1]], Bin)

# Quick way to find the solution when no interruptible load is provided:
# @constraint(model, x1[2:ST] .== 0)

# variables from decision rules
X_bar = Dict((stt, sttt, i) => @variable(model, [1:P[stt], 1:r[sttt, i]]) for stt in 2:ST for sttt in zetaset[stt] for i in 1:K[sttt])
X_hat = Dict((stt, sttt, i) => @variable(model, [1:P[stt], 1:g[sttt, i]]) for stt in 2:ST for sttt in zetaset[stt] for i in 1:K[sttt])
Y_hat = Dict((stt, sttt, i) => @variable(model, [1:Q[stt], 1:g[sttt, i]], Int) for stt in 2:ST for sttt in zetaset[stt] for i in 1:K[sttt])

f = Dict((st, sttt, i) => @variable(model, [1:N[st]]) for st in 2:ST for sttt in 1:st for i in 1:K[sttt])

# dual variables
Delta = Dict((st, sttt, i) => @variable(model, [1:N[st]]) for st in 2:ST for sttt in 1:st for i in 1:K[sttt])
Phi = Dict(st => @variable(model, [1:N[st], 1:Mu[st]], lower_bound = 0.0) for st in 2:ST)

# variables required to introduce penalty on recourse coefficients
tx_bar = Dict((stt, sttt, i) => @variable(model, [1:size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I), 1:r[sttt, i]]) for stt in 2:ST for sttt in zetaset[stt] for i in 1:K[sttt])
tx_hat = Dict((stt, sttt, i) => @variable(model, [1:size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I), 1:g[sttt, i]]) for stt in 2:ST for sttt in zetaset[stt] for i in 1:K[sttt])

@constraint(model, [stt in 2:ST, sttt in zetaset[stt], i in K[sttt]], tx_bar[stt, sttt, i] .>= X_bar[stt, sttt, i][1:size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)])
@constraint(model, [stt in 2:ST, sttt in zetaset[stt], i in K[sttt]], tx_bar[stt, sttt, i] .>= -X_bar[stt, sttt, i][1:size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)])
@constraint(model, [stt in 2:ST, sttt in zetaset[stt], i in K[sttt]], tx_hat[stt, sttt, i] .>= X_hat[stt, sttt, i][1:size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)])
@constraint(model, [stt in 2:ST, sttt in zetaset[stt], i in K[sttt]], tx_hat[stt, sttt, i] .>= -X_hat[stt, sttt, i][1:size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I)])

# limits on Y_hat as Y_hat \in {0, ±1, ±2}
@constraint(model, [st in 2:ST, stt in zetaset[st], i in 1:K[stt]], Y_hat[st, stt, i] .<= 2)
@constraint(model, [st in 2:ST, stt in zetaset[st], i in 1:K[stt]], Y_hat[st, stt, i] .>= -2)


@constraint(model, A[1, 1, 1]*x1 + D[1, 1, 1]*y1 .<= b[1, 1, 1])

# 30b
@constraint(model, [t in 2:ST, st in 1:t, i in 1:K[st]], f[t, st, i] .== A[t, st, i]*x1 + D[t, st, i]*y1 - b[t, st, i])

# 28a
@constraint(model, [t in 2:ST], Phi[t]*u[t] + sum(sum(Delta[t, st, i] for i in 1:K[st]) for st in 1:t) .<= 0)

# 28b
@constraint(model, [st in 2:ST, sttt in 1:st, i in 1:K[sttt], j in 1:r[sttt, i], k in 1:min(2, sttt)],
(f[st, sttt, i] - Phi[st] * W[st, sttt, i]) * sum(v_bar[sttt, i, j][k]) - Delta[st, sttt, i]
+ sum((A_tilde[st, stt] * X_bar[stt, sttt, i] * v_bar[sttt, i, j][k] + A_tilde[st, stt] * X_hat[stt, sttt, i] * v_hat[sttt, i, j][k]
+ D_tilde[st, stt] * Y_hat[stt, sttt, i] * v_hat[sttt, i, j][k]) for stt in inverse_zeta(sttt, st, zetaset)) .<= 0)

## PDhat and PChat have a constant term (X_bar[st, 1, 1] + X_hat[st, 1, 1]) which adds more flexibility to the value of PDtilde and PCtilde
# these constraints ensure that this flexibility is removed but need to add this restriction only for PDtilde and PCtilde and not y_mt
@constraint(model, [st in 2:ST], X_bar[st, 1, 1][1:size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I), :] .== 0)
@constraint(model, [st in 2:ST], X_hat[st, 1, 1][1:size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I), :] .== 0)

if recourse_type == 1
    @constraint(model, [st in 2:ST], Y_hat[st, 1, 1][1:(sum(sum(length(R[i, m]) for m in M[i]) for i in I)), :] .== 0)
end

for st in 2:ST
    for stt in zetaset[st]
        for i in K[stt]
            for p in 1:P[st]
                for j in 1:r[stt, i]
                    set_start_value(X_bar[st, stt, i][p, j], 0)
                end
                for j in 1:g[stt, i]
                    set_start_value(X_hat[st, stt, i][p, j], 0)
                end
            end
            for q in 1:Q[st]
                for j in 1:g[stt, i]
                    set_start_value(Y_hat[st, stt, i][q, j], 0)
                end
            end
            for n in 1:N[st]
                set_start_value(Delta[st, stt, i][n], 0)
            end
        end
    end
    for n in 1:N[st]
        for m in Mu[st]
            set_start_value(Phi[st][n,m], 0)
        end
    end
end


@objective(model, Min, x1[1] + sum(sum(sum(sum(penalty[st]*tx_bar[st, sttt, i])/r[sttt,i] for i in K[sttt]) for sttt in zetaset[st]) for st in 2:ST) + sum(sum(sum(sum(penalty[st]*tx_hat[st, sttt, i])/g[sttt,i] for i in K[sttt]) for sttt in zetaset[st]) for st in 2:ST))

optimize!(model)