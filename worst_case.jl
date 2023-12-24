PDbar = Dict((t, st, j, i, m, r) => zeros(2) for t in 2:ST for st in 2:t for j in J for i in I for m in M[i] for r in R[i, m])
PDhat = Dict((t, st, j, i, m, r) => zeros(g[t, 1]) for t in 2:ST for st in 2:t for j in J for i in I for m in M[i] for r in R[i, m])
yhat = Dict((t, st, i, m, r) => zeros(g[t, 1]) for t in 2:ST for st in 2:t for i in I for m in M[i] for r in R[i, m])
PDnom = Dict((t, j, i, m, r) => 0.0 for t in 2:ST for j in J for i in I for m in M[i] for r in R[i, m])


for t in 2:ST
    local count = 1
    for j in J
        for i in I
            for m in M[i]
                for r in R[i, m]
                    if recourse_type == 1
                        PDnom[t, j, i, m, r] = value.(x1)[1 + tfin + sum(size(M[i])[1] for i in I) + (t-2)*size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count]
                    else
                        PDnom[t, j, i, m, r] = value.(x1)[1 + tfin + (1+tfin)*sum(size(M[i])[1] for i in I) + (t-2)*size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + count]
                    end
                    count += 1
                end
            end
        end
    end
end

for t in 2:ST
    for st in zetaset[t]
        local count = 1
        for j in J
            for i in I
                for m in M[i]
                    for r in R[i, m]
                        PDbar[t, st, j, i, m, r] = value.(X_bar[t, st, 1])[count, :]
                        PDhat[t, st, j, i, m, r] = value.(X_hat[t, st, 1])[count, :]
                        count += 1
                    end
                end
            end
        end
    end
end

if recourse_type == 1
    for t in 2:ST
        for st in zetaset[t]
            local count = 1
            for i in I
                for m in M[i]
                    for r in R[i, m]
                        yhat[t, st, i, m, r] = value.(Y_hat[t, st, 1])[count, :]
                        count += 1
                    end
                end
            end
        end
    end
end

using JuMP, Gurobi

wc = Model(Gurobi.Optimizer)
@variable(wc, λ[t in 2:ST, i in K[t], j in 1:r[t, i], k in 1:2] >= 0)

@constraint(wc, [t in 2:ST, i in 1:K[t]], sum(sum(λ[t, i, j, k] for k = 1:2) for j in 1:r[t,i]) == 1)
@constraint(wc, [t in 2:ST, i in 1:K[t], j in 1:r[t,i], k in 1:2], λ[t, i, j, k] <= 1)

@expressions(wc, begin
    w[t in 2:ST, i in 1:K[t]], sum(sum(λ[t, i, j, k]*breakpoints[t, i][j-1+k] for k in 1:2) for j in 1:r[t, i])
end)

@constraint(wc, sum(w[t,1] for t in 2:ST)>=load_reduction_requests)

@constraint(wc, [t in 2:ST, i in 1:K[t]], w[t, i] <= 1)
@constraint(wc, [t in 2:ST, i in 1:K[t]], sum(sum(w[k, i] for i in 1:K[k]) for k = 2:t) <= Γ[t-1])

@expressions(wc, begin
    w_bar[t in 2:ST, i in 1:K[t]], sum(sum(λ[t, i, j, k]*v_bar[t, i, j][k] for k in 1:2) for j in 1:r[t,i])
end)
@expressions(wc, begin
    w_hat[t in 2:ST, i in 1:K[t]], sum(sum(λ[t, i, j, k]*v_hat[t, i, j][k] for k in 1:2) for j in 1:r[t,i])
end)

#PD_imrjt
@expressions(wc, begin
    PDimrjt[t in 2:ST, j in J, i in I, m in M[i], r in R[i, m]], sum(sum(PDbar[t, st, j, i, m, r][k]*w_bar[st, 1][k] for k in 1:2) + sum(PDhat[t, st, j, i, m, r][k]*w_hat[st, 1][k] for k in 1:g[t, 1]) for st in max(2, t-zeta_):t)
end)


#y_imrt
if recourse_type == 1
    @expressions(wc, begin
        yimrt[t in 2:ST, i in I, m in M[i], r in R[i, m]], sum(sum(yhat[t, st, i, m, r][k]*w_hat[st, 1][k] for k in 1:g[t, 1]) for st in 2:t)
    end)
end

@constraint(wc, [t in 2:ST], value.(y1[t-1]) -w[t, 1] >= 0)

if recourse_type == 1
    @objective(wc, Max, sum(alphaEC[t-1]*sum(sum(sum(sum(γ[i, m, r, j]*PDimrjt[t, j, i, m, r] for r in R[i, m]) for m in M[i]) for i in I) for j in J) + alphaEC[t-1]*sum(sum(sum(δ[i, m, r]*yimrt[t, i, m, r] for r in R[i, m]) for m in M[i]) for i in I) for t in 2:ST))
else
    @objective(wc, Max, sum(alphaEC[t-1]*sum(sum(sum(sum(γ[i, m, r, j]*PDimrjt[t, j, i, m, r] for r in R[i, m]) for m in M[i]) for i in I) for j in J) for t in 2:ST))
end
optimize!(wc)

#  Product flows
PD_tilde = Dict((t, j, i, m, r) => 0.0 for t in 1:tfin for j in J for i in I for m in M[i] for r in R[i, m])
y = Dict((t, i, m, r) => 0.0 for t in 1:tfin for i in I for m in M[i] for r in R[i, m])
y_tilde = Dict((t, i, m, r) => 0.0 for t in 1:tfin for i in I for m in M[i] for r in R[i, m])

for t in 1:tfin
    local count = 1
    for j in J
        for i in I
            for m in M[i]
                for r in R[i, m]
                    PD_tilde[t, j, i, m, r] = sum(sum(PDbar[t+1, st, j, i, m, r]'*value.(w_bar[st, k]) + PDhat[t+1, st, j, i, m, r]'*value.(w_hat[st, k]) for k in 1:K[st]) for st in max(2, t+1-zeta_):(t+1))
                end
            end
        end
    end
end

count = 1
for t in 1:tfin
    for i in I
        for m in M[i]
            for r in R[i, m]
                y[t, i, m, r] = value.(y1)[tfin + count]
                y_tilde[t, i, m, r] = sum(sum(yhat[t+1, st, i, m, r]'*value.(w_hat[st, k]) for k in 1:K[st]) for st in max(2, t-zeta_):(t+1))
                global count += 1
            end
        end
    end
end


# create vectors of electricity consumed by each process and recourse in electricity consumption
EC = Dict(i => zeros(tfin) for i in I)
EC_tilde = Dict(i => zeros(tfin) for i in I)
IL = Vector{Float64}(undef, tfin)

for i in I
    for t in 1:tfin
        EC[i][t] = sum(sum(δ[i, m, r]*y[t, i, m, r] + sum(γ[i, m, r, j]*PDnom[t+1, j, i, m, r] for j in J) for r in R[i, m]) for m in M[i])
        EC_tilde[i][t] = sum(sum(δ[i, m, r]*y_tilde[t, i, m, r] + sum(γ[i, m, r, j]*PD_tilde[t, j, i, m, r] for j in J) for r in R[i, m]) for m in M[i])
    end
end

EC_nom = sum(EC[i][:] for i in I)

for t in 1:tfin
    IL[t] = value.(x1)[t+1]
end


# Inventory levels################################################
IV = Dict(j => zeros(tfin) for j in J)
countj = 1
for j in J
    IV[j][1] = IVini[j] + sum(sum(sum(PDnom[2,j,i,m,r] for r in R[i,m]) for m in M[i]) for i in Ihat[j]; init = 0) - sum(sum(sum(PDnom[2,j,i,m,r] for r in R[i,m]) for m in M[i]) for i in Ibar[j]; init = 0) - Demand[j,1]
    # IV[j][1] =  value.(x1)[1 +
    # tfin +
    # (1+tfin)*sum(size(M[i])[1] for i in I) +
    # tfin * size(J)[1] * sum(sum(length(R[i, m]) for m in M[i]) for i in I) + countj] + sum(sum(sum(PDnom[2,j,i,m,r] for r in R[i,m]) for m in M[i]) for i in Ihat[j]; init = 0) - sum(sum(sum(PDnom[2,j,i,m,r] for r in R[i,m]) for m in M[i]) for i in Ibar[j]; init = 0) - Demand[j,1]
    global countj += 1
    for t in 2:tfin
        IV[j][t] = IV[j][t-1] + sum(sum(sum(PDnom[t+1,j,i,m,r] for r in R[i,m]) for m in M[i]) for i in Ihat[j]; init = 0) - sum(sum(sum(PDnom[t+1,j,i,m,r] for r in R[i,m]) for m in M[i]) for i in Ibar[j]; init = 0) - Demand[j,t]
    end
end

IVtilde = Dict(j => zeros(tfin) for j in J)
for j in J
    IVtilde[j][1] = sum(sum(sum(PD_tilde[1,j,i,m,r] for r in R[i,m]) for m in M[i]) for i in Ihat[j]; init = 0) - sum(sum(sum(PD_tilde[1,j,i,m,r] for r in R[i,m]) for m in M[i]) for i in Ibar[j]; init = 0)
    for t in 2:tfin
        IVtilde[j][t] = IVtilde[j][t-1] + sum(sum(sum(PD_tilde[t,j,i,m,r] for r in R[i,m]) for m in M[i]) for i in Ihat[j]; init = 0) - sum(sum(sum(PD_tilde[t,j,i,m,r] for r in R[i,m]) for m in M[i]) for i in Ibar[j]; init = 0)
    end
end


for j in J
    IV[j][:] = round.(IV[j][:], digits = 1)
    IVtilde[j][:] = round.(IVtilde[j][:], digits = 1)
end

for i in I
    EC[i][:] = round.(EC[i][:], digits = 1)
    EC_tilde[i][:] = round.(EC_tilde[i][:], digits = 1)
end

EC_nom = round.(EC_nom, digits = 1)
IL = round.(IL, digits = 1)
