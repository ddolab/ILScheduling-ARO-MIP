function vertices_generation(breakpoints, T, K, R)
    vertices = Dict((t, i, j) => [zeros(R[t, i]), zeros(R[t, i])] for t in 2:T for i in 1:K[t] for j in 1:R[t, i])
    vertices[1, 1, 1] = [[1], [1]]
    for t in 2:T
        for i in 1:K[t]
            for j in 1:R[t, i]
                left = zeros(R[t, i])
                right = zeros(R[t, i])
                if j == 1
                    left[1] = breakpoints[t, i][1]
                    right[1] = breakpoints[t, i][2]
                elseif j == 2
                    left[1] = breakpoints[t, i][2]
                    right[1] = breakpoints[t, i][2]
                    right[2] = breakpoints[t, i][3] - breakpoints[t, i][2]
                else
                    left[1] = breakpoints[t, i][2]
                    right[1] = breakpoints[t, i][2]
                    right[2] = breakpoints[t, i][3] - breakpoints[t, i][2]
                    for jj in 3:j
                        left[jj-1] = breakpoints[t, i][jj] - breakpoints[t, i][jj-1]
                        right[jj] = breakpoints[t, i][jj+1] - breakpoints[t, i][jj]
                    end
                end
                vertices[(t, i, j)][1] = left
                vertices[(t, i, j)][2] = right
            end
        end
    end
    return vertices
end

function G_generation(T, K, R)
    G = Dict((t, i, j) => [zeros(R[t, i]-1), zeros(R[t, i]-1)] for t in 2:T for i in 1:K[t] for j in 1:R[t, i])
    G[(1, 1, 1)] = [[1], [1]]
    for t in 2:T
        for i in 1:K[t]
            if R[t, i] == 1
                G[(t, i, 1)] = [[1], [1]]
            else
                for j in 1:R[t, i]
                    left = zeros(R[t, i] - 1)
                    right = zeros(R[t, i] - 1)
                    if j >= 2
                        for jj in 1:(j-1)
                            left[jj] = 1
                            right[jj] = 1
                        end
                    end
                    G[(t, i, j)][1] = left
                    G[(t, i, j)][2] = right
                end
            end
        end
    end
    return G
end


# need a function that creates zeta sets
function zeta_sets(stages, zeta_)
    # stages = [1,2,....ST]
    Zeta_Sets = Dict()
    for t in stages
        if t == 1
            Zeta_Sets[t] = 1:1
        # elseif t-1 - zeta[t-1] == 1
        #     Zeta_Sets[t] = 1:t
        # elseif t-1 - zeta[t-1] > 1
        #     Zeta_Sets[t] = append!([1], t-zeta[t-1]:t)
        else Zeta_Sets[t] = append!([1], max(t-zeta_, 2):t)
        end
    end
    return Zeta_Sets
end

# needed in the robust counterpart
function inverse_zeta(sttt, st, zetaset)
    inv_zeta = []
    for t in 2:st
        for stt in zetaset[t]
            if sttt == stt
                append!(inv_zeta, t)
            end
        end
    end
    return inv_zeta
end


# to find the big-M parameter Omega
function dECmax(i,m,r,J,F,a,n,gamma)
    M = Model(Gurobi.Optimizer)
    set_optimizer_attribute(M, MOI.Silent(), true)

    @variable(M, PD1[j in J])
    @variable(M, PD2[j in J])
    
    for f in F[i,m,r]
        @constraint(M, sum((PD1[j]-a[i,m,r,f,j])*n[i,m,r,f,j] for j in J) >= 0)
        @constraint(M, sum((PD2[j]-a[i,m,r,f,j])*n[i,m,r,f,j] for j in J) >= 0)
    end
    
    @objective(M, Max,sum(gamma[i,m,r,j]*(PD1[j]-PD2[j]) for j in J))
    optimize!(M)
    
    return objective_value(M)
end

