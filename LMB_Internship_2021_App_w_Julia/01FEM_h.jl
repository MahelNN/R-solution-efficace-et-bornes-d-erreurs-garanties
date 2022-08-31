######################################################################
## @Autor: Mahel Ndoumba N.                                         ##
## LMB internship 2021 : Efficient resolution and guaranteed bounds ##
##                         for eigenvalue problems in electronic    ##
##                                structure calculations            ##
######################################################################
using Arpack
include("00Exact_sol_h.jl")
#
##################### Subroutines for the finite element mehtod ################
#On a sub interval I of [-R,R], return one of the r^th piece from the P^1 hat fonction
function shape_func(x, r, j, h, R)
    X_v = -R:h:R
    ϕ = zeros(2)

    ϕ[1] = (X_v[j+1] - x) / h
    ϕ[2] = (x - X_v[j]) / h
    # local shape function evaluated in x
    return ϕ[r]
end


# First derivative of the shape function, return the r^th piece
function deriv_shape_func(x, r, j, h, R)
    X_v = -R:h:R
    β = zeros(2)

    β[1] = (X_v[j+1] - x) / h
    β[2] = (x - X_v[j]) / h

    # local derivate of the shape function
    dϕ = zeros(2)

    dϕ[1] = -1 / h
    dϕ[2] = 1 / h
    return dϕ[r]
end


# Mass and Rigidity matrices
function Mat_dis_MK(N, h)
    M = zeros(N, N)
    a = h / 6

    [M[i, i] = 4 for i = 2:N-1]
    [M[i, i+1] = 1 for i = 1:N-1]
    [M[i, i-1] = 1 for i = 2:N]

    M[1, 1] = 2
    M[N, N] = 2

    K = zeros(N, N)
    b = 1 / h

    [K[i, i] = 2 for i = 2:N-1]
    [K[i, i+1] = -1 for i = 1:N-1]
    [K[i, i-1] = -1 for i = 2:N]

    K[1, 1] = 1
    K[N, N] = 1
    return [sparse(a .* M), sparse(b .* K)]
end


# Return the matrix at a cusp position
function Mat_dis_MV(x, N, h, R)
    X_v = -R:h:R
    MV = zeros(N, N)

    for j = 1:N-1
        δ(x) = (X_v[j] <= x) * (x < X_v[j+1])
        ϕ_r(r) = shape_func(x, r, j, h, R) * δ(x)

        MV[j, j+1] = ϕ_r(1) * ϕ_r(2)
        MV[j+1, j] = ϕ_r(1) * ϕ_r(2)

        MV[j, j] += ϕ_r(1)^2
        MV[j+1, j+1] += ϕ_r(2)^2
    end
    return sparse(MV)
end

# Return the L2 norm
function Norm_L2_sol_num(U, Q)
    norm_L2 = U' * Q * U
    return √(norm_L2)
end


# Finite element solver
function FEM_1D1(N, R, z_1, z_2, μ_1, μ_2)
    h = 2 * R / N
    X_v = -R:h:R
    d_glo = N + 1
    # Mass and regidity matrices
    M_glo, K_glo = Mat_dis_MK(d_glo, h)
    # Matrix of the potential
    MV_glo =
        z_1 .* Mat_dis_MV(μ_1, d_glo, h, R) +
        z_2 .* Mat_dis_MV(μ_2, d_glo, h, R)

    A_glo = K_glo .- MV_glo

    A_h = A_glo[2:(d_glo-1), 2:(d_glo-1)]
    M_h = M_glo[2:(d_glo-1), 2:(d_glo-1)]

    # A_h and M_h were defined as sparse matrices
    # so we transform them as regular/plain matrices
    A_h = Array(A_h)
    M_h = Array(M_h)

    # eigen solver
    λ_ef, U_ef = eigen(A_h, M_h)

    U_dis = [0.0; U_ef[:, 1]; 0.0]
    # Take the positive components of the numerical solution
    U_dis = Pos_vec(U_dis)
    return [λ_ef[1], U_dis]
end


# L2 and H1 error norms
function Err_L2_H1(U_dis, h, n, R, λ, z_1, z_2, μ_1, μ_2)
    X_v = -R:h:R
    d_loc = 2
    l2_err = 0
    h1_err = 0

    for j = 1:n  #% Loop on the cells
        x1 = (X_v[j+1] + X_v[j]) / 2 - sqrt(3) / (2 * sqrt(5)) * h
        x2 = (X_v[j+1] + X_v[j]) / 2
        x3 = (X_v[j+1] + X_v[j]) / 2 + sqrt(3) / (2 * sqrt(5)) * h

        u1 = Exact_sol(x1, λ, z_1, z_2, μ_1, μ_2)
        u2 = Exact_sol(x2, λ, z_1, z_2, μ_1, μ_2)
        u3 = Exact_sol(x3, λ, z_1, z_2, μ_1, μ_2)

        du1 = dExact_sol(x1, λ, z_1, z_2, μ_1, μ_2)
        du2 = dExact_sol(x2, λ, z_1, z_2, μ_1, μ_2)
        du3 = dExact_sol(x3, λ, z_1, z_2, μ_1, μ_2)

        U_loc = zeros(d_loc, 1)
        Phi1 = zeros(d_loc, 1)
        Phi2 = zeros(d_loc, 1)
        Phi3 = zeros(d_loc, 1)

        dPhi1 = zeros(d_loc, 1)
        dPhi2 = zeros(d_loc, 1)
        dPhi3 = zeros(d_loc, 1)

        for r_loc = 1:d_loc  #% Loop on the local dof
            r_glo = (j - 1) + r_loc # % Local-to-global map
            U_loc[r_loc, 1] = U_dis[r_glo, 1]

            Phi1[r_loc, 1] = shape_func(x1, r_loc, j, h, R)
            Phi2[r_loc, 1] = shape_func(x2, r_loc, j, h, R)
            Phi3[r_loc, 1] = shape_func(x3, r_loc, j, h, R)

            dPhi1[r_loc, 1] = deriv_shape_func(x1, r_loc, j, h, R)
            dPhi2[r_loc, 1] = deriv_shape_func(x2, r_loc, j, h, R)
            dPhi3[r_loc, 1] = deriv_shape_func(x3, r_loc, j, h, R)
        end

        l2_err =
            l2_err .+ 5 * h / 18 .* (u1 .- U_loc' * Phi1)^2 .+
            4 * h / 9 .* (u2 .- U_loc' * Phi2)^2 .+
            5 * h / 18 .* (u3 .- U_loc' * Phi3)^2

        h1_err =
            h1_err .+ 5 * h / 18 .* (du1 .- U_loc' * dPhi1)^2 .+
            4 * h / 9 .* (du2 .- U_loc' * dPhi2)^2 .+
            5 * h / 18 .* (du3 .- U_loc' * dPhi3)^2
    end
    return [√(l2_err[1]), √(l2_err[1] + h1_err[1])]
end


# Return the relative error due to the trocature of the solution on [-R,R]
function Err_R(R, λ, z_1, z_2, μ_1, μ_2)
    α = √(-λ)
    # Paramètres Aᵢ, Bᵢ for the exact soluiton
    B1, A2, B2, A3 = Coef_exact_sol(λ, z_1, z_2, μ_1, μ_2)

    err_R = (B1 + A3) * exp(-R * α) / α
    return err_R
end


# Interpolate a solution on a finier mesh
function FE_interpolation(U_N, R, scale)
    N = size(U_N)[1]
    h0 = 2 * R / (N - 1) # original mesh size
    X_0 = -R:h0:R

    h = h0 / scale # new mesh size

    d_glo = ceil(Int64, 2 * R / h) + 1 #dof

    X_v = -R:h:R

    U_N1 = zeros(d_glo) # new discretized vector solution

    # Define the new vector components
    for j = 1:N-1
        for k = 1:d_glo
            coef_kp1 = ((k - 1) % scale) / scale
            coef_k = 1 - coef_kp1
            if (X_0[j] <= X_v[k]) * (X_v[k] < X_0[j+1])
                U_N1[k] = (coef_k * U_N[j] + coef_kp1 * U_N[j+1])
            else # copy only the right components
                continue
            end
        end
    end

    ϕ_r(r, l) = shape_func.(X_v, r, l, h, R)

    for k = 1:d_glo-1
        X_j = zeros(d_glo)

        [X_j[l] = 1 for l = k:k+1] # interpolate only on the relevant cells

        U_N1 =
            U_N1 .+ U_N1[k] .* ϕ_r(1, k) .* X_j .+ U_N1[k+1] .* ϕ_r(2, k) .* X_j
    end
    return U_N1
end


# Compute the residual and its dual norm
function Res(λ_N, U_N, R, z_1, z_2, μ_1, μ_2)
    N = size(U_N)[1]
    h0 = 2 * R / (N - 1) # original mesh size
    X_0 = -R:h0:R

    scale = 2
    h = h0 / scale # new nesh size

    d_glo = ceil(Int64, 2 * R / h) + 1

    X_v = -R:h:R

    U_N1 = FE_interpolation(U_N, R, scale)

    M_glo, K_glo = Mat_dis_MK(d_glo, h)
    MV_glo =
        z_1 .* Mat_dis_MV(μ_1, d_glo, h, R) +
        z_2 .* Mat_dis_MV(μ_2, d_glo, h, R)

    U_N1 ./= √(U_N1' * M_glo * U_N1) # L^2 normalization

    A_glo = K_glo - MV_glo # we use the V norm and not the H^1 one

    F_glo = (K_glo - MV_glo - λ_N .* M_glo) * U_N1

    A_h = A_glo[2:d_glo-1, 2:d_glo-1]
    F_h = F_glo[2:d_glo-1]

    # The vector that realise the supremum in the résidual solves
    V_res = \(A_h, F_h)

    V_res = [0.0; V_res; 0.0]
    # Computation of the dual norm of the residual
    dual_norm_res = (V_res' * A_glo * V_res)
    return [√dual_norm_res, V_res]
end


##
EoF = "✔"
