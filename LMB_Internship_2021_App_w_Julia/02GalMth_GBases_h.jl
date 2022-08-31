######################################################################
## @Autor: Mahel Ndoumba N.                                         ##
## LMB internship 2021 : Efficient resolution and guaranteed bounds ##
##                         for eigenvalue problems in electronic    ##
##                                structure calculations            ##
######################################################################
using Optim
include("01FEM_h.jl")

##
# Gaussian functions used as basis functions
function Gaussian_m_σ(x, m, σ)
    C_μ = 1 / (σ * √(2π))
    α_μ = 1 / (2 * σ^2)
    return C_μ * exp(-α_μ * (x - m)^2)
end

function deriv_Gaussian_m_σ(x, m, σ)
    C_μ = 1 / (σ * √(2π))
    α_μ = 1 / (2 * σ^2)
    return -2 * α_μ * (x - m) * C_μ * exp(-α_μ * (x - m)^2)
end


# Parameters used in the formula for Gaussian functions product rule
function Parameters_gauss_prod(m_1, σ_1, m_2, σ_2)
    C_μ = 1 / (σ_1 * √(2π))
    α_μ = 1 / (2 * σ_1^2)
    C_ν = 1 / (σ_2 * √(2π))
    α_ν = 1 / (2 * σ_2^2)

    α_μν = (α_μ * α_ν) / (α_μ + α_ν)
    C_μν = C_μ * C_ν * exp(-α_μν * (m_1 - m_2)^2)
    γ_μν = α_μ + α_ν
    m_μν = (α_μ * m_1 + α_ν * m_2) / (α_μ + α_ν)
    return [C_μν, γ_μν, m_μν, α_μ, α_ν]
end


# Compute the product of two Gaussian function using their Gaussian parameters
function Gaussian_product(x, m_1, σ_1, m_2, σ_2)
    C_μν, γ_μν, m_μν, α_μ, α_ν = Parameters_gauss_prod(m_1, σ_1, m_2, σ_2)
    return C_μν * exp(-γ_μν * (x - m_μν)^2)
end


# Compute the integrale of two Gaussian functions
function Int_gauss_prod(m_1, σ_1, m_2, σ_2)
    C_μν, γ_μν, m_μν, α_μ, α_ν = Parameters_gauss_prod(m_1, σ_1, m_2, σ_2)
    return C_μν * √(π / γ_μν)
end

# Compute the integrale of the derivative of two Gaussian functions
function Int_deriv_gauss_prod(m_1, σ_1, m_2, σ_2)
    C_μν, γ_μν, m_μν, α_μ, α_ν = Parameters_gauss_prod(m_1, σ_1, m_2, σ_2)
    S_μν = 4 * α_μ * α_ν * C_μν * √(π / γ_μν)
    return S_μν * (1 / 2 / γ_μν + (m_μν - m_1) * (m_μν - m_2))
end


# Use the Gaussian functions as basis functions
function Basis_fct(x, Vec_μ, Vec_σ, r)
    N = length(Vec_σ)
    @assert (0 < r) * (r <= N)
    # X_v = -R:h:R
    ϕ_r = Gaussian_m_σ(x, Vec_μ[r], Vec_σ[r])
    return ϕ_r
end


function deriv_Basis_fct(x, Vec_μ, Vec_σ, r)
    N = length(Vec_σ)
    @assert (0 < r) * (r <= N)
    # X_v = -R:h:R
    dϕ_r = deriv_Gaussian_m_σ(x, Vec_μ[r], Vec_σ[r])
    return dϕ_r
end


# Interpolated the numerical solution on the basis set
function Interpolated_num_sol(x, U_i, Vec_μ, Vec_σ)
    N = length(Vec_σ)
    # U(x) = ∑_{i=1:N}Uᵢϕᵢ(x)
    U_sol = U_i[1] .* Basis_fct(x, Vec_μ, Vec_σ, 1)
    for i = 2:N
        ϕ_i = Basis_fct(x, Vec_μ, Vec_σ, i)
        U_sol = U_sol .+ U_i[i] .* ϕ_i
    end
    return U_sol
end


# Galerkin method with gaussian bases solver
function GM_BG_VD1D(Vec_μ, Vec_σ, μ_1, μ_2, z_1, z_2, R, h_sol)
    # Initialization
    N = length(Vec_μ)
    X_dis = -R:h_sol:R

    # Discretization matrices
    B = zeros(N, N)
    A = zeros(N, N)
    MV = zeros(N, N)

    for i = 1:N
        for j = 1:N
            # Gaussian function parameters
            μ_i = Vec_μ[i]
            μ_j = Vec_μ[j]
            σ_i = Vec_σ[i]
            σ_j = Vec_σ[j]

            # B_ij = ∫ϕᵢϕⱼ
            B[i, j] = Int_gauss_prod(μ_i, σ_i, μ_j, σ_j)
            # Potential Vu = z₁(ϕᵢϕⱼ)(μ₁) + z₂(ϕᵢϕⱼ)(μ₂)

            delta(x) = Gaussian_product(x, μ_i, σ_i, μ_j, σ_j)
            MV[i, j] = z_1 * delta(μ_1) + z_2 * delta(μ_2)

            # A_ij = ∫ϕᵢ'ϕⱼ' - z₁(ϕᵢϕⱼ)(μ₁) - z₂(ϕᵢϕⱼ)(μ₂)
            A[i, j] = Int_deriv_gauss_prod(μ_i, σ_i, μ_j, σ_j) - MV[i, j]
        end
    end
    # eigen solver
    λ_h, U = eigen(A, B)
    # interpolate the numerical solution
    U_dis = [Interpolated_num_sol(x_i, U[:, 1], Vec_μ, Vec_σ) for x_i in X_dis]
    # Take the numerical solution with positive components
    U_dis = Pos_vec(real.(U_dis))
    return [λ_h[1], U_dis]
end


# Vector of nucleus positions
function Vec_position(Nμ_1, Nμ_2, μ_1, μ_2)
    Vec_μ_1 = μ_1 .* ones(Nμ_1)
    Vec_μ_2 = μ_2 .* ones(Nμ_2)
    return [Vec_μ_1; Vec_μ_2]
end


# Greedy algorithm to find the optimum vector of standard deviation
# with respect to the lower first approximated eigenvalue
function Opt_vec_σ_eig(Nμ_1, Nμ_2, R, z_1, z_2, μ_1, μ_2, h_sol)
    function f(σ) # objective for the initial optimization
        return GM_BG_VD1D(Vec_μ, [σ; σ], μ_1, μ_2, z_1, z_2, R, h_sol)[1]
    end
    function f2(σ) # objective function
        return GM_BG_VD1D(Vec_μ, [Vec_σ; σ], μ_1, μ_2, z_1, z_2, R, h_sol)[1]
    end
    # Find the first two standard deviation at once
    Vec_μ = Vec_position(1, 1, μ_1, μ_2)

    σ_opt = Optim.minimizer(optimize(f, 0.5, 5))

    Vec_σ = [σ_opt; σ_opt]

    Vec_σ1 = [σ_opt]

    Vec_σ2 = [σ_opt]

    nμ_1 = 1

    nμ_2 = 1

    λ_Nm1, U_Nm1 = GM_BG_VD1D(Vec_μ, Vec_σ, μ_1, μ_2, z_1, z_2, R, h_sol)
    # Alternate the necleus position to find the optimal standard deviation
    # one at a time
    for j = 1:max(Nμ_1, Nμ_2)
        if (nμ_1 < Nμ_1)
            nμ_1 += 1

            Vec_μ = [Vec_μ; μ_1]

            σ_opt = Optim.minimizer(optimize(f2, 1e-5, 7))

            Vec_σ = [Vec_σ; σ_opt]

            Vec_σ1 = [Vec_σ1; σ_opt]
        end

        if (nμ_2 < Nμ_2)
            nμ_2 += 1

            Vec_μ = [Vec_μ; μ_2]

            σ_opt = Optim.minimizer(optimize(f2, 1e-5, 7))

            Vec_σ = [Vec_σ; σ_opt]

            Vec_σ2 = [Vec_σ2; σ_opt]
        end
    end
    return [Vec_σ1; Vec_σ2]
end


# Greedy algorithm to find the optimum vector of standard deviation
# with respect to the dual norm of the residual
function Opt_vec_σ_res(Nμ_1, Nμ_2, R, z_1, z_2, μ_1, μ_2, h_sol)
    function f(σ) # objective for the initial optimization
        λ_N, U_N = GM_BG_VD1D(Vec_μ, [σ; σ], μ_1, μ_2, z_1, z_2, R, h_sol)
        return (Res(λ_N, U_N, R, z_1, z_2, μ_1, μ_2)[1])
    end
    function f2(σ) #objective function
        λ_N, U_N = GM_BG_VD1D(Vec_μ, [Vec_σ; σ], μ_1, μ_2, z_1, z_2, R, h_sol)
        return (Res(λ_N, U_N, R, z_1, z_2, μ_1, μ_2)[1])
    end

    Vec_μ = Vec_position(1, 1, μ_1, μ_2)

    σ_opt = Optim.minimizer(optimize(f, 1, 3))

    Vec_σ = [σ_opt; σ_opt]

    Vec_σ1 = [σ_opt]

    Vec_σ2 = [σ_opt]

    nμ_1 = 1

    nμ_2 = 1

    for j = 1:max(Nμ_1, Nμ_2)
        if (nμ_1 < Nμ_1)
            nμ_1 += 1

            Vec_μ = [Vec_μ; μ_1]

            σ_opt = Optim.minimizer(optimize(f2, 1e-5, 7))

            Vec_σ = [Vec_σ; σ_opt]

            Vec_σ1 = [Vec_σ1; σ_opt]
        end

        if (nμ_2 < Nμ_2)
            nμ_2 += 1

            Vec_μ = [Vec_μ; μ_2]

            σ_opt = Optim.minimizer(optimize(f2, 1e-5, 7))

            Vec_σ = [Vec_σ; σ_opt]

            Vec_σ2 = [Vec_σ2; σ_opt]
        end
    end
    return [Vec_σ1; Vec_σ2]
end

##
EoF = "✓"
