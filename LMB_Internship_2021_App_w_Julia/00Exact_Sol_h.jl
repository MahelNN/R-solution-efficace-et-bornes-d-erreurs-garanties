######################################################################
## @Autor: Mahel Ndoumba N.                                         ##
## LMB internship 2021 : Efficient resolution and guaranteed bounds ##
##                         for eigenvalue problems in electronic    ##
##                                structure calculations            ##
######################################################################
using LinearAlgebra, Distributions, SparseArrays
using Plots, LaTeXStrings

pyplot()
################## Subroutines exact solution ##################################
##
# Matrix from the eigenvalue problem
function sys_Mλ(λ, z_1, z_2, μ_1, μ_2)
    #Initialization
    α = √(-λ)
    M = zeros(4, 4)
    #ligne 1
    M[1, 1] = 1
    M[1, 2] = -exp(-2 * α * μ_1)
    M[1, 3] = -1
    #ligne 2
    M[2, 1] = α - z_1
    M[2, 2] = α * exp(-2 * α * μ_1)
    M[2, 3] = -α
    #ligne 3
    M[3, 2] = 1
    M[3, 3] = exp(2 * α * μ_2)
    M[3, 4] = -1
    #ligne 4
    M[4, 2] = -α
    M[4, 3] = α * exp(2 * α * μ_2)
    M[4, 4] = α - z_2
    return M
end


# Dichotomy method to find the root of a function
function dicho(f, a, b, iter)
    #Initialization
    if f(a) * f(b) <= 0
        for i = 1:iter
            c = (a + b) / 2
            if f(a) * f(c) <= 0 <= 0
                b = c
            else
                a = c
            end
        end
    else
        println(" [$a,$b] is not a good interval")
    end
    return [a, b]
end


# Find λ_exact such that det(sys_Mλ) = 0
function λ_exact(z_1, z_2, μ_1, μ_2)
    λ_min = -0.5 * (z_1^2 + z_2^2)
    # interval in which we will look for λ_exact
    x_λ = λ_min:0.01:-1e-6
    eps = z_1 * z_2

    Expl_λ(λ) = λ * (det(sys_Mλ(λ, z_1, z_2, μ_1, μ_2)) <= 1e-3)
    potential_λ = Expl_λ.(x_λ) # we only need a rough gest of the first value such that det(M) = 0
    # accurate λ_exact with Dichotomy method

    x_min = potential_λ[(potential_λ).!=0][1]
    ϵ = 0.1
    a = x_min - ϵ
    b = x_min + ϵ

    detMλ(λ) = det(sys_Mλ(λ, z_1, z_2, μ_1, μ_2))
    while detMλ(a) * detMλ(b) > 0.0
        a += ϵ
        b += ϵ
    end

    iter = 150
    return dicho(detMλ, a, b, iter)[1]
end


function Pos_vec(U)
    U_abs = (mean(U) > 0) ? U : -U
    return U_abs
end


# Coefficients for the exact solution
function Coef_exact_sol(λ, z_1, z_2, μ_1, μ_2)
    M = sys_Mλ(λ, z_1, z_2, μ_1, μ_2)
    E = eigen(M)
    ind_λ_0 = findall(x -> abs(x) == minimum(abs.(E.values)), abs.(E.values))

    coef_AB = (norm.(E.vectors[:, ind_λ_0]))
    return coef_AB
end


# Constant of L^2 normalization for the exact solution
function Norm_L2_exa_sol(λ, z_1, z_2, μ_1, μ_2)
    # Initialization
    α = √(-λ)
    # Parameters Aᵢ, Bᵢ for the exact solution
    B1, A2, B2, A3 = Coef_exact_sol(λ, z_1, z_2, μ_1, μ_2)

    normL2_u1 = B1^2 * (exp(2 * α * μ_1)) / 2 / α
    normL2_u2 =
        A2^2 * (exp(-2 * α * μ_1) - exp(-2 * α * μ_2)) / 2 / α -
        B2^2 * (exp(2 * α * μ_1) - exp(2 * α * μ_2)) / 2 / α +
        2 * A2 * B2 * (μ_2 - μ_1)
    normL2_u3 = A3^2 * (exp(-2 * α * μ_2)) / 2 / α

    norm_L2_u = normL2_u1 + normL2_u2 + normL2_u3
    return √(norm_L2_u)
end


function Exact_sol(x, λ, z_1, z_2, μ_1, μ_2)
    # Initialization
    α = √(-λ)
    gen_fom(x, A, B) = A * exp(-α * x) + B * exp(α * x)
    # Parameters Aᵢ, Bᵢ for the exact solution
    B1, A2, B2, A3 = Coef_exact_sol(λ, z_1, z_2, μ_1, μ_2)
    # Sur [-∞,μ₁] u(x) = B₁*exp(√-λ*x)
    u1 = gen_fom(x, 0.0, B1) * (x <= μ_1)
    # # Sur [μ₁,μ₂]  u(x) = A₂*exp(-√-λ*x) + B₂*exp(√-λ*x)
    u2 = gen_fom(x, A2, B2) * (μ_1 < x) * (x <= μ_2)
    # # Sur [μ₂,+∞] u(x) = A₃*exp(-√-λ*x)
    u3 = gen_fom(x, A3, 0.0) * (μ_2 < x)

    norm_L2_u = Norm_L2_exa_sol(λ, z_1, z_2, μ_1, μ_2)
    u_exa = (u1 + u2 + u3) / norm_L2_u
    return u_exa
end


# first derivative of the exact solution
function dExact_sol(x, λ, z_1, z_2, μ_1, μ_2)
    # Initialization
    α = √(-λ)
    dgen_fom(x, A, B) = -α * A * exp(-α * x) + α * B * exp(α * x)
    # Paramètres Aᵢ, Bᵢ pour la solution exacte
    B1, A2, B2, A3 = Coef_exact_sol(λ, z_1, z_2, μ_1, μ_2)
    # Sur [-∞,μ₁] u(x) = B₁*exp(√-λ*x)
    du1 = dgen_fom(x, 0.0, B1) * (x <= μ_1)
    # # Sur [μ₁,μ₂]  u(x) = A₂*exp(-√-λ*x) + B₂*exp(√-λ*x)
    du2 = dgen_fom(x, A2, B2) * (μ_1 < x) * (x <= μ_2)
    # # Sur [μ₂,+∞] u(x) = A₃*exp(-√-λ*x)
    du3 = dgen_fom(x, A3, 0.0) * (μ_2 < x)

    norm_L2_u = Norm_L2_exa_sol(λ, z_1, z_2, μ_1, μ_2)

    du_exa = (du1 + du2 + du3) / norm_L2_u

    return du_exa
end
##
EoF = "✓"
