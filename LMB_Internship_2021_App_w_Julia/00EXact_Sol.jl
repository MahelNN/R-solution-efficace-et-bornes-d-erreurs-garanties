include("00Exact_Sol_h.jl")
################################ Main ##########################################
# Initialization
function Init()
    μ_1 = -1
    μ_2 = 1

    z_1 = 1
    z_2 = 1
    R = 3
    return [z_1, z_2, μ_1, μ_2, R]
end

z_1, z_2, μ_1, μ_2, R = Init()

λ_1 = λ_exact(z_1, z_2, μ_1, μ_2)
x = -R:0.01:R

# Exact solution vector
U_sol = Exact_sol.(x, λ_1, z_1, z_2, μ_1, μ_2)

# Plots of 3 setups of solution 
# Plot parameters
legend = :bottom
xtickfontsize = 14
ytickfontsize = 14
legendfontsize = 9.5


p1 = plot(
    x,
    U_sol,
    label = "u_exa: z_1 = $z_1, z_2 = $z_2 ",
    legend = legend,
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    framestyle = :box,
)
# savefig("sol_exa11")

# 2nd setup
z_1, z_2 = [2, 1]

λ_1 = λ_exact(z_1, z_2, μ_1, μ_2)
x = -R:0.01:R

U_sol = Exact_sol.(x, λ_1, z_1, z_2, μ_1, μ_2)

p2 = plot(
    x,
    U_sol,
    label = "u_exa: z_1 = $z_1, z_2 = $z_2 ",
    legend = legend,
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    framestyle = :box,
)
# savefig("sol_exa21")

# 3rd setup
z_1, z_2 = [1, 2]

λ_1 = λ_exact(z_1, z_2, μ_1, μ_2)
x = -R:0.01:R

U_sol = Exact_sol.(x, λ_1, z_1, z_2, μ_1, μ_2)

p3 = plot!(
    x,
    U_sol,
    label = "u_exa: z_1 = $z_1, z_2 = $z_2 ",
    legend = legend,
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    framestyle = :box,
)
# savefig("sol_exa12")

# Display plots
display(p1)
display(p2)
display(p3)
