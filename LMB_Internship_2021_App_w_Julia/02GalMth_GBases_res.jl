######################################################################
## @Autor: Mahel Ndoumba N.                                         ##
## LMB internship 2021 : Efficient resolution and guaranteed bounds ##
##                         for eigenvalue problems in electronic    ##
##                                structure calculations            ##
######################################################################

include("02GalMth_GBases_h.jl")
##
function Init()
    z_1 = 1
    z_2 = 1

    μ_1 = -1
    μ_2 = 1

    λ_1 = λ_exact(z_1, z_2, μ_1, μ_2)

    R = 10
    return [z_1, z_2, μ_1, μ_2, λ_1, R]
end

z_1, z_2, μ_1, μ_2, λ_1, R = Init()

x = -R:0.1:R

U_exa = Exact_sol.(x, λ_1, z_1, z_2, μ_1, μ_2)
# number of basis function per nucleus position
Nμ = 7

Nμ1_sup_bases = 0

Nμ2_sup_bases = 0
# number of supplementary basis functions on either side
Nμ_1 = Nμ + Nμ1_sup_bases

Nμ_2 = Nμ + Nμ2_sup_bases
# Vector of necleus positions
Vec_μ = Vec_position(Nμ_1, Nμ_2, μ_1, μ_2)

# Find the optimal vector of standard deviations
Vec_σ = Opt_vec_σ_res(Nμ_1, Nμ_2, R, z_1, z_2, μ_1, μ_2, 0.1)

@show Vec_σ

# Refinement of the numerical solution
h_sol = 0.1

λ_N, U_N = GM_BG_VD1D(Vec_μ, Vec_σ, μ_1, μ_2, z_1, z_2, R, h_sol)

@show λ_N

# plot the numerical solution
X_dis = -R:h_sol:R

p2 = plot(X_dis, U_N, label = "num sol; Nμ_1 = $Nμ_1");
p1 = plot!(x, U_exa, label = "exact solution ")


##
# Evolution of the numerical solution with respect to the numbers
#of Gaussian basis function
# Plot initialization
legend = :outerright

xtickfontsize = 14
ytickfontsize = 14
legendfontsize = 14

# parameter initialization
n_max = 8

z_1, z_2, μ_1, μ_2, λ_1, R = Init()

R = 10

nμ1_sup_bases = 0

nμ2_sup_bases = 0
# array for the plot objects
ps = Array{Any}(nothing, n_max)

h_sol = 0.1
X_exa = -R:0.1:R

# Exact solution vector
S_exa = Exact_sol.(X_exa, λ_1, z_1, z_2, μ_1, μ_2)

for n = 1:n_max
    nμ_1 = n + nμ1_sup_bases#*(n == 11)

    nμ_2 = n + nμ2_sup_bases

    n_μ = nμ_1 + nμ_2

    Vec_μ = Vec_position(nμ_1, nμ_2, μ_1, μ_2)

    Vec_σ1 = Opt_vec_σ_res(nμ_1, nμ_2, R, z_1, z_2, μ_1, μ_2, 0.1)

    display("σ")
    @show Vec_σ1[nμ_1], Vec_σ1[end]
    # @show Vec_σ1

    λ_gmd1, U_gmd1 = GM_BG_VD1D(Vec_μ, Vec_σ1, μ_1, μ_2, z_1, z_2, R, h_sol)

    display("λ_N")
    @show λ_gmd1

    #Plot the numerical solutions
    X_gm = -R:h_sol:R

    ps[n] = plot(
        [X_gm, X_exa],
        [U_gmd1, S_exa],
        label = ["uₙ" "u"],
        color = [:blue :red],
        linestyle = [:dashdot :solid],
        xlabel = "x",
        ylabel = "u(x)",
        title = "Gaussian Bases: ($nμ_1 + $nμ_2 ) basis fcts;  R=$R",
        legends = legend,
        xtickfontsize = xtickfontsize,
        ytickfontsize = xtickfontsize,
        legendfontsize = legendfontsize,
        framestyle = :box,
    )
    # savefig("05BGauss_Res_$n_μ")
    # savefig("05BGauss_Res_1_$n_μ")
    display(ps[n])
end



##
# Error analysis
# parameter initialization
dof_max = 11

nμ1_sup_bases = 0
# nμ1_sup_bases[1]

nμ2_sup_bases = 0

Vec_dof = zeros(dof_max)

z_1, z_2, μ_1, μ_2, λ_1, R = Init()

R = 10

h_err = 0.1

X_err = -R:h_err:R

S_exa = Exact_sol.(X_err, λ_1, z_1, z_2, μ_1, μ_2)

# error arrays initialization
λ_Err_gm = zeros(dof_max)

U_H1err_gm = zeros(dof_max)

U_L2err_gm = zeros(dof_max)

for n = 1:dof_max
    nμ_1 = n + nμ1_sup_bases * (n == dof_max)

    nμ_2 = n + nμ2_sup_bases#*(n==dof_max)

    n_μ = nμ_1 + nμ_2

    Vec_dof[n] = n_μ

    Vec_μ = Vec_position(nμ_1, nμ_2, μ_1, μ_2)

    Vec_σ1 = Opt_vec_σ_res(nμ_1, nμ_2, R, z_1, z_2, μ_1, μ_2, 0.1)

    λ_gmd, U_gmd = GM_BG_VD1D(Vec_μ, Vec_σ1, μ_1, μ_2, z_1, z_2, R, h_err)

    @show λ_gmd

    λ_Err_gm[n] = λ_gmd - λ_1

    U_L2err_gm[n], U_H1err_gm[n] =
        Err_L2_H1(U_gmd, h_err, size(U_gmd)[1] - 1, R, λ_1, z_1, z_2, μ_1, μ_2)
end




##
# Plot initialization
titlestring1 = latexstring("λ_N - λ")
titlestring2 = latexstring("||U_N - U||_{H^1}")
titlestring3 = latexstring("||U_N - U||_{L^2}")
label_dof = [
    latexstring(
        "R = $R; dof_{max} = ($(dof_max+nμ1_sup_bases), $(dof_max+nμ2_sup_bases)) ",
    ) for R in R'
]

marker = [:auto :auto :auto :none :none]
legend = :outerright

xtickfontsize = 14
ytickfontsize = 14
legendfontsize = 9

dofticks = (1.0:0.5:6)
λgmticks = 10 .^ (-0.5:-0.25:-3.0)
H1gmticks = 10 .^ (-0.0:-0.125:-2.0)
L2gmticks = 10 .^ (-0.0:-0.25:-3.0)

# Plots
p11 = plot(
    (Vec_dof[1:end]) .^ 0.5,
    λ_Err_gm[1:end],
    c = :red,
    label = label_dof,
    marker = :hex,
    yscale = :log10,
    xticks = dofticks,
    yticks = λgmticks,
    title = titlestring1,
    xlabel = L" √(N )",
    ylabel = L" log(err)",
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    framestyle = :box,
)
# savefig("06Err_lambda_Galerkin_res.png")
# savefig("06Err_lambda_Galerkin_res2.png")


p12 = plot(
    Vec_dof[1:end] .^ 0.5,
    [U_H1err_gm[1:end]],
    c = :green,
    label = label_dof,
    marker = :d,
    yscale = :log10,
    xticks = dofticks,
    yticks = H1gmticks,
    legend = :best,
    title = titlestring2,
    xlabel = L" √(N )",
    ylabel = L" log(err)",
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    # legends=:outerright,
    framestyle = :box,
)
# savefig("06Err_UH1_Galerkin_res.png")
# savefig("06Err_UH1_Galerkin_res2.png")


p13 = plot(
    Vec_dof[1:end] .^ 0.5,
    U_L2err_gm[1:end],
    c = :blue,
    label = label_dof,
    marker = :x,
    yscale = :log10,
    xticks = dofticks,
    yticks = L2gmticks,
    title = titlestring3,
    xlabel = L" √(N )",
    ylabel = L" log(err)",
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    # legends=:outerright,
    framestyle = :box,
)
# savefig("06Err_UL2_Galerkin_res.png")
# savefig("06Err_UL2_Galerkin_res2.png")


p14 = plot(
    [Vec_dof[1:end] .^ 0.5],
    [
        (U_H1err_gm[1:end]),
        (U_L2err_gm[1:end]),
        (λ_Err_gm[1:end]),
        1.75 .* exp.(-Vec_dof[1:end] .^ 0.5),
        0.5 .* exp.(-Vec_dof[1:end] .^ 0.5),
    ],
    c = [:blue :green :red :gray :brown :black],
    label = [L"H^1-norm" L"L^2-norm " "err_λ" L"y = 1.75 \cdot e^{-x^{0.5}}" L"y = 0.5 \cdot e^{-√x}}"],
    marker = marker,
    linestyle = [:solid :solid :solid :dashdot :dash],
    yscale = :log10,
    xticks = dofticks,
    yticks = H1gmticks,
    # yticks = 10 .^(-0.:-0.5:-2),
    title = latexstring("σ~ optimization~ over~ ||Res(λ_N,u_N)||_{-1}; R = $R"),
    xlabel = L"√(N) ",
    ylabel = L" log(err)",
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    ylims = (L2gmticks[end-2], L2gmticks[2]),
    legends = :bottomleft,
    # legends=legend,
    framestyle = :box,
)
# savefig("06Err_Galerkin_res.png")
# savefig("06Err_Galerkin_res2.png")


# Display plots
display(p11)
display(p12)
display(p13)
display(p14)
##
