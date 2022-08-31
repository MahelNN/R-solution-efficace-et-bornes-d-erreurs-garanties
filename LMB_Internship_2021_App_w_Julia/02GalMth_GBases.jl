######################################################################
## @Autor: Mahel Ndoumba N.                                         ##
## LMB internship 2021 : Efficient resolution and guaranteed bounds ##
##                         for eigenvalue problems in electronic    ##
##                                structure calculations            ##
######################################################################
include("02GalMth_GBases_h.jl")
##
# Initialization
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

x = -R:0.01:R

U_exa = Exact_sol.(x, λ_1, z_1, z_2, μ_1, μ_2)

p1 = plot(x, U_exa, label = "exact solution ");
# number of basis function per nucleus position
Nμ = 19
# number of supplementary basis functions on either side
Nμ1_sup_bases = 2

Nμ2_sup_bases = 0
# number of basis functions on either side
Nμ_1 = Nμ + Nμ1_sup_bases

Nμ_2 = Nμ + Nμ2_sup_bases
# Vector of necleus positions
Vec_μ = Vec_position(Nμ_1, Nμ_2, μ_1, μ_2)

# Find the optimal vector of standard deviations
Vec_σ = Opt_vec_σ_eig(Nμ_1, Nμ_2, R, z_1, z_2, μ_1, μ_2, 0.01)

@show Vec_σ

# Refinement of the numerical solution
h_sol = 0.01

λ_N, U_N = GM_BG_VD1D(Vec_μ, Vec_σ, μ_1, μ_2, z_1, z_2, R, h_sol)

@show λ_N

X_dis = -R:h_sol:R
# plot the numerical solution
p2 = plot!(X_dis, U_N, label = "num sol; Nμ_1 = $Nμ_1")


##
# Evolution of the numerical solution with respect to the numbers
#of Gaussian basis function
# Plot initialization
legend = :outerright

xtickfontsize = 14
ytickfontsize = 14
legendfontsize = 14

# parameter initialization
n_max = 15

z_1, z_2, μ_1, μ_2, λ_1, R = Init()

R = 4

nμ1_sup_bases = 0

nμ2_sup_bases = 0
# array for the plot objects
ps = Array{Any}(nothing, n_max)

h_sol = 0.01
X_exa = -R:0.01:R

# Exact solution vector
S_exa = Exact_sol.(X_exa, λ_1, z_1, z_2, μ_1, μ_2)

for n = 1:n_max
    nμ_1 = n + nμ1_sup_bases * (n != 1) * (n == 2)

    nμ_2 = n + nμ2_sup_bases * (n != 1) * ((n == 6) || (n == 7) || (n == 8))

    n_μ = nμ_1 + nμ_2

    Vec_μ = Vec_position(nμ_1, nμ_2, μ_1, μ_2)

    Vec_σ1 = Opt_vec_σ_eig(nμ_1, nμ_2, R, z_1, z_2, μ_1, μ_2, 0.01)


    display("σ = ")
    @show Vec_σ1[nμ_1], Vec_σ1[end]
    # @show Vec_σ1

    λ_gmd1, U_gmd1 = GM_BG_VD1D(Vec_μ, Vec_σ1, μ_1, μ_2, z_1, z_2, R, h_sol)

    display("λ_N = ")
    @show λ_gmd1

    #Plot numerical solutions
    X_gm = -R:h_sol:R

    ps[n] = plot(
        [X_gm, X_exa],
        [U_gmd1, S_exa],
        label = ["uₙ" "u"],
        color = [:blue :red],
        linestyle = [:dashdot :solid],
        xlabel = "x",
        ylabel = "u(x)",
        title = "Gaussian Bases ($n_μ basis fct); R=$R",
        legends = legend,
        xtickfontsize = xtickfontsize,
        ytickfontsize = xtickfontsize,
        legendfontsize = legendfontsize,
        framestyle = :box,
    )
    # savefig("03BGauss_eig_$n_μ")
    display(ps[n])
end

##
# Error analysis
# parameter initialization
dof_max = 10

nμ1_sup_bases = 0

nμ2_sup_bases = 0

Vec_dof = zeros(dof_max)

z_1, z_2, μ_1, μ_2, λ_1, R = Init()

R = 20

h_err = 0.01

X_err = -R:h_err:R

S_exa = Exact_sol.(X_err, λ_1, z_1, z_2, μ_1, μ_2)

# error arrays initialization
λ_Err_gm = zeros(dof_max)

U_H1err_gm = zeros(dof_max)

U_L2err_gm = zeros(dof_max)

nμ1_sup_bases = 1

nμ2_sup_bases = 0

for n = 1:dof_max
    nμ_1 = 2 * n + nμ1_sup_bases#*(n!=1)

    nμ_2 = 2 * n + nμ2_sup_bases

    n_μ = nμ_1 + nμ_2

    Vec_dof[n] = n_μ

    Vec_μ = Vec_position(nμ_1, nμ_2, μ_1, μ_2)

    Vec_σ2 = Opt_vec_σ_eig(nμ_1, nμ_2, R, z_1, z_2, μ_1, μ_2, 0.0051)

    λ_gmd, U_gmd = GM_BG_VD1D(Vec_μ, Vec_σ2, μ_1, μ_2, z_1, z_2, R, h_err)

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
label_dof = latexstring("dof_{max} = $dof_max; R = $R")
marker = [:auto :auto :auto :none :none]
legend = :outerright

xtickfontsize = 14
ytickfontsize = 14
legendfontsize = 10

dofticks = 1:0.5:6
λgmticks = 10 .^ (-0.5:-0.25:-6.0)
H1gmticks = 10 .^ (-0.0:-0.25:-2.0)
L2gmticks = 10 .^ (-0.0:-0.5:-6.0)


# Plots
p11 = plot(
    (Vec_dof) .^ 0.5,
    λ_Err_gm,
    c = :red,
    label = label_dof,
    marker = :hex,
    yscale = :log10,
    xticks = dofticks,
    yticks = λgmticks,
    title = titlestring1,
    xlabel = L" √(N )",
    ylabel = L" log(err)",
    # legends = legend,
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    framestyle = :box,
)
# savefig("04Err_lambda_Galerkin_eig.png")


p12 = plot(
    (Vec_dof) .^ 0.5,
    U_H1err_gm,
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
    # legends = legend,
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    framestyle = :box,
)
# savefig("04Err_UH1_Galerkin_eig.png")


p13 = plot(
    (Vec_dof) .^ 0.5,
    U_L2err_gm,
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
# savefig("04Err_UL2_Galerkin_eig.png")


# relative error curves comparison plots
p14 = plot(
    [
        Vec_dof .^ 0.5,
        Vec_dof .^ 0.5,
        Vec_dof .^ 0.5,
        Vec_dof .^ 0.5,
        Vec_dof[1:end] .^ 0.5,
    ],
    [
        (U_H1err_gm),
        (U_L2err_gm),
        (λ_Err_gm),
        0.75 .* exp.(-Vec_dof .^ 0.25),
        0.5 .* exp.(-Vec_dof .^ 0.5)[1:end],
    ],
    c = [:blue :green :red :gray :brown],
    label = [L"H^1-norm" L"L^2-norm " "err_λ" L"y = \frac{3}{4}e^{-x^{0.25}}" L"y = \frac{1}{2}e^{-√x}}"],
    marker = marker,
    linestyle = [:solid :solid :solid :dashdot :dash],
    yscale = :log10,
    xticks = 1:6,
    yticks = L2gmticks,
    title = "Optimization over the eigenvalues; R=$R",
    xlabel = L" √(N )",
    ylabel = L" log(err)",
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    legends = :bottomleft,
    framestyle = :box,
)
# savefig("04Err_Galerkin_eig.png")


# Display plots
display(p11)
display(p12)
display(p13)
display(p14)
##
