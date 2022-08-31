
include("01FEM_h.jl")
##
# Initialization
function Init()
    z_1 = 1
    z_2 = 1

    μ_1 = -1
    μ_2 = 1

    λ_1 = λ_exact(z_1, z_2, μ_1, μ_2)

    R = 6
    return [z_1, z_2, μ_1, μ_2, λ_1, R]
end

z_1, z_2, μ_1, μ_2, λ_1, R = Init()

x = -R:0.01:R

# Exact solution vector and its plot
U_sol = Exact_sol.(x, λ_1, z_1, z_2, μ_1, μ_2)

p1 = plot(x, U_sol, label = "exact solution ");

# Numeical solution and its plot
N = 65

result = FEM_1D1(N, R, z_1, z_2, μ_1, μ_2)

h = 2 * R / N

x_fem = -R:h:R
p2 = plot!(x_fem, result[2], label = "num sol");

# display plots
display(p2)


##
# Evolution of the numerucal solution
z_1, z_2, μ_1, μ_2, λ_1, R = Init()

l_max = 6;
R = 10

for l = 1:l_max
    # # % Number of cells
    n = 2^(l + 1)
    # # % Meshsize
    h = 2 * R / n
    # % Positions of the vertices
    X_v = -R:h:R
    λ_dis, U_dis = FEM_1D1(n, R, z_1, z_2, μ_1, μ_2)

    # % Plot the exact and discrete solutions
    X_exa = -R:0.01:R

    X_dis = -R:h:R

    S_exa = Exact_sol.(X_exa, λ_1, z_1, z_2, μ_1, μ_2)

    p_i = plot(
        [X_dis, X_exa],
        [U_dis, S_exa],
        label = ["uₙ" "u"],
        c = [:blue :red],
        marker = [:circle :none],
        markercolor = :transparent,
        title = "Elements finis; dof = $(n+1) ; R = $R (V_eff = δ)",
        xlabel = "x",
        ylabel = "u(x)",
        framestyle = :box,
    )
    display(p_i)
    # savefig("00Elt_finis_ef_$n")
end


##
# Error analysis
#Initialization
z_1, z_2, μ_1, μ_2, λ_1, R = Init()

N_R = 5

Vec_R = zeros(N_R)

N_h = N_R
Vec_h = zeros(N_h)

λ_Err_ef = zeros(N_h, N_R)
U_H1err_ef = zeros(N_h, N_R)
U_L2err_ef = zeros(N_h, N_R)
Vec_res_ef = zeros(N_h, N_R)

for l = 1:N_h
    d_loc = 2
    h = 1 / 2^(l)
    Vec_h[l] = h

    for r = 1:N_R
        R = 3 * (r + 1) * μ_2
        Vec_R[r] = R
        N = ceil(Int64, 2 * R / h)

        X_v = -R:h:R

        S_exa = Exact_sol.(X_v, λ_1, z_1, z_2, μ_1, μ_2)
        λ_ef, U_ef = FEM_1D1(N, R, z_1, z_2, μ_1, μ_2)

        # eigenvalue relative error
        λ_Err_ef[l, r] = (λ_ef - λ_1)
        # eigen vector L^2 and H^1-norm relative error
        U_L2err_ef[l, r], U_H1err_ef[l, r] =
            Err_L2_H1(U_ef, h, N, R, λ_1, z_1, z_2, μ_1, μ_2)
    end
end



##
# Plot initialization
titlestring1 = latexstring("λ_N - λ")
titlestring2 = latexstring("||U_N - U||_{H^1}")
titlestring3 = latexstring("||U_N - U||_{L^2}")

label_R = [" R = $R" for R in round.(Vec_R', digits = 2)]
label_h = [" h = $h" for h in round.(Vec_h', digits = 3)]

marker = [:d :star :square :x :h :+ :v]
legend = :outerright

xtickfontsize = 14
ytickfontsize = 14
legendfontsize = 11

# personalized ticks for the plots
hticks = 10 .^ (-3:0.5:2.25)
# Rticks = 10 .^(0.25:0.125:1.75)

H1ticks = 10 .^ (3:-0.5:-3.5)
H1pticks = 10 .^ (3:-0.5:-3.5)
# H1ticks = 10 .^(-1:-0.5:-3.25)
# H1pticks = 10 .^(-0.5:-0.5:-3.5)

L2ticks = 10 .^ (-1.0:-1:-7)
L2pticks = 10 .^ (-0.0:-1:-8)

λticks = 10 .^ (3:-1.0:-7)
λpticks = 10 .^ (3:-1.0:-7)

p51 = plot(
    Vec_h,
    U_H1err_ef,
    label = label_R,
    marker = marker,
    title = titlestring2,
    xlabel = "log(h)",
    ylabel = "log(err)",
    scale = :log10,
    xticks = hticks,
    yticks = H1ticks,
    legend = legend,
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    framestyle = :box,
)
# savefig("01Err_H1_h_ef")

p52 = plot(
    Vec_h,
    U_L2err_ef,
    label = label_R,
    marker = marker,
    title = titlestring3,
    xlabel = "log(h)",
    ylabel = "log(err)",
    scale = :log10,
    xticks = hticks,
    yticks = L2ticks,
    legend = legend,
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    framestyle = :box,
)
# savefig("01Err_L2_h_ef")


p53 = plot(
    Vec_h,
    λ_Err_ef,
    marker = marker,
    title = titlestring1,
    label = label_R,
    xlabel = "log(h)",
    ylabel = "log(err)",
    scale = :log10,
    xticks = hticks,
    yticks = λticks,
    legend = legend,
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    framestyle = :box,
)
# savefig("01Err_lambda_h_ef")


p54 = plot(
    Vec_R,
    U_H1err_ef',
    label = label_h,
    marker = marker,
    title = titlestring2,
    xlabel = "R",
    ylabel = "log(err)",
    yscale = :log10,
    yticks = H1pticks,
    legend = legend,
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    framestyle = :box,
)
# savefig("01Err_H1_R_ef")


p55 = plot(
    Vec_R,
    U_L2err_ef',
    label = label_h,
    marker = marker,
    title = titlestring3,
    xlabel = "R",
    ylabel = "log(err)",
    yscale = :log10,
    # xticks = Rticks,
    yticks = L2pticks,
    legend = legend,
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    framestyle = :box,
)
# savefig("01Err_L2_R_ef")


p56 = plot(
    Vec_R,
    λ_Err_ef',
    title = titlestring1,
    marker = marker,
    label = label_h,
    xlabel = "R",
    ylabel = "log(err)",
    yscale = :log10,
    yticks = λpticks,
    legend = legend,
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    framestyle = :box,
)
# savefig("01Err_lambda_R_ef")


Y57 = [
    U_H1err_ef[:, end],
    U_L2err_ef[:, end],
    λ_Err_ef[:, end],
    0.3 .* Vec_h,
    0.005 .* Vec_h .^ 2,
]
linestyle = [:solid :solid :solid :dash :dot]
p57 = plot(
    Vec_h,
    Y57,
    c = [:red :green :blue :gray :brown],
    marker = [:d :x :o :none :none],
    markercolor = [:red :green :transparent],
    linestyle = linestyle,
    title = ["Comparison (R = $R)" for R in Vec_R[end]],
    label = ["H¹" "L²" "err_λ" "slope=1" "slope=2"],
    xlabel = "log(h)",
    ylabel = "log(err)",
    scale = :log10,
    xticks = hticks,
    yticks = L2ticks,
    legend = :bottomright,
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    framestyle = :box,
)
# savefig("02Err_L2H1Lambda_h_ef")


# Error function with respect to the radius of the domain
err_R = Err_R.(Vec_R, λ_1, z_1, z_2, μ_1, μ_2)

Y58 = [
    U_H1err_ef'[:, end],
    U_L2err_ef'[:, end],
    λ_Err_ef'[:, end],
    0.4 .* err_R,
    err_R[1:end-3] .^ 2,
    # err_R .^ 2,
]


p58 = plot(
    [Vec_R, Vec_R, Vec_R, Vec_R, Vec_R[1:end-3]],
    # [Vec_R, Vec_R, Vec_R, Vec_R, Vec_R],
    Y58,
    c = [:red :green :blue :gray :brown],
    marker = [:d :x :o :none :none],
    linestyle = linestyle,
    title = ["Comparison (h = $h)" for h in Vec_h[end]],
    label = [L"H¹" L"L²" "err_λ" L"Slope=-1" L"Slope=-2"],
    xlabel = "R",
    ylabel = "log(err)",
    yscale = :log10,
    yticks = 10 .^ (3:-1.0:-20),
    legend = :bottomleft,
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    framestyle = :box,
)
# savefig("02Err_L2H1Lambda_R_ef")

# Display plots
display(p51)
display(p52)
display(p53)
display(p54)
display(p55)
display(p56)
display(p57)
display(p58)


##
# The residual and its dual norm for finite elements scheme
z_1, z_2, μ_1, μ_2, λ_1, R = Init()
R = 10

N_R = 4

Vec_R = zeros(N_R)

N_h = N_R
Vec_h = zeros(N_h)

U_H1err_ef = zeros(N_h)
res_norm_hm1 = zeros(N_h)

for l = 1:N_h
    d_loc = 2
    h = 1 / 2^(l)
    Vec_h[l] = h

    N = ceil(Int64, 2 * R / h)

    X_v = -R:h:R

    S_exa = Exact_sol.(X_v, λ_1, z_1, z_2, μ_1, μ_2)

    λ_ef, U_ef = FEM_1D1(N, R, z_1, z_2, μ_1, μ_2)

    U_H1err_ef[l] = Err_L2_H1(U_ef, h, N, R, λ_1, z_1, z_2, μ_1, μ_2)[2]

    res = Res(λ_ef, U_ef, R, z_1, z_2, μ_1, μ_2)

    res_norm_hm1[l] = res[1]

    Vec_res = res[2]

    p_i = plot(
        [X_v],
        [Vec_res],
        label = "Residual",
        title = "Finite elements; dof = $(N+1) ; R = $R (V_eff = δ)",
        xlabel = "x",
        legend = :top,
        xtickfontsize = xtickfontsize,
        ytickfontsize = xtickfontsize,
        legendfontsize = legendfontsize,
        framestyle = :box,
    )
    display(p_i)
end


##
# Plots initialization
titlestring1 = latexstring("||res(U_N, λ_N)||_{V'}")
titlestring2 = latexstring("||U_N - U||_{H^1}")

label_R = [" R = $R" for R in round.(R', digits = 2)]
label_h = [" h = $h" for h in round.(Vec_h', digits = 3)]
legend = :bottomright

# Personaliezd ticks
hticks = 10 .^ (-0.0:-0.25:-3)
H1ticks = 10 .^ (-1.0:-0.125:-6)

p59 = plot(
    Vec_h,
    U_H1err_ef,
    label = label_R[end],
    marker = marker,
    title = titlestring2,
    xlabel = "log(h)",
    ylabel = "log(err)",
    scale = :log10,
    xticks = hticks,
    yticks = H1ticks,
    legend = legend,
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    framestyle = :box,
)


p510 = plot(
    Vec_h,
    res_norm_hm1,
    marker = marker,
    title = titlestring1,
    label = label_R,
    xlabel = "log(h)",
    ylabel = L"||res(U_N, λ_N)||_{-1}",
    scale = :log10,
    xticks = hticks,
    yticks = H1ticks,
    legend = legend,
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    framestyle = :box,
)



Y3 = [
    U_H1err_ef[:, end],
    res_norm_hm1[:, end],
    0.15 .* Vec_h,
    # Vec_h .^-1,
]
linestyle = [:solid :solid :dash :dot]
p511 = plot(
    Vec_h,
    Y3,
    c = [:red :green :gray :brown],
    marker = [:d :x :none :none],
    markercolor = [:red :green :transparent :transparent],
    linestyle = linestyle,
    title = "Comparison ",
    label = ["H¹-norm" L"||res(U_N, λ_N)||_{V'}" "slope=1"],
    xlabel = "log(h)",
    ylabel = "log(err)",
    scale = :log10,
    xticks = hticks,
    yticks = H1ticks,
    legend = legend,
    xtickfontsize = xtickfontsize,
    ytickfontsize = xtickfontsize,
    legendfontsize = legendfontsize,
    framestyle = :box,
)

# Display plots
display(p59)
display(p510)
display(p511)

##
