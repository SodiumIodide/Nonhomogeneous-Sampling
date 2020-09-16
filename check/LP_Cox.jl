using LinearAlgebra
using Plots
using CSV
using DataFrames
using Random
using ProgressMeter

global const lam0_init = 1.0 / 10.0
global const lam0_final = 197.0 / 10.0
global const lam0_var = 1.0
global const c0 = 0.0
global const XS00 = 10.0 / 99.0
global const XSs00 = c0 * XS00

global const lam1_init = 1.0 / 10.0
global const lam1_final = 21.0 / 10.0
global const lam1_var = 0.1
global const c1 = 0.0
global const XS10 = 100.0 / 11.0
global const XSs10 = c1 * XS10

global const x_L = 0.0  # left boundary x-coordinate
global const x_R = 10.0  # right boundary x-coordinate
global const N = 50  # # of segments composing system
global const cps = 20  # cells per segment

global const num_calcs = convert(Int64, 1e4)

# Corresponds to [1 2 3]
global const sample_list = [:uniform :exponential :gaussian]
global const sample_type0 = 1
global const sample_type1 = 1

global rng = MersenneTwister(1234)

function lambda_sample(sample_type::Int64, init::Float64, final::Float64, var::Float64)::Float64
    @inline function lambda_sample_uniform(init::Float64, final::Float64)::Float64
        local delta::Float64 = final - init
        return rand(rng) * delta + init
    end

    @inline function lambda_sample_exponential(init::Float64, final::Float64)::Float64
        local mean::Float64 = (final + init) / 2.0
        return -mean * log(rand(rng));
    end

    @inline function lambda_sample_gaussian(init::Float64, final::Float64, var::Float64)::Float64
        # Box-Muller Algorithm
        local result::Float64 = -1.0
        local mean::Float64 = (final + init) / 2.0
        while result < 0.0
            local radius::Float64 = sqrt(-2.0 * log(rand(rng)))
            local angle::Float64 = 2.0 * pi * rand(rng)
            local gaussian::Float64 = radius * sin(angle)
            result = mean + sqrt(var) * gaussian
        end
        return result
    end

    if sample_list[sample_type] == :uniform
        return lambda_sample_uniform(init, final)
    elseif sample_list[sample_type] == :exponential
        return lambda_sample_exponential(init, final)
    elseif sample_list[sample_type] == :gaussian
        return lambda_sample_gaussian(init, final, var)
    else
        println("Not a recognized sample type")
        exit()
    end
end

mutable struct RunningStat
    m_n::Int64
    m_oldM::Float64
    m_newM::Float64
    m_oldS::Float64
    m_newS::Float64
    m_min::Float64
    m_max::Float64
end

RunningStat()::RunningStat = RunningStat(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

@inline function clear!(rs::RunningStat)::Nothing
    rs = RunningStat()
    return nothing
end

function push!(rs::RunningStat, x::Float64)::Nothing
    if !isnan(x)
        rs.m_n += 1

        # See Knuth TAOCP vol 2, 3rd edition, page 232
        if rs.m_n == 1
            rs.m_oldM = rs.m_newM = x
            rs.m_oldS = 0.0
            rs.m_min = rs.m_max = x
        else
            rs.m_newM = rs.m_oldM + (x - rs.m_oldM) / rs.m_n
            rs.m_newS = rs.m_oldS + (x - rs.m_oldM) * (x - rs.m_newM)

            # Set up for next iteration
            rs.m_oldM = rs.m_newM
            rs.m_oldS = rs.m_newS

            # Min and max
            if x < rs.m_min
                rs.m_min = x
            elseif x > rs.m_max
                rs.m_max = x
            end
        end
    end
    return nothing
end

@inline function num(rs::RunningStat)::Int64
    return rs.m_n
end

@inline function mean(rs::RunningStat)::Float64
    return (rs.m_n > 0) ? rs.m_newM : 0.0
end

@inline function variance(rs::RunningStat)::Float64
    return (rs.m_n > 1) ? (rs.m_newS / (rs.m_n - 1)) : 0.0
end

@inline function least(rs::RunningStat)::Float64
    return rs.m_min
end

@inline function greatest(rs::RunningStat)::Float64
    return rs.m_max
end

# Formation of x vectors
global const Mn = N+1  # # of segments bounds points
global const mT = cps * N
global const cT = mT + 1
global const xn = collect(Float64, range(x_L, stop=x_R, length=N+1))  # vectof of segment bound values
global x_temp = zeros(Float64, cps+1, N)
global x_temp[:, 1] = collect(Float64, range(x_L, stop=xn[2], length=cps+1))
for i = 2:N
    global x_temp[:, i] = collect(Float64, range(xn[i], stop=xn[i+1], length=cps+1))
end
global x = zeros(Float64, mT+1, 1)
global x[1:cps+1] = x_temp[:,1]
for i = 2:N
    global x[((i-1)*cps+2):(i*cps+1)] = x_temp[2:(cps+1),i]
end

global result_flux = Array{RunningStat, 2}(undef, (2, cT))
global result_vf = Array{RunningStat, 2}(undef, (2, cT))
for i = 1:cT
    global result_flux[1, i] = RunningStat()
    global result_flux[2, i] = RunningStat()
    global result_vf[1, i] = RunningStat()
    global result_vf[2, i] = RunningStat()
end

println(num_calcs, " samples being run")

@showprogress 1 for calc_number in 1:num_calcs
    #if calc_number % 100 == 0
    #    print("\rCalculation number ", calc_number, " / ", num_calcs)
    #end

    # Material 0 Data
    local lam00::Float64 = lambda_sample(sample_type0, lam0_init, lam0_final, lam0_var)  # lambda_0

    # Material 1 Data
    local lam10::Float64 = lambda_sample(sample_type1, lam1_init, lam1_final, lam1_var)

    #local lam00::Float64 = (lam0_init + lam0_final) / 2.0
    #local lam10::Float64 = (lam1_init + lam1_final) / 2.0

    local XS0::Array{Float64, 2} = zeros(Float64, N, 1)
    local XS1::Array{Float64, 2} = zeros(Float64, N, 1)
    local lam0n::Array{Float64, 2} = zeros(Float64, N, 1)
    local lam1n::Array{Float64, 2} = zeros(Float64, N, 1)
    local lamcn::Array{Float64, 2} = zeros(Float64, N, 1)
    local XSs0::Array{Float64, 2} = zeros(Float64, N, 1)
    local XSs1::Array{Float64, 2} = zeros(Float64, N, 1)
    local XSR0::Array{Float64, 2} = zeros(Float64, N, 1)
    local XSR1::Array{Float64, 2} = zeros(Float64, N, 1)

    for n::Int64 = 1:N
        XS0[n] = XS00
        XS1[n] = XS10
        XSs0[n] = XSs00
        XSs1[n] = XSs10
        lam0n[n] = lam00
        lam1n[n] = lam10
        XSR0[n] = XS0[n] - XSs0[n] / 2.0 + 1.0 / lam0n[n]
        XSR1[n] = XS1[n] - XSs1[n] / 2.0 + 1.0 / lam1n[n]
        lamcn[n] = 1.0 ./ (1.0./lam0n[n] + 1.0./lam1n[n])
    end

    local Theta_n::Array{Float64, 2} = zeros(Float64, 4, Mn)
    local Theta::Array{Float64, 2} = zeros(Float64, 4, cT)
    local p_n::Array{Float64, 2} = zeros(Float64, 2, Mn)
    local p::Array{Float64, 2} = zeros(Float64, 2, cT)
    local psi_n::Array{Float64, 2} = zeros(Float64, 4, Mn)
    local psi::Array{Float64, 2} = zeros(Float64, 4, cT)

    # Boundary conditions

    # Use the calculated constant values
    p_n[1, 1] = lam0n[1] / (lam0n[1] + lam1n[1])
    p_n[2, 1] = 1.0 - p_n[1, 1]

    # Important: Order of Theta_i^{+/-,n} is Different from Report
    # Theta = < Theta_0+, Theta_1+, Theta_0-, Theta_1- >

    # Unit incident flux - values are vf-weighted
    Theta[1, 1] = p_n[1, 1]  # Theta_0+
    Theta[2, 1] = p_n[2, 1]  # Theta_1+
    Theta[3, cT] = 0.0  # Theta_0-
    Theta[4, cT] = 0.0  # Theta_1-

    local M::Array{Float64, 3} = zeros(Float64, 4, 4, N)
    local k::Array{Float64, 2} = zeros(Float64, 4, N)
    local xi0p::Array{Float64, 2} = zeros(Float64, 4, N)  # xi_0+
    local xi0m::Array{Float64, 2} = zeros(Float64, 4, N)  # xi_1+
    local xi1p::Array{Float64, 2} = zeros(Float64, 4, N)  # xi_0-
    local xi1m::Array{Float64, 2} = zeros(Float64, 4, N)  # xi_1-

    # Obtain eigenvalues and eigenvectors:
    for n::Int64 = 1:N
        M[:, :, n] = [ (-XSR0[n]) (1.0/lam1n[n]) (XSs0[n]/2.0) (0.0);  # xi_0+
                       (1.0/lam0n[n]) (-XSR1[n]) (0.0) (XSs1[n]/2.0);  # xi_1+
                       (-XSs0[n]/2.0) (0.0) (XSR0[n]) (-1.0/lam1n[n]);  # xi_0-
                       (0.0) (-XSs1[n]/2.0) (-1.0/lam0n[n]) (XSR1[n])]  # xi_1-
        local d::Array{Float64, 1}, xixi::Array{Float64, 2} = eigen(M[:, :, n])
        xixi = transpose(xixi)
        xi0p[:, n] = xixi[:, 1]
        xi1p[:, n] = xixi[:, 2]
        xi0m[:, n] = xixi[:, 3]
        xi1m[:, n] = xixi[:, 4]
        k[:, n] = d
    end

    local C::Array{Float64, 2} = zeros(Float64, 4*N, 4*N)
    local CC::Array{Float64, 3} = zeros(Float64, 4, 8, N-1)
    local B::Array{Float64, 2} = zeros(Float64, 4*N, 1)

    # Equations for the left boundary
    local en::Array{Float64, 1} = [exp(x_L*k[1,1]), exp(x_L*k[2,1]), exp(x_L*k[3,1]), exp(x_L*k[4,1])]
    C[1:2, 1:4] = [(xi0p[1,1]*en[1]) (xi0p[2,1]*en[2]) (xi0p[3,1]*en[3]) (xi0p[4,1]*en[4]);
                (xi1p[1,1]*en[1]) (xi1p[2,1]*en[2]) (xi1p[3,1]*en[3]) (xi1p[4,1]*en[4])]
    B[1:2] = [Theta[1,1] Theta[2,1]]
    # Equations for the right boundary
    local eN::Array{Float64, 1} = [exp(x_R*k[1,N]), exp(x_R*k[2,N]), exp(x_R*k[3,N]), exp(x_R*k[4,N])]
    C[4*N-1:4*N,4*N-3:4*N] = [(xi0m[1,N]*eN[1]) (xi0m[2,N]*eN[2]) (xi0m[3,N]*eN[3]) (xi0m[4,N]*eN[4]);
                            (xi1m[1,N]*eN[1]) (xi1m[2,N]*eN[2]) (xi1m[3,N]*eN[3]) (xi1m[4,N]*eN[4])]
    # The rest of the internal segment boundaries (4N-2 of them, each w/ 4 equations)
    for n::Int64 = 2:N
        local enmo::Array{Float64, 1} = [exp(xn[n]*k[1,n-1]), exp(xn[n]*k[2,n-1]), exp(xn[n]*k[3,n-1]), exp(xn[n]*k[4,n-1])]
        en = [exp(xn[n]*k[1,n]), exp(xn[n]*k[2,n]), exp(xn[n]*k[3,n]), exp(xn[n]*k[4,n])]

        CC[1,:,n-1] = [(xi0p[1,n-1]*enmo[1]) (xi0p[2,n-1]*enmo[2]) (xi0p[3,n-1]*enmo[3]) (xi0p[4,n-1]*enmo[4]) (-xi0p[1,n]*en[1]) (-xi0p[2,n]*en[2]) (-xi0p[3,n]*en[3]) (-xi0p[4,n]*en[4])]

        CC[2,:,n-1] = [(xi1p[1,n-1]*enmo[1]) (xi1p[2,n-1]*enmo[2]) (xi1p[3,n-1]*enmo[3]) (xi1p[4,n-1]*enmo[4]) (-xi1p[1,n]*en[1]) (-xi1p[2,n]*en[2]) (-xi1p[3,n]*en[3]) (-xi1p[4,n]*en[4])]

        CC[3,:,n-1] = [(xi0m[1,n-1]*enmo[1]) (xi0m[2,n-1]*enmo[2]) (xi0m[3,n-1]*enmo[3]) (xi0m[4,n-1]*enmo[4]) (-xi0m[1,n]*en[1]) (-xi0m[2,n]*en[2]) (-xi0m[3,n]*en[3]) (-xi0m[4,n]*en[4])]

        CC[4,:,n-1] = [(xi1m[1,n-1]*enmo[1]) (xi1m[2,n-1]*enmo[2]) (xi1m[3,n-1]*enmo[3]) (xi1m[4,n-1]*enmo[4]) (-xi1m[1,n]*en[1]) (-xi1m[2,n]*en[2]) (-xi1m[3,n]*en[3]) (-xi1m[4,n]*en[4])]
    end

    # Arrange CC's sheets into proper coefficient matrix positions
    for n::Int64 = 1:N-1
        C[4*n-1:4*n-1+3, 4*n-3:4*n-3+7] = CC[:, :, n]
    end

    # Obtain a_j^n coefficients from boundary condition matrix
    local aa::Array{Float64, 2} = C\B

    # Rearrange a_j^n coefficients for easier index-referencing later:
    local a::Array{Float64, 2} = zeros(Float64, 4, N)
    for n::Int64 = 1:N
        a[:,n] = aa[4*n-3:4*n-3+3]
    end

    local jj::Int64 = 0
    for i::Int64 = 1:cT
        # Determine which segment we are in
        for j::Int64 = 1:N
            if x[i] == xn[j]
                jj += 1
                break
            end
        end

        local axi0p::Array{Float64, 1} = [a[1,jj]*xi0p[1,jj], a[2,jj]*xi0p[2,jj], a[3,jj]*xi0p[3,jj], a[4,jj]*xi0p[4,jj]]
        local axi1p::Array{Float64, 1} = [a[1,jj]*xi1p[1,jj], a[2,jj]*xi1p[2,jj], a[3,jj]*xi1p[3,jj], a[4,jj]*xi1p[4,jj]]
        local axi0m::Array{Float64, 1} = [a[1,jj]*xi0m[1,jj], a[2,jj]*xi0m[2,jj], a[3,jj]*xi0m[3,jj], a[4,jj]*xi0m[4,jj]]
        local axi1m::Array{Float64, 1} = [a[1,jj]*xi1m[1,jj], a[2,jj]*xi1m[2,jj], a[3,jj]*xi1m[3,jj], a[4,jj]*xi1m[4,jj]]
        en = [exp(x[i]*k[1,jj]), exp(x[i]*k[2,jj]), exp(x[i]*k[3,jj]), exp(x[i]*k[4,jj])]

        Theta[1,i] = axi0p[1]*en[1] + axi0p[2]*en[2] + axi0p[3]*en[3] + axi0p[4]*en[4]  # Theta_0+
        Theta[2,i] = axi1p[1]*en[1] + axi1p[2]*en[2] + axi1p[3]*en[3] + axi1p[4]*en[4]  # Theta_1+
        Theta[3,i] = axi0m[1]*en[1] + axi0m[2]*en[2] + axi0m[3]*en[3] + axi0m[4]*en[4]  # Theta_0-
        Theta[4,i] = axi1m[1]*en[1] + axi1m[2]*en[2] + axi1m[3]*en[3] + axi1m[4]*en[4]  # Theta_1-
    end

    # Calculation of p_i^n and p_i (x):
    for n::Int64 = 2:Mn
        local dxn::Float64 = xn[n] - xn[n-1]

        # Calculation of p_i^n, i=0,1
        local a_vf::Float64 = dxn/lamcn[n-1]
        local b_vf::Float64 = p_n[2,n-1]*exp(-a_vf)
        local c_vf::Float64 = lamcn[n-1]/lam0n[n-1] * (1.0-exp(-a_vf))
        p_n[2,n] = b_vf + c_vf
        p_n[1,n] = 1.0 - p_n[2,n]
    end
    p[1,1] = p_n[1,1]
    p[2,1] = p_n[2,1]
    local nn::Int64 = 1
    for i::Int64 = 2:cT
        # Calculation of p_i(x), i=0,1
        local ae::Float64 = (x[i] - xn[nn])/lamcn[nn]
        local be::Float64 = p_n[2,nn]*exp(-ae)
        local ce::Float64 = lamcn[nn]/lam0n[nn]*(1.0 - exp(-ae))
        p[2,i] = be + ce
        p[1,i] = 1.0 - p[2,i]

        for k::Int64 = nn+1:Mn
            if x[i] == xn[k]
                nn += 1
                break
            end
        end
    end
    p[1,cT] = p_n[1,Mn]
    p[2,cT] = p_n[2,Mn]

    for i::Int64 = 1:cT
        psi[1,i] = Theta[1,i]/p[1,i]  # psi_0+
        psi[2,i] = Theta[2,i]/p[2,i]  # psi_1+
        psi[3,i] = Theta[3,i]/p[1,i]  # psi_0-
        psi[4,i] = Theta[4,i]/p[2,i]  # psi_1-
    end

    local flux::Array{Float64, 2} = zeros(Float64, 2, cT)
    for i::Int64 = 1:cT
        flux[1,i] = (psi[1,i] + psi[3,i])  # flux_0
        flux[2,i] = (psi[2,i] + psi[4,i])  # flux_1
    end

    # Average results
    for i::Int64 = 1:cT
        push!(result_flux[1, i], flux[1, i])
        push!(result_flux[2, i], flux[2, i])
        push!(result_vf[1, i], p[1, i])
        push!(result_vf[2, i], p[2, i])
    end
end

println("\nCalculation done, writing output...")

if sample_list[sample_type0] == :uniform
    file_name0 = "uniform"
elseif sample_list[sample_type0] == :exponential
    file_name0 = "exponential"
elseif sample_list[sample_type0] == :gaussian
    file_name0 = "gaussian"
end

if sample_list[sample_type1] == :uniform
    file_name1 = "uniform"
elseif sample_list[sample_type1] == :exponential
    file_name1 = "exponential"
elseif sample_list[sample_type1] == :gaussian
    file_name1 = "gaussian"
end

"""
plot(x, transpose(mean.(result_flux)), title="Flux", label=["Material 0" "Material 1"], color=[:red :blue], yaxis=:log10)
ylabel!("Flux")
xlabel!("x")
ylims!((1e-3, 1e1))
savefig(string("cox_output_", file_name0, "_", file_name1, ".png"))

plot(x, transpose(mean.(result_vf)), reuse=false, title="Volume Fraction", label=["Material 0" "Material 1"], color=[:red :blue])
ylabel!("Volume Fraction")
xlabel!("x")
ylims!((0.0, 1.0))
savefig(string("cox_volfrac_", file_name0, "_", file_name1, ".png"))
"""

tabular = DataFrame(x=vec(x), flux_0=vec(transpose(mean.(result_flux[1, :]))), flux_0_var=vec(transpose(variance.(result_flux[1, :]))), flux_1=vec(transpose(mean.(result_flux[2, :]))), flux_1_var=vec(transpose(variance.(result_flux[2, :]))), vf_0=vec(transpose(mean.(result_vf[1, :]))), vf_0_var=vec(transpose(variance.(result_vf[1, :]))), vf_1=vec(transpose(mean.(result_vf[2, :]))), vf_1_var=vec(transpose(variance.(result_vf[2, :]))))

CSV.write(string(file_name0, "_", file_name1, ".csv"), tabular)

println("Done")
