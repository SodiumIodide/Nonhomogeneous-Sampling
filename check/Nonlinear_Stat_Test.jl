using Random
using ProgressMeter

global const lam0_init = 1.0 / 10.0
global const lam0_final = 197.0 / 10.0
global const lam0_var = 1.0

global const lam1_init = 1.0 / 10.0
global const lam1_final = 21.0 / 10.0
global const lam1_var = 0.1

global const num_calcs = convert(Int64, 1e7)

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
        return -mean * log(rand(rng))
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
    elseif smaple_list[sample_type] == :gaussian
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

global vf_0 = RunningStat()
global vf_1 = RunningStat()

@showprogress 1 for calc_number in 1:num_calcs
    local lam00::Float64 = lambda_sample(sample_type0, lam0_init, lam0_final, lam0_var)
    local lam10::Float64 = lambda_sample(sample_type1, lam1_init, lam1_final, lam1_var)

    local p_0::Float64 = lam00 / (lam00 + lam10)
    local p_1::Float64 = 1.0 - p_0

    push!(vf_0, p_0)
    push!(vf_1, p_1)
end

println("Volume fraction 0: ", mean(vf_0))
println("Volume fraction 1: ", mean(vf_1))
