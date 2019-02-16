#!/usr/bin/env julia

include("GeometryGen.jl")
using .GeometryGen
include("MeshMap.jl")
using .MeshMap
include("RunningStatistics.jl")
using .RunningStatistics
using Random
using Future
using LinearAlgebra
using DataFrames
using CSV
using Base.Threads
include("Constants.jl")

function main()::Nothing
    # Iteration condition
    local gen_array::Array{MersenneTwister, 1} = let m::MersenneTwister = MersenneTwister(1234)
        [m; accumulate(Future.randjump, fill(big(10)^20, nthreads() - 1), init=m)]
    end

    # Computational values
    local times::Vector{Float64} = [(x * delta_t + t_init) * sol for x in 1:num_t]

    function sigma_a(opacity_term::Float64, temp::Float64)::Float64
        return opacity_term / temp^3
    end

    function c_v(spec_heat_term::Float64, temp::Float64)::Float64
        return spec_heat_term
    end

    function balance_intensity(opacity::Float64, past_intensity::Float64, past_temp::Float64)::Float64
        local term_1::Float64 = delta_t * sol^2 * opacity * arad * past_temp^4
        local term_2::Float64 = (1.0 - delta_t * sol * opacity) * past_intensity

        return term_1 + term_2
    end

    function balance_temp(opacity::Float64, spec_heat::Float64, density::Float64, past_intensity::Float64, past_temp::Float64)::Float64
        local term_1::Float64 = delta_t / (density * spec_heat) * opacity * past_intensity
        local term_2::Float64 = - delta_t / (density * spec_heat) * sol * opacity * arad * past_temp^4
        local term_3::Float64 = past_temp

        return term_1 + term_2 + term_3
    end

    # Parallel arrays
    local stat_1_intensity::Array{RunningStat, 2} = [RunningStat(0, 0.0, 0.0, 0.0, 0.0) for i in 1:nthreads(), j in 1:num_t]
    local stat_2_intensity::Array{RunningStat, 2} = [RunningStat(0, 0.0, 0.0, 0.0, 0.0) for i in 1:nthreads(), j in 1:num_t]
    local stat_1_temp::Array{RunningStat, 2} = [RunningStat(0, 0.0, 0.0, 0.0, 0.0) for i in 1:nthreads(), j in 1:num_t]
    local stat_2_temp::Array{RunningStat, 2} = [RunningStat(0, 0.0, 0.0, 0.0, 0.0) for i in 1:nthreads(), j in 1:num_t]

    println(string("Proceeding with ", nthreads(), " computational threads..."))

    local printlock::SpinLock = SpinLock()

    # Outer loop
    @threads for i = 1:max_iterations
        # Prevent random number clashing with discrete generators
        local (t_delta::Vector{Float64}, t_arr::Vector{Float64}, materials::Vector{Int32}, num_cells::Int64) = get_geometry(chord_1, chord_2, t_max, num_divs, rng=gen_array[threadid()])

        local intensity::Vector{Float64} = zeros(num_cells)
        local temp::Vector{Float64} = zeros(num_cells)

        # First loop uses initial conditions
        local (intensity_value::Float64, temp_value::Float64) = (init_intensity, init_temp)

        # Inner loop - Explicit Method
        for (index, material) in enumerate(materials)
            local delta_t_unstruct::Float64 = t_delta[index]
            local (opacity_term::Float64, spec_heat_term::Float64, dens::Float64) = (material == 1) ? (opacity_1, spec_heat_1, dens_1) : (opacity_2, spec_heat_2, dens_2)

            local opacity::Float64 = sigma_a(opacity_term, temp_value)
            local spec_heat::Float64 = c_v(spec_heat_term, temp_value)

            local new_intensity_value::Float64 = balance_intensity(opacity, intensity_value, temp_value)
            local new_temp_value::Float64 = balance_temp(opacity, spec_heat, dens, intensity_value, temp_value)

            (intensity_value, temp_value) = (new_intensity_value, new_temp_value)

            intensity[index] = intensity_value  # erg/cm^2-s
            temp[index] = temp_value  # eV
        end

        local material_intensity_array::Array{Float64, 2} = material_calc(intensity, t_delta, num_cells, materials, delta_t, num_t, convert(Int32, 2))
        local material_temp_array::Array{Float64, 2} = material_calc(temp, t_delta, num_cells, materials, delta_t, num_t, convert(Int32, 2))

        for k in 1:num_t
            if (material_intensity_array[k, 1] != 0.0)
                push(stat_1_intensity[threadid(), k], material_intensity_array[k, 1])  # erg/cm^2-s
                push(stat_1_temp[threadid(), k], material_temp_array[k, 1])  # eV
            else
                push(stat_2_intensity[threadid(), k], material_intensity_array[k, 2])  # erg/cm^2-s
                push(stat_2_temp[threadid(), k], material_temp_array[k, 2])  # eV
            end
        end

        # Need to reference Core namespace for thread-safe printing
        if (i % num_say == 0)
            lock(printlock) do
                Core.println(string("Iteration Number ", i))
            end
        end
    end

    local num_1::Vector{Float64} = vec(sum(convert.(Float64, num.(stat_1_intensity)), dims=1))
    local num_2::Vector{Float64} = vec(sum(convert.(Float64, num.(stat_2_intensity)), dims=1))

    function compute_mean(mat_r::Array{RunningStat, 2}, num_vec::Vector{Float64})::Vector{Float64}
        return vec(sum(mean.(mat_r) .* convert.(Float64, num.(mat_r)), dims=1) ./ num_vec')
    end

    local material_1_intensity::Vector{Float64} = compute_mean(stat_1_intensity, num_1)
    local material_2_intensity::Vector{Float64} = compute_mean(stat_2_intensity, num_2)
    local material_1_temp::Vector{Float64} = compute_mean(stat_1_temp, num_1)
    local material_2_temp::Vector{Float64} = compute_mean(stat_2_temp, num_2)

    function compute_variance(mat_r::Array{RunningStat, 2}, num_vec::Vector{Float64}, mean_vec::Vector{Float64})::Vector{Float64}
        local prefix::Vector{Float64} = vec((num_vec .- 1.0).^(-1))
        local first_sum::Vector{Float64} = vec(sum((convert.(Float64, num.(mat_r)) .- 1.0) .* variance.(mat_r), dims=1))
        local second_sum::Vector{Float64} = vec(sum(convert.(Float64, num.(mat_r)) .* (mean.(mat_r) .- mean_vec').^2, dims=1))

        return prefix .* (first_sum .+ second_sum)
    end

    local variance_1_intensity::Vector{Float64} = compute_variance(stat_1_intensity, num_1, material_1_intensity)
    local variance_2_intensity::Vector{Float64} = compute_variance(stat_2_intensity, num_2, material_2_intensity)
    local variance_1_temp::Vector{Float64} = compute_variance(stat_1_temp, num_1, material_1_temp)
    local variance_2_temp::Vector{Float64} = compute_variance(stat_2_temp, num_2, material_2_temp)

    local tabular::DataFrame = DataFrame(time=times, intensity1=material_1_intensity, varintensity1=variance_1_intensity, temperature1=material_1_temp, vartemperature1=variance_1_temp, intensity2=material_2_intensity, varintensity2=variance_2_intensity, temperature2=material_2_temp, vartemperature2=variance_2_temp)

    CSV.write("out/nonlinear/data/nonlinear_exp.csv", tabular)

    return nothing
end

main()
