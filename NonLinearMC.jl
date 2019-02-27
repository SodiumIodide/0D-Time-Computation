#!/usr/bin/env julia

include("RunningStatistics.jl")
using Random
using DataFrames
using CSV
using Base.Threads
using Future
include("Constants.jl")

function main()::Nothing
    # Iteration condition
    local gen_array::Array{MersenneTwister, 1} = let m::MersenneTwister = MersenneTwister(1234)
        [m; accumulate(Future.randjump, fill(big(10)^20, nthreads() - 1), init=m)]
    end

    # Computational values
    local times::Vector{Float64} = [(x * delta_t + t_init) * sol for x in 1:num_t]

    # Probability for material sampling
    local prob_1::Float64 = chord_1 / (chord_1 + chord_2)
    local change_prob_1::Float64 = 1.0 / chord_1 * delta_t
    local change_prob_2::Float64 = 1.0 / chord_2 * delta_t
    if ((change_prob_1 > 1.0) || (change_prob_2 > 1.0))
        println("The value for delta_t is too large for sampling")
        return nothing
    end

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
    local stat_1_intensity::Array{RunningStatistics.RunningStat, 2} = RunningStatistics.threadarray(num_t)
    local stat_2_intensity::Array{RunningStatistics.RunningStat, 2} = RunningStatistics.threadarray(num_t)
    local stat_1_temp::Array{RunningStatistics.RunningStat, 2} = RunningStatistics.threadarray(num_t)
    local stat_2_temp::Array{RunningStatistics.RunningStat, 2} = RunningStatistics.threadarray(num_t)

    println(string("Proceeding with ", nthreads(), " computational threads..."))

    local printlock::SpinLock = SpinLock()

    # Outer loop
    @threads for i = 1:max_iterations
        local rand_num::Float64

        # First loop uses initial conditions
        local (intensity_value::Float64, temp_value::Float64) = (init_intensity, init_temp)

        # Sample initial starting material
        rand_num = rand(gen_array[threadid()], Float64)
        local material_num::Int32 = (rand_num < prob_1) ? 1 : 2

        # Inner loop
        for (index, time) in enumerate(times)
            local (opacity_term::Float64, spec_heat_term::Float64, dens::Float64, change_prob::Float64) = (material_num == 1) ? (opacity_1, spec_heat_1, dens_1, change_prob_1) : (opacity_2, spec_heat_2, dens_2, change_prob_2)

            # Sample whether material changes
            rand_num = rand(gen_array[threadid()], Float64)

            if (rand_num > change_prob)
                local opacity::Float64 = sigma_a(opacity_term, temp_value)
                local spec_heat::Float64 = c_v(spec_heat_term, temp_value)

                local new_intensity_value::Float64 = balance_intensity(opacity, intensity_value, temp_value)
                local new_temp_value::Float64 = balance_temp(opacity, spec_heat, dens, intensity_value, temp_value)

                (intensity_value, temp_value) = (new_intensity_value, new_temp_value)

                if (material_num == 1)
                    RunningStatistics.push(stat_1_intensity[index, threadid()], intensity_value)  # erg/cm^2-s
                    RunningStatistics.push(stat_1_temp[index, threadid()], temp_value)  # eV
                else
                    RunningStatistics.push(stat_2_intensity[index, threadid()], intensity_value)  # erg/cm^2-s
                    RunningStatistics.push(stat_2_temp[index, threadid()], temp_value)  # eV
                end
            else
                material_num = (material_num == 1) ? 2 : 1
            end
        end

        # Need to reference Core namespace for thread-safe printing
        if (i % num_say == 0)
            lock(printlock) do
                Core.println(string("Iteration Number ", i))
            end
        end
    end

    local num_1::Vector{Float64} = vec(sum(convert.(Float64, RunningStatistics.num.(stat_1_intensity)), dims=2))
    local num_2::Vector{Float64} = vec(sum(convert.(Float64, RunningStatistics.num.(stat_2_intensity)), dims=2))

    local material_1_intensity::Vector{Float64} = RunningStatistics.compute_mean(stat_1_intensity, num_1)
    local material_2_intensity::Vector{Float64} = RunningStatistics.compute_mean(stat_2_intensity, num_2)
    local material_1_temp::Vector{Float64} = RunningStatistics.compute_mean(stat_1_temp, num_1)
    local material_2_temp::Vector{Float64} = RunningStatistics.compute_mean(stat_2_temp, num_2)

    local variance_1_intensity::Vector{Float64} = RunningStatistics.compute_variance(stat_1_intensity, num_1, material_1_intensity)
    local variance_2_intensity::Vector{Float64} = RunningStatistics.compute_variance(stat_2_intensity, num_2, material_2_intensity)
    local variance_1_temp::Vector{Float64} = RunningStatistics.compute_variance(stat_1_temp, num_1, material_1_temp)
    local variance_2_temp::Vector{Float64} = RunningStatistics.compute_variance(stat_2_temp, num_2, material_2_temp)

    # Save maximum and minimum data
    local max_intensity_1::Float64 = RunningStatistics.compute_max(stat_1_intensity)
    local min_intensity_1::Float64 = RunningStatistics.compute_min(stat_1_intensity)
    local max_intensity_2::Float64 = RunningStatistics.compute_max(stat_2_intensity)
    local min_intensity_2::Float64 = RunningStatistics.compute_min(stat_2_intensity)
    local max_temp_1::Float64 = RunningStatistics.compute_max(stat_1_temp)
    local min_temp_1::Float64 = RunningStatistics.compute_min(stat_1_temp)
    local max_temp_2::Float64 = RunningStatistics.compute_max(stat_2_temp)
    local min_temp_2::Float64 = RunningStatistics.compute_min(stat_2_temp)

    local tabular::DataFrame = DataFrame(time=times, intensity1=material_1_intensity, varintensity1=variance_1_intensity, temperature1=material_1_temp, vartemperature1=variance_1_temp, intensity2=material_2_intensity, varintensity2=variance_2_intensity, temperature2=material_2_temp, vartemperature2=variance_2_temp)

    local minmax::DataFrame = DataFrame(maxint1=max_intensity_1, minint1=min_intensity_1, maxtemp1=max_temp_1, mintemp1=min_temp_1, maxint2=max_intensity_2, minint2=min_intensity_2, maxtemp2=max_temp_2, mintemp2=min_temp_2)

    CSV.write("out/nonlinear/data/nonlinearmc.csv", tabular)
    CSV.write("out/nonlinear/pdf_data/minmax.csv", minmax)

    return nothing
end

main()
