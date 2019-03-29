#!/usr/bin/env julia

include("RunningStatistics.jl")
include("PhysicsFunctions.jl")
using Random
using DataFrames
using CSV
using Base.Threads
using Future
include("Constants.jl")

function main()::Nothing
    set_zero_subnormals(true)

    # Iteration condition
    local gen_array::Array{MersenneTwister, 1} = let m::MersenneTwister = MersenneTwister(1234)
        @fastmath [m; accumulate(Future.randjump, fill(big(10)^20, nthreads() - 1), init=m)]
    end

    # Computational values
    local times::Vector{Float64} = @fastmath [(x * delta_t + t_init) * sol for x in 1:num_t]

    # Probability for material sampling
    local prob_1::Float64 = @fastmath chord_1 / (chord_1 + chord_2)
    local change_prob_1::Float64 = @fastmath 1.0 / chord_1 * delta_t
    local change_prob_2::Float64 = @fastmath 1.0 / chord_2 * delta_t
    if @fastmath((change_prob_1 > 1.0) || (change_prob_2 > 1.0))
        println("The value for delta_t is too large for sampling")
        return nothing
    end

    # Parallel arrays
    local stat_1_intensity::Array{RunningStatistics.RunningStat, 2} = RunningStatistics.threadarray(num_t)
    local stat_2_intensity::Array{RunningStatistics.RunningStat, 2} = RunningStatistics.threadarray(num_t)
    local stat_1_temp::Array{RunningStatistics.RunningStat, 2} = RunningStatistics.threadarray(num_t)
    local stat_2_temp::Array{RunningStatistics.RunningStat, 2} = RunningStatistics.threadarray(num_t)

    println("Proceeding with ", nthreads(), " computational threads...")

    local printlock::SpinLock = SpinLock()

    # Outer loop
    @threads for i = 1:max_iterations
        local rand_num::Float64

        # First loop uses initial conditions
        local (intensity_value::Float64, temp_value::Float64) = (init_intensity, init_temp)

        # Sample initial starting material
        rand_num = @inbounds @fastmath rand(gen_array[threadid()], Float64)
        local material_num::Int32 = @fastmath (rand_num < prob_1) ? 1 : 2

        # Inner loop
        for (index, time) in enumerate(times)
            local (opacity_term::Float64, spec_heat_term::Float64, dens::Float64, change_prob::Float64) = @fastmath (material_num == 1) ? (opacity_1, spec_heat_1, dens_1, change_prob_1) : (opacity_2, spec_heat_2, dens_2, change_prob_2)

            # Sample whether material changes
            rand_num = @inbounds @fastmath rand(gen_array[threadid()], Float64)

            if (rand_num > change_prob)
                local opacity::Float64 = PhysicsFunctions.sigma_a(opacity_term, temp_value)
                local spec_heat::Float64 = PhysicsFunctions.c_v(spec_heat_term, temp_value)

                local new_intensity_value::Float64 = PhysicsFunctions.balance_intensity(opacity, intensity_value, temp_value)
                local new_temp_value::Float64 = PhysicsFunctions.balance_temp(opacity, spec_heat, dens, intensity_value, temp_value)

                (intensity_value, temp_value) = (new_intensity_value, new_temp_value)

                if @fastmath(material_num == 1)
                    @inbounds RunningStatistics.push(stat_1_intensity[index, threadid()], intensity_value)  # erg/cm^2-s
                    @inbounds RunningStatistics.push(stat_1_temp[index, threadid()], temp_value)  # eV
                else
                    @inbounds RunningStatistics.push(stat_2_intensity[index, threadid()], intensity_value)  # erg/cm^2-s
                    @inbounds RunningStatistics.push(stat_2_temp[index, threadid()], temp_value)  # eV
                end
            else
                material_num = @fastmath (material_num == 1) ? 2 : 1
            end
        end

        # Need to reference Core namespace for thread-safe printing
        if @fastmath(i % num_say == 0)
            lock(printlock) do
                Core.println("History Number ", i)
            end
        end
    end

    # Total tallies for each quantity
    local num_1::Vector{Float64} = RunningStatistics.total(stat_1_intensity)
    local num_2::Vector{Float64} = RunningStatistics.total(stat_2_intensity)

    # Mean values
    local material_1_intensity::Vector{Float64} = RunningStatistics.compute_mean(stat_1_intensity, num_1)
    local material_2_intensity::Vector{Float64} = RunningStatistics.compute_mean(stat_2_intensity, num_2)
    local material_1_temp::Vector{Float64} = RunningStatistics.compute_mean(stat_1_temp, num_1)
    local material_2_temp::Vector{Float64} = RunningStatistics.compute_mean(stat_2_temp, num_2)

    # Variance values
    local variance_1_intensity::Vector{Float64} = RunningStatistics.compute_variance(stat_1_intensity, num_1, material_1_intensity)
    local variance_2_intensity::Vector{Float64} = RunningStatistics.compute_variance(stat_2_intensity, num_2, material_2_intensity)
    local variance_1_temp::Vector{Float64} = RunningStatistics.compute_variance(stat_1_temp, num_1, material_1_temp)
    local variance_2_temp::Vector{Float64} = RunningStatistics.compute_variance(stat_2_temp, num_2, material_2_temp)

    # Maximum and minimum values
    local max_intensity_1::Vector{Float64} = RunningStatistics.compute_max(stat_1_intensity)
    local min_intensity_1::Vector{Float64} = RunningStatistics.compute_min(stat_1_intensity)
    local max_intensity_2::Vector{Float64} = RunningStatistics.compute_max(stat_2_intensity)
    local min_intensity_2::Vector{Float64} = RunningStatistics.compute_min(stat_2_intensity)
    local max_temp_1::Vector{Float64} = RunningStatistics.compute_max(stat_1_temp)
    local min_temp_1::Vector{Float64} = RunningStatistics.compute_min(stat_1_temp)
    local max_temp_2::Vector{Float64} = RunningStatistics.compute_max(stat_2_temp)
    local min_temp_2::Vector{Float64} = RunningStatistics.compute_min(stat_2_temp)

    local tabular::DataFrame = DataFrame(time=times, intensity1=material_1_intensity, varintensity1=variance_1_intensity, maxintensity1=max_intensity_1, minintensity1=min_intensity_1, temperature1=material_1_temp, vartemperature1=variance_1_temp, maxtemperature1=max_temp_1, mintemperature1=min_temp_1, intensity2=material_2_intensity, varintensity2=variance_2_intensity, maxintensity2=max_intensity_2, minintensity2=min_intensity_2, temperature2=material_2_temp, vartemperature2=variance_2_temp, maxtemperature2=max_temp_2, mintemperature2=min_temp_2)

    CSV.write("out/nonlinear/data/nonlinearmc.csv", tabular)

    return nothing
end

main()
