#!/usr/bin/env julia

include("RunningStatistics.jl")
include("PhysicsFunctions.jl")
using Random
using DataFrames
using CSV
include("Constants.jl")

function main()::Nothing
    # Iteration condition
    local cont_calc_outer::Bool = true
    local generator::MersenneTwister = MersenneTwister(1234)

    # Computational values
    local iteration_number::Int64 = 0
    local times::Vector{Float64} = [(x * delta_t + t_init) * sol for x in 1:num_t]

    # Probability for material sampling
    local prob_1::Float64 = chord_1 / (chord_1 + chord_2)
    local change_prob_1::Float64 = 1.0 / chord_1 * delta_t
    local change_prob_2::Float64 = 1.0 / chord_2 * delta_t
    if ((change_prob_1 > 1.0) || (change_prob_2 > 1.0))
        println("The value for delta_t is too large for sampling")
        return nothing
    end

    # Parallel arrays
    local stat_1_intensity::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:num_t]
    local stat_2_intensity::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:num_t]
    local stat_1_temp::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:num_t]
    local stat_2_temp::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:num_t]

    # Outer loop
    while (cont_calc_outer)
        local rand_num::Float64

        # First loop uses initial conditions
        local (intensity_value::Float64, temp_value::Float64) = (init_intensity, init_temp)

        # Sample initial starting material
        rand_num = rand(generator, Float64)
        local material_num::Int32 = (rand_num < prob_1) ? 1 : 2

        # Inner loop
        for (index, time) in enumerate(times)
            local (opacity_term::Float64, spec_heat_term::Float64, dens::Float64, change_prob::Float64) = (material_num == 1) ? (opacity_1, spec_heat_1, dens_1, change_prob_1) : (opacity_2, spec_heat_2, dens_2, change_prob_2)

            # Sample whether material changes
            rand_num = rand(generator, Float64)

            if (rand_num > change_prob)
                local opacity::Float64 = PhysicsFunctions.sigma_a(opacity_term, temp_value)
                local spec_heat::Float64 = PhysicsFunctions.c_v(spec_heat_term, temp_value)

                local new_intensity_value::Float64 = PhysicsFunctions.balance_intensity(opacity, delta_t, intensity_value, temp_value)
                local new_temp_value::Float64 = PhysicsFunctions.balance_temp(opacity, spec_heat, dens, delta_t, intensity_value, temp_value)

                (intensity_value, temp_value) = (new_intensity_value, new_temp_value)

                if (material_num == 1)
                    RunningStatistics.push(stat_1_intensity[index], intensity_value)  # erg/cm^2-s
                    RunningStatistics.push(stat_1_temp[index], temp_value)  # eV
                else
                    RunningStatistics.push(stat_2_intensity[index], intensity_value)  # erg/cm^2-s
                    RunningStatistics.push(stat_2_temp[index], temp_value)  # eV
                end
            else
                material_num = (material_num == 1) ? 2 : 1
            end
        end

        iteration_number += 1

        if (iteration_number % num_say == 0)
            println("History Number ", iteration_number)
        end

        if (iteration_number > max_iterations)
            cont_calc_outer = false
        end
    end

    # Mean values
    local material_1_intensity::Vector{Float64} = RunningStatistics.mean.(stat_1_intensity)
    local material_2_intensity::Vector{Float64} = RunningStatistics.mean.(stat_2_intensity)
    local material_1_temp::Vector{Float64} = RunningStatistics.mean.(stat_1_temp)
    local material_2_temp::Vector{Float64} = RunningStatistics.mean.(stat_2_temp)

    # Variance values
    local variance_1_intensity::Vector{Float64} = RunningStatistics.variance.(stat_1_intensity)
    local variance_2_intensity::Vector{Float64} = RunningStatistics.variance.(stat_2_intensity)
    local variance_1_temp::Vector{Float64} = RunningStatistics.variance.(stat_1_temp)
    local variance_2_temp::Vector{Float64} = RunningStatistics.variance.(stat_2_temp)

    # Maximum and minimum values
    local max_intensity_1::Vector{Float64} = RunningStatistics.greatest.(stat_1_intensity)
    local min_intensity_1::Vector{Float64} = RunningStatistics.least.(stat_1_intensity)
    local max_intensity_2::Vector{Float64} = RunningStatistics.greatest.(stat_2_intensity)
    local min_intensity_2::Vector{Float64} = RunningStatistics.least.(stat_2_intensity)
    local max_temp_1::Vector{Float64} = RunningStatistics.greatest.(stat_1_temp)
    local min_temp_1::Vector{Float64} = RunningStatistics.least.(stat_1_temp)
    local max_temp_2::Vector{Float64} = RunningStatistics.greatest.(stat_2_temp)
    local min_temp_2::Vector{Float64} = RunningStatistics.least.(stat_2_temp)

    local tabular::DataFrame = DataFrame(time=times, intensity1=material_1_intensity, varintensity1=variance_1_intensity, maxintensity1=max_intensity_1, minintensity1=min_intensity_1, temperature1=material_1_temp, vartemperature1=variance_1_temp, maxtemperature1=max_temp_1, mintemperature1=min_temp_1, intensity2=material_2_intensity, varintensity2=variance_2_intensity, maxintensity2=max_intensity_2, minintensity2=min_intensity_2, temperature2=material_2_temp, vartemperature2=variance_2_temp, maxtemperature2=max_temp_2, mintemperature2=min_temp_2)

    CSV.write("out/nonlinear/data/nonlinearmc.csv", tabular)

    return nothing
end

main()
