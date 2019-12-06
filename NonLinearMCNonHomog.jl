#!/usr/bin/env julia

include("RunningStatistics.jl")
include("PhysicsFunctions.jl")
include("NonHomogChord.jl")
using Random
using DataFrames
using CSV
using ProgressMeter
include("Constants.jl")

function main()::Nothing
    set_zero_subnormals(true)

    # Generator
    local generator::MersenneTwister = MersenneTwister(1234)

    # Computational values
    local times::Vector{Float64} = @fastmath [x * delta_t + t_init for x in 1:num_t]

    # Parallel arrays
    local stat_1_intensity::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:num_t]
    local stat_2_intensity::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:num_t]
    local stat_1_temp::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:num_t]
    local stat_2_temp::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:num_t]
    local stat_1_opacity::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:num_t]
    local stat_2_opacity::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:num_t]

    # Problem parameters
    # Linear model slope
    local (slope_1::Float64, slope_2::Float64) = (0.0, 0.0)
    # Quadratic model equation parameters: a*x^2 + b*x + c
    local (param_a_1::Float64, param_b_1::Float64, param_c_1::Float64) = (0.0, 0.0, 0.0)
    local (param_a_2::Float64, param_b_2::Float64, param_c_2::Float64) = (0.0, 0.0, 0.0)
    # The value that the chord possesses for the maximum value computed (for rejection purposes)
    # As exponential distribution is computed using the inverse of the chord-length, the minimum value is chosen
    local (limiting_value_chord_1::Float64, limiting_value_chord_2::Float64) = (0.0, 0.0)
    if (quad)
        (param_a_1, param_b_1, param_c_1, param_a_2, param_b_2, param_c_2) = NonHomogChord.get_quad_params(start_chord_1, end_chord_1, start_chord_2, end_chord_2, t_max)
        limiting_value_chord_1 = minimum([start_chord_1, end_chord_1, NonHomogChord.mid_value_1_func(start_chord_1, end_chord_1)])
        limiting_value_chord_2 = minimum([start_chord_2, end_chord_2, NonHomogChord.mid_value_2_func(start_chord_2, end_chord_2)])
    else
        slope_1 = @fastmath (end_chord_1 - start_chord_1) / (t_max - t_init)
        slope_2 = @fastmath (end_chord_2 - start_chord_2) / (t_max - t_init)
        limiting_value_chord_1 = min(start_chord_1, end_chord_1)
        limiting_value_chord_2 = min(start_chord_2, end_chord_2)
    end

    # Outer loop
    @showprogress 1 for iteration_number = 1:max_iterations
        local rand_num::Float64

        # First loop uses initial conditions
        local (intensity_value::Float64, temp_value::Float64) = (init_intensity, init_temp)

        # Sample initial starting material
        # For nonhomogeneous chords, initial time is at 0.0 and so the probability is equivalent to the constant term ratio
        local prob_1::Float64 = @fastmath start_chord_1 / (start_chord_1 + start_chord_2)
        rand_num = @fastmath rand(generator, Float64)
        material_num = @fastmath (rand_num < prob_1) ? 1 : 2

        local time_value::Float64 = t_init

        # Inner loop
        for index = 1:num_t
            local (opacity_term::Float64, spec_heat_term::Float64, dens_term::Float64) = @fastmath (material_num == 1) ? (opacity_1, spec_heat_1, dens_1) : (opacity_2, spec_heat_2, dens_2)
            local (change_prob::Float64)

            if (quad)
                @inbounds @fastmath change_prob = (material_num == 1) ? delta_t / NonHomogChord.quad_chord(param_a_1, param_b_1, param_c_1, time_value) : delta_t / NonHomogChord.quad_chord(param_a_2, param_b_2, param_c_2, time_value)
            else
                @inbounds @fastmath change_prob = (material_num == 1) ? delta_t / NonHomogChord.linear_chord(start_chord_1, slope_1, time_value) : delta_t / NonHomogChord.linear_chord(start_chord_2, slope_2, time_value)
            end

            # Sample whether material changes
            rand_num = @fastmath rand(generator, Float64)

            if @fastmath (rand_num > change_prob)
                local opacity::Float64 = PhysicsFunctions.sigma_a(opacity_term, temp_value)
                local spec_heat::Float64 = PhysicsFunctions.c_v(spec_heat_term, temp_value)
                local dens::Float64 = PhysicsFunctions.rho(dens_term, temp_value)

                local new_intensity_value::Float64 = PhysicsFunctions.balance_intensity(opacity, delta_t, intensity_value, temp_value)
                local new_temp_value::Float64 = PhysicsFunctions.balance_temp(opacity, spec_heat, dens, delta_t, intensity_value, temp_value)

                (intensity_value, temp_value) = (new_intensity_value, new_temp_value)

                if @fastmath (material_num == 1)
                    @inbounds RunningStatistics.push(stat_1_intensity[index], intensity_value)  # erg/cm^2-s
                    @inbounds RunningStatistics.push(stat_1_temp[index], temp_value)  # eV
                    @inbounds RunningStatistics.push(stat_1_opacity[index], opacity)  # cm^-1
                else
                    @inbounds RunningStatistics.push(stat_2_intensity[index], intensity_value)  # erg/cm^2-s
                    @inbounds RunningStatistics.push(stat_2_temp[index], temp_value)  # eV
                    @inbounds RunningStatistics.push(stat_2_opacity[index], opacity)  # cm^-1
                end
            else
                material_num = @fastmath (material_num == 1) ? 2 : 1
            end
            time_value += delta_t  # s
        end
    end

    @fastmath times .*= sol

    # Mean values
    local material_1_intensity::Vector{Float64} = RunningStatistics.mean.(stat_1_intensity)
    local material_2_intensity::Vector{Float64} = RunningStatistics.mean.(stat_2_intensity)
    local material_1_temp::Vector{Float64} = RunningStatistics.mean.(stat_1_temp)
    local material_2_temp::Vector{Float64} = RunningStatistics.mean.(stat_2_temp)
    local material_1_opacity::Vector{Float64} = RunningStatistics.mean.(stat_1_opacity)
    local material_2_opacity::Vector{Float64} = RunningStatistics.mean.(stat_2_opacity)

    # Variance values
    local variance_1_intensity::Vector{Float64} = RunningStatistics.variance.(stat_1_intensity)
    local variance_2_intensity::Vector{Float64} = RunningStatistics.variance.(stat_2_intensity)
    local variance_1_temp::Vector{Float64} = RunningStatistics.variance.(stat_1_temp)
    local variance_2_temp::Vector{Float64} = RunningStatistics.variance.(stat_2_temp)
    local variance_1_opacity::Vector{Float64} = RunningStatistics.variance.(stat_1_opacity)
    local variance_2_opacity::Vector{Float64} = RunningStatistics.variance.(stat_2_opacity)

    # Maximum and minimum values
    local max_intensity_1::Vector{Float64} = RunningStatistics.greatest.(stat_1_intensity)
    local min_intensity_1::Vector{Float64} = RunningStatistics.least.(stat_1_intensity)
    local max_intensity_2::Vector{Float64} = RunningStatistics.greatest.(stat_2_intensity)
    local min_intensity_2::Vector{Float64} = RunningStatistics.least.(stat_2_intensity)
    local max_temp_1::Vector{Float64} = RunningStatistics.greatest.(stat_1_temp)
    local min_temp_1::Vector{Float64} = RunningStatistics.least.(stat_1_temp)
    local max_temp_2::Vector{Float64} = RunningStatistics.greatest.(stat_2_temp)
    local min_temp_2::Vector{Float64} = RunningStatistics.least.(stat_2_temp)
    local max_opacity_1::Vector{Float64} = RunningStatistics.greatest.(stat_1_opacity)
    local min_opacity_1::Vector{Float64} = RunningStatistics.least.(stat_1_opacity)
    local max_opacity_2::Vector{Float64} = RunningStatistics.greatest.(stat_2_opacity)
    local min_opacity_2::Vector{Float64} = RunningStatistics.least.(stat_2_opacity)

    local tabular::DataFrame = DataFrame(time=times, intensity1=material_1_intensity, varintensity1=variance_1_intensity, maxintensity1=max_intensity_1, minintensity1=min_intensity_1, temperature1=material_1_temp, vartemperature1=variance_1_temp, maxtemperature1=max_temp_1, mintemperature1=min_temp_1, opacity1=material_1_opacity, varopacity1=variance_1_opacity, maxopacity1=max_opacity_1, minopacity1=min_opacity_1, intensity2=material_2_intensity, varintensity2=variance_2_intensity, maxintensity2=max_intensity_2, minintensity2=min_intensity_2, temperature2=material_2_temp, vartemperature2=variance_2_temp, maxtemperature2=max_temp_2, mintemperature2=min_temp_2, opacity2=material_2_opacity, varopacity2=variance_2_opacity, maxopacity2=max_opacity_2, minopacity2=min_opacity_2)

    if (quad)
        CSV.write("out/nonlinear/data/nonlinearmc_nonhomog_quad.csv", tabular)
    else
        CSV.write("out/nonlinear/data/nonlinearmc_nonhomog_linear.csv", tabular)
    end

    return nothing
end

main()
