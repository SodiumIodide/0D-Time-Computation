#!/usr/bin/env julia

include("Constants.jl")
using Random
using DataFrames
using CSV
using Base.Threads

function main()::Nothing
    # Iteration condition
    local generator::MersenneTwister = MersenneTwister(1234)

    # Computational values
    local iteration_number::Atomic{Int64} = Atomic{Int64}(0)
    local material_1_intensity::Vector{Float64} = zeros(num_t)
    local material_2_intensity::Vector{Float64} = zeros(num_t)
    local material_1_hits::Vector{Float64} = zeros(num_t)
    local material_1_temp::Vector{Float64} = zeros(num_t)
    local material_2_temp::Vector{Float64} = zeros(num_t)
    local material_2_hits::Vector{Float64} = zeros(num_t)
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

    # Variables for variance computation
    local sum_1_intensity_square::Vector{Float64} = zeros(num_t)
    local sum_2_intensity_square::Vector{Float64} = zeros(num_t)
    local sum_1_temp_square::Vector{Float64} = zeros(num_t)
    local sum_2_temp_square::Vector{Float64} = zeros(num_t)
    local variance_1_intensity::Vector{Float64} = zeros(num_t)
    local variance_2_intensity::Vector{Float64} = zeros(num_t)
    local variance_1_temp::Vector{Float64} = zeros(num_t)
    local variance_2_temp::Vector{Float64} = zeros(num_t)

    println(string("Proceeding with ", nthreads(), " computational threads..."))

    # Locking conditions
    #local array_lock_1::SpinLock = SpinLock()
    #local array_lock_2::SpinLock = SpinLock()
    #local random_lock::SpinLock = SpinLock()

    # Outer loop
    #@threads for i = 1:max_iterations
    for i = 1:max_iterations
        local rand_num::Float64

        # First loop uses initial conditions
        local (intensity_value::Float64, temp_value::Float64) = (init_intensity, init_temp)

        # Sample initial starting material
        #lock(random_lock)
        rand_num = rand(generator, Float64)
        #unlock(random_lock)
        local material_num::Int32 = (rand_num < prob_1) ? 1 : 2

        # Inner loop
        for (index, time) in enumerate(times)
            local (opacity_term::Float64, spec_heat_term::Float64, dens::Float64, change_prob::Float64) = (material_num == 1) ? (opacity_1, spec_heat_1, dens_1, change_prob_1) : (opacity_2, spec_heat_2, dens_2, change_prob_2)

            # Sample whether material changes
            #lock(random_lock) do
            rand_num = rand(generator, Float64)
            #unlock(random_lock)

            if (rand_num > change_prob)
                local opacity::Float64 = sigma_a(opacity_term, temp_value)
                local spec_heat::Float64 = c_v(spec_heat_term, temp_value)

                local new_intensity_value::Float64 = balance_intensity(opacity, intensity_value, temp_value)
                local new_temp_value::Float64 = balance_temp(opacity, spec_heat, dens, intensity_value, temp_value)

                (intensity_value, temp_value) = (new_intensity_value, new_temp_value)

                if (material_num == 1)
                    #lock(array_lock_1)
                    material_1_intensity[index] += intensity_value  # erg/cm^2-s
                    material_1_temp[index] += temp_value  # eV
                    sum_1_intensity_square[index] += intensity_value^2
                    sum_1_temp_square[index] += temp_value^2
                    material_1_hits[index] += 1.0
                    #unlock(array_lock_1)
                else
                    #lock(array_lock_2)
                    material_2_intensity[index] += intensity_value  # erg/cm^2-s
                    material_2_temp[index] += temp_value  # eV
                    sum_2_intensity_square[index] += intensity_value^2
                    sum_2_temp_square[index] += temp_value^2
                    material_2_hits[index] += 1.0
                    #unlock(array_lock_2)
                end
            else
                material_num = (material_num == 1) ? 2 : 1
            end
        end

        atomic_add!(iteration_number, 1)

        if (iteration_number[] % num_say == 0)
            println(string("Iteration Number ", iteration_number[]))
        end
    end

    local max_iterations_f::Float64 = convert(Float64, max_iterations)
    local variance_prefix::Float64 = 1.0 / (max_iterations_f * (max_iterations_f - 1.0))
    variance_1_intensity = variance_prefix .* (max_iterations_f .* sum_1_intensity_square - material_1_intensity.^2)
    variance_2_intensity = variance_prefix .* (max_iterations_f .* sum_2_intensity_square - material_2_intensity.^2)
    variance_1_temp = variance_prefix .* (max_iterations_f .* sum_1_temp_square - material_1_temp.^2)
    variance_2_temp = variance_prefix .* (max_iterations_f .* sum_2_temp_square - material_2_temp.^2)

    material_1_intensity ./= material_1_hits
    material_2_intensity ./= material_2_hits
    material_1_temp ./= material_1_hits
    material_2_temp ./= material_2_hits

    tabular::DataFrame = DataFrame(time=times, intensity1=material_1_intensity, varintensity1=variance_1_intensity, temperature1=material_1_temp, vartemperature1=variance_1_temp, intensity2=material_2_intensity, varintensity2=variance_2_intensity, temperature2=material_2_temp, vartemperature2=variance_2_temp)

    CSV.write("out/nonlinear/data/nonlinearmc.csv", tabular)

    return nothing
end

main()
