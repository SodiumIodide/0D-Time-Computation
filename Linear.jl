#!/usr/bin/env julia

# Linear.jl

include("Constants.jl")
using LinearAlgebra
using DataFrames
using CSV

function main()::Nothing
    # Initial conditions
    local init_intensity_1::Float64 = init_intensity  # erg/cm^2-s
    local init_intensity_2::Float64 = init_intensity  # erg/cm^2-s
    local init_energy_1::Float64 = arad * init_temp^4  # erg/cm^3
    local init_energy_2::Float64 = arad * init_temp^4  # erg/cm^3

    # Compute constant terms for matrix
    local term_1::Float64 = 1.0 + delta_t * sol * opacity_1 + delta_t / chord_1
    local term_2::Float64 = - delta_t / chord_2
    local term_3::Float64 = - delta_t * sol^2 * opacity_1
    local term_4::Float64 = - delta_t / chord_1
    local term_5::Float64 = 1.0 + delta_t * sol * opacity_2 + delta_t / chord_2
    local term_6::Float64 = - delta_t * sol^2 * opacity_2
    local term_7::Float64 = - delta_t * factor_1 * opacity_1
    local term_8::Float64 = 1.0 + delta_t * factor_1 * opacity_1 * sol + delta_t / chord_1
    local term_9::Float64 = - delta_t / chord_2
    local term_10::Float64 = - delta_t * factor_2 * opacity_2
    local term_11::Float64 = - delta_t / chord_1
    local term_12::Float64 = 1.0 + delta_t * factor_2 * opacity_2 * sol + delta_t / chord_2

    # Set up problem definitions
    # This matrix is a constant parameter
    local coupled_matrix::Array{Float64,2} = [
        term_1 term_2 term_3 0.0;
        term_4 term_5 0.0 term_6;
        term_7 0.0 term_8 term_9;
        0.0 term_10 term_11 term_12
    ]

    # Vectors in problem
    local vector::Vector{Float64} = Vector{Float64}(undef, 4)
    local previous_vector::Array{Float64,1} = [
        *(init_intensity_1, volfrac_1)
        *(init_intensity_2, volfrac_2)
        *(init_energy_1, volfrac_1)
        *(init_energy_2, volfrac_2)
    ]

    # Preallocate vectors
    local intensity_1::Vector{Float64} = Vector{Float64}(undef, num_t)
    local intensity_2::Vector{Float64} = Vector{Float64}(undef, num_t)
    local energy_1::Vector{Float64} = Vector{Float64}(undef, num_t)
    local energy_2::Vector{Float64} = Vector{Float64}(undef, num_t)

    for t = 1:num_t
        vector = coupled_matrix \ previous_vector
        previous_vector = vector
        intensity_1[t] = previous_vector[1]
        intensity_2[t] = previous_vector[2]
        energy_1[t] = previous_vector[3]
        energy_2[t] = previous_vector[4]
    end

    local times::Vector{Float64} = [(x * delta_t + t_init) * sol for x in 1:num_t]

    # Take the implicit volume fraction into account
    intensity_1 /= volfrac_1  # erg/cm^2-s
    intensity_2 /= volfrac_2  # erg/cm^2-s
    energy_1 /= volfrac_1  # erg/cm^3
    energy_2 /= volfrac_2  # erg/cm^3

    local temp_1::Vector{Float64} = (energy_1 ./ arad).^(1.0 / 4.0)  # eV
    local temp_2::Vector{Float64} = (energy_2 ./ arad).^(1.0 / 4.0)  # eV

    local y_min::Float64 = floor(log10(min(min(intensity_1...), min(intensity_2...))))
    local y_max::Float64 = ceil(log10(max(max(intensity_1...), max(intensity_2...))))

    tabular::DataFrame = DataFrame(time=times, intensity1=intensity_1, temperature1=temp_1, intensity2=intensity_2, temperature2=temp_2)

    CSV.write("out/nonlinear/data/linear.csv", tabular)

    return nothing
end

main()
