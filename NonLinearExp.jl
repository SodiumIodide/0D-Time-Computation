#!/usr/bin/env julia

include("GeometryGen.jl")
include("MeshMap.jl")
include("Constants.jl")
using Random
using LinearAlgebra
using DataFrames
using CSV

function main()::Nothing
    # Iteration condition
    local cont_calc_outer::Bool = true
    local generator::MersenneTwister = MersenneTwister(1234)

    # Computational values
    local iteration_number::Int64 = 0
    local material_1_intensity::Vector{Float64} = zeros(num_t)
    local material_2_intensity::Vector{Float64} = zeros(num_t)
    local unconditional_intensity::Vector{Float64} = zeros(num_t)
    local material_1_temp::Vector{Float64} = zeros(num_t)
    local material_2_temp::Vector{Float64} = zeros(num_t)
    local unconditional_temp::Vector{Float64} = zeros(num_t)
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

    # Varibles for variance computation
    local sum_1_intensity_square::Vector{Float64} = zeros(num_t)
    local sum_2_intensity_square::Vector{Float64} = zeros(num_t)
    local sum_1_temp_square::Vector{Float64} = zeros(num_t)
    local sum_2_temp_square::Vector{Float64} = zeros(num_t)
    local variance_1_intensity::Vector{Float64} = zeros(num_t)
    local variance_2_intensity::Vector{Float64} = zeros(num_t)
    local variance_1_temp::Vector{Float64} = zeros(num_t)
    local variance_2_temp::Vector{Float64} = zeros(num_t)

    # Outer loop
    while (cont_calc_outer)
        local (t_delta::Vector{Float64}, t_arr::Vector{Float64}, materials::Vector{Int32}, num_cells::Int64) = GeometryGen.get_geometry(chord_1, chord_2, t_max, num_divs, rng=generator)

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

        local material_intensity_array::Array{Float64, 2} = MeshMap.material_calc(intensity, t_delta, num_cells, materials, delta_t, num_t, convert(Int32, 2))
        local material_temp_array::Array{Float64, 2} = MeshMap.material_calc(temp, t_delta, num_cells, materials, delta_t, num_t, convert(Int32, 2))
        local material_intensity_1_square::Vector{Float64} = material_intensity_array[:, 1].^2
        local material_intensity_2_square::Vector{Float64} = material_intensity_array[:, 2].^2
        local material_temp_1_square::Vector{Float64} = material_temp_array[:, 1].^2
        local material_temp_2_square::Vector{Float64} = material_temp_array[:, 2].^2

        material_1_intensity += material_intensity_array[:, 1]
        material_2_intensity += material_intensity_array[:, 2]
        material_1_temp += material_temp_array[:, 1]
        material_2_temp += material_temp_array[:, 2]
        sum_1_intensity_square += material_intensity_1_square
        sum_2_intensity_square += material_intensity_2_square
        sum_1_temp_square += material_temp_1_square
        sum_2_temp_square += material_temp_2_square

        iteration_number += 1

        if (iteration_number % num_say == 0)
            println(string("Iteration Number ", iteration_number))
        end

        if (iteration_number > max_iterations)
            cont_calc_outer = false
        end
    end

    local max_iterations_f::Float64 = convert(Float64, max_iterations)
    local variance_prefix::Float64 = 1.0 / (max_iterations_f * (max_iterations_f - 1.0))
    variance_1_intensity = variance_prefix .* (max_iterations_f .* sum_1_intensity_square - material_1_intensity.^2)
    variance_2_intensity = variance_prefix .* (max_iterations_f .* sum_2_intensity_square - material_2_intensity.^2)
    variance_1_temp = variance_prefix .* (max_iterations_f .* sum_1_temp_square - material_1_temp.^2)
    variance_2_temp = variance_prefix .* (max_iterations_f .* sum_2_temp_square - material_2_temp.^2)

    material_1_intensity ./= (max_iterations_f * volfrac_1)
    material_2_intensity ./= (max_iterations_f * volfrac_2)
    material_1_temp ./= (max_iterations_f * volfrac_1)
    material_2_temp ./= (max_iterations_f * volfrac_2)

    tabular::DataFrame = DataFrame(time=times, intensity1=material_1_intensity, varintensity1=variance_1_intensity, temperature1=material_1_temp, vartemperature1=variance_1_temp, intensity2=material_2_intensity, varintensity2=variance_2_intensity, temperature2=material_2_temp, vartemperature2=variance_2_temp)

    CSV.write("out/nonlinear/data/nonlinear.csv", tabular)

    return nothing
end

main()
