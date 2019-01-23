#!/usr/bin/env julia

include("GeometryGen.jl")
include("MeshMap.jl")
include("Constants.jl")
using Random
using LinearAlgebra
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

    function balance_a(intensity::Float64, temp::Float64, delta_t::Float64, opacity::Float64, past_intensity::Float64)::Float64
        local term_1::Float64 = intensity
        local term_2::Float64 = - delta_t * sol^2 * opacity * arad * temp^4
        local term_3::Float64 = delta_t * sol * opacity * intensity
        local term_4::Float64 = - past_intensity

        return term_1 + term_2 + term_3 + term_4
    end

    function balance_b(intensity::Float64, temp::Float64, delta_t::Float64, opacity::Float64, spec_heat::Float64, density::Float64, past_temp::Float64)::Float64
        local term_1::Float64 = temp
        local term_2::Float64 = - delta_t / (density * spec_heat) * opacity * intensity
        local term_3::Float64 = delta_t / (density * spec_heat) * sol * opacity * arad * temp^4
        local term_4::Float64 = - past_temp

        return term_1 + term_2 + term_3 + term_4
    end

    function make_jacobian(intensity::Float64, temp::Float64, delta_t::Float64, opacity_term::Float64, density::Float64, spec_heat_term::Float64)::Array{Float64, 2}
        local term_1::Float64 = 1.0 + delta_t * sol * opacity_term / temp^3
        local term_2::Float64 = - delta_t * sol^2 * opacity_term * arad - 3.0 * delta_t * sol * opacity_term / temp^4 * intensity
        local term_3::Float64 = - delta_t / (density * spec_heat_term) * opacity_term / temp^3
        local term_4::Float64 = 1.0 + 3.0 * delta_t * opacity_term / (density * spec_heat_term * temp^4) * intensity + delta_t / (density * spec_heat_term) * sol * opacity_term * arad

        return [
            term_1 term_2;
            term_3 term_4
        ]
    end

    function relative_change(delta::Vector{Float64}, original_terms::Vector{Float64})::Float64
        local current_norm::Float64 = sqrt(sum([x^2 for x in delta]))
        local original_norm::Float64 = sqrt(sum([x^2 for x in original_terms]))

        return current_norm / original_norm
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

    print(string("Proceeding with ", nthreads(), " computational threads...\n"))

    # Locking conditions
    local array_lock::SpinLock = SpinLock()
    local random_lock::SpinLock = SpinLock()

    # Outer loop
    @threads for i = 1:max_iterations
        local intensity_value::Float64
        local temp_value::Float64
        local t_delta::Vector{Float64}
        local t_arr::Vector{Float64}
        local materials::Vector{Int32}
        local num_cells::Int64

        # Prevent random number clashing with lock as argument
        (t_delta, t_arr, materials, num_cells) = GeometryGen.get_geometry(chord_1, chord_2, t_max, num_divs, rng=generator, lock=random_lock)

        local intensity::Vector{Float64} = zeros(num_cells)
        local temp::Vector{Float64} = zeros(num_cells)

        # First loop uses initial conditions
        (intensity_value, temp_value) = (init_intensity, init_temp)

        # Inner loop
        for (index, material) in enumerate(materials)
            local delta_t_unstruct::Float64 = t_delta[index]
            local (opacity_term::Float64, spec_heat_term::Float64, dens::Float64) = (material == 1) ? (opacity_1, spec_heat_1, dens_1) : (opacity_2, spec_heat_2, dens_2)

            local original_terms::Vector{Float64} = [
                intensity_value,
                temp_value
            ]
            local old_terms::Vector{Float64} = deepcopy(original_terms)
            local new_terms::Vector{Float64} = deepcopy(original_terms)

            local past_intensity::Float64 = intensity_value
            local past_temp::Float64 = temp_value

            local error::Float64 = 1.0

            # Newtonian loop
            while (error >= tolerance)
                local opacity::Float64 = sigma_a(opacity_term, old_terms[2])
                local spec_heat::Float64 = c_v(spec_heat_term, old_terms[2])

                local jacobian::Array{Float64, 2} = make_jacobian(old_terms[1], old_terms[2], delta_t_unstruct, opacity_term, dens, spec_heat_term)
                local func_vector::Vector{Float64} = [
                    balance_a(old_terms[1], old_terms[2], delta_t_unstruct, opacity, past_intensity),
                    balance_b(old_terms[1], old_terms[2], delta_t_unstruct, opacity, spec_heat, dens, past_temp)
                ]

                local delta::Vector{Float64} = jacobian \ - func_vector

                new_terms = delta + old_terms
                old_terms = new_terms

                error = relative_change(delta, original_terms)
            end
            intensity_value = new_terms[1]
            temp_value = new_terms[2]
            intensity[index] = intensity_value  # erg/cm^2-s
            temp[index] = temp_value  # eV
        end

        local material_intensity_array::Array{Float64, 2} = MeshMap.material_calc(intensity, t_delta, num_cells, materials, delta_t, num_t, convert(Int32, 2))
        local material_temp_array::Array{Float64, 2} = MeshMap.material_calc(temp, t_delta, num_cells, materials, delta_t, num_t, convert(Int32, 2))
        local material_intensity_1_square::Vector{Float64} = material_intensity_array[:, 1].^2
        local material_intensity_2_square::Vector{Float64} = material_intensity_array[:, 2].^2
        local material_temp_1_square::Vector{Float64} = material_temp_array[:, 1].^2
        local material_temp_2_square::Vector{Float64} = material_temp_array[:, 2].^2

        # Lock the threaded system to increment values
        lock(array_lock) do
            material_1_intensity += material_intensity_array[:, 1]
            material_2_intensity += material_intensity_array[:, 2]
            material_1_temp += material_temp_array[:, 1]
            material_2_temp += material_temp_array[:, 2]
            sum_1_intensity_square += material_intensity_1_square
            sum_2_intensity_square += material_intensity_2_square
            sum_1_temp_square += material_temp_1_square
            sum_2_temp_square += material_temp_2_square
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

    material_1_intensity ./= (max_iterations_f * volfrac_1)
    material_2_intensity ./= (max_iterations_f * volfrac_2)
    material_1_temp ./= (max_iterations_f * volfrac_1)
    material_2_temp ./= (max_iterations_f * volfrac_2)

    tabular::DataFrame = DataFrame(time=times, intensity1=material_1_intensity, varintensity1=variance_1_intensity, temperature1=material_1_temp, vartemperature1=variance_1_temp, intensity2=material_2_intensity, varintensity2=variance_2_intensity, temperature2=material_2_temp, vartemperature2=variance_2_temp)

    CSV.write("out/nonlinear/data/nonlinear.csv", tabular)

    return nothing
end

main()
