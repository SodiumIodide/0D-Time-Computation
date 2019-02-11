#!/usr/bin/env julia

using LinearAlgebra
using DataFrames
using CSV
include("Constants.jl")

function main()::Nothing
    # Counter
    local total_iterations::Int64 = 0

    # Computational values
    local material_intensity::Vector{Float64} = zeros(num_t)
    local material_temp::Vector{Float64} = zeros(num_t)
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

    # Loop through time discretizations
    local intensity_value::Float64 = init_intensity
    local temp_value::Float64 = init_temp
    local opacity_term::Float64 = volfrac_1 * opacity_1 + volfrac_2 * opacity_2
    local spec_heat_term::Float64 = volfrac_1 * spec_heat_1 + volfrac_2 * spec_heat_2
    local dens::Float64 = volfrac_1 * dens_1 + volfrac_2 * dens_2
    for index = 1:num_t
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

            local jacobian::Array{Float64, 2} = make_jacobian(old_terms[1], old_terms[2], delta_t, opacity_term, dens, spec_heat_term)
            local func_vector::Array{Float64, 1} = [
                balance_a(old_terms[1], old_terms[2], delta_t, opacity, past_intensity),
                balance_b(old_terms[1], old_terms[2], delta_t, opacity, spec_heat, dens, past_temp)
            ]

            local delta::Vector{Float64} = jacobian \ - func_vector

            new_terms = delta + old_terms
            old_terms = new_terms

            error = relative_change(delta, original_terms)
            total_iterations += 1
        end
        intensity_value = new_terms[1]
        temp_value = new_terms[2]
        material_intensity[index] = intensity_value  # erg/cm^2-s
        material_temp[index] = temp_value  # eV
    end

    local avg_iterations::Float64 = convert(Float64, total_iterations) / convert(Float64, num_t)
    println(string("Average number of Newton iterations: ", avg_iterations))

    function plot_limits_log10(series::Vector{Float64})::Tuple{Float64, Float64}
        local new_series::Vector{Float64} = log10.(series)
        local minimum_value::Float64 = floor(minimum(new_series))
        local maximum_value::Float64 = ceil(maximum(new_series))

        return (10.0^minimum_value, 10.0^maximum_value)
    end

    tabular::DataFrame = DataFrame(time=times, intensity=material_intensity, temperature=material_temp)

    CSV.write("out/nonlinear/data/homog_nonlinear.csv", tabular)

    return nothing
end

main()
