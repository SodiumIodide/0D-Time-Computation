#!/usr/bin/env julia

include("PhysicsFunctions.jl")
using LinearAlgebra
using DataFrames
using CSV
include("Constants.jl")

function main()::Nothing
    set_zero_subnormals(true)

    # Counter
    local total_iterations::Int64 = 0

    # Computational values
    local material_intensity::Vector{Float64} = zeros(num_t)
    local material_temp::Vector{Float64} = zeros(num_t)
    local times::Vector{Float64} = @fastmath [(x * delta_t + t_init) * sol for x in 1:num_t]

    # Loop through time discretizations
    local intensity_value::Float64 = init_intensity
    local temp_value::Float64 = init_temp
    local opacity_term::Float64 = @fastmath volfrac_1 * opacity_1 + volfrac_2 * opacity_2
    local spec_heat_term::Float64 = @fastmath volfrac_1 * spec_heat_1 + volfrac_2 * spec_heat_2
    local dens::Float64 = @fastmath volfrac_1 * dens_1 + volfrac_2 * dens_2
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
        while @fastmath (error >= tolerance)
            local opacity::Float64 = @inbounds PhysicsFunctions.sigma_a(opacity_term, old_terms[2])
            local spec_heat::Float64 = @inbounds PhysicsFunctions.c_v(spec_heat_term, old_terms[2])

            local jacobian::Array{Float64, 2} = @inbounds PhysicsFunctions.make_jacobian(old_terms[1], old_terms[2], delta_t, opacity_term, dens, spec_heat_term)
            local func_vector::Array{Float64, 1} = [
                @inbounds PhysicsFunctions.balance_a(old_terms[1], old_terms[2], delta_t, opacity, past_intensity)
                @inbounds PhysicsFunctions.balance_b(old_terms[1], old_terms[2], delta_t, opacity, spec_heat, dens, past_temp)
            ]

            local delta::Vector{Float64} = @inbounds @fastmath jacobian \ - func_vector

            new_terms = @fastmath delta + old_terms
            old_terms = deepcopy(new_terms)

            error = PhysicsFunctions.relative_change(delta, original_terms)
            @fastmath total_iterations += 1
        end
        intensity_value = @inbounds new_terms[1]
        temp_value = @inbounds new_terms[2]
        @inbounds material_intensity[index] = intensity_value  # erg/cm^2-s
        @inbounds material_temp[index] = temp_value  # eV
    end

    local avg_iterations::Float64 = @fastmath convert(Float64, total_iterations) / convert(Float64, num_t)
    println("Average number of Newton iterations: ", avg_iterations)

    function plot_limits_log10(series::Vector{Float64})::Tuple{Float64, Float64}
        local new_series::Vector{Float64} = @fastmath @. log10(series)
        local minimum_value::Float64 = @fastmath floor(minimum(new_series))
        local maximum_value::Float64 = @fastmath ceil(maximum(new_series))

        return @fastmath (10.0^minimum_value, 10.0^maximum_value)
    end

    local tabular::DataFrame = DataFrame(time=times, intensity=material_intensity, temperature=material_temp)

    CSV.write("out/nonlinear/data/homog_nonlinear.csv", tabular)

    return nothing
end

main()
