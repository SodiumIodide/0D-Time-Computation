#!/usr/bin/env julia

include("PhysicsFunctions.jl")
include("RunningStatistics.jl")
using LinearAlgebra
using DataFrames
using CSV
using ProgressMeter
include("Constants.jl")

function main()::Nothing
    set_zero_subnormals(true)

    # Computational values
    local times::Vector{Float64} = @fastmath [(x * delta_t + t_init) * sol for x in 1:num_t]

    # Material property arrays
    local intensity_1::Vector{Float64} = zeros(num_t)
    local intensity_2::Vector{Float64} = zeros(num_t)
    local temperature_1::Vector{Float64} = zeros(num_t)
    local temperature_2::Vector{Float64} = zeros(num_t)

    # Begin with initial conditions
    local intensity_value_1::Float64 = init_intensity
    local intensity_value_2::Float64 = init_intensity
    local temp_value_1::Float64 = init_temp
    local temp_value_2::Float64 = init_temp

    @showprogress 1 for index=1:@fastmath(num_t)
        local original_terms::Vector{Float64} = vec([
            intensity_value_1
            intensity_value_2
            temp_value_1
            temp_value_2
        ])

        local old_terms::Vector{Float64} = deepcopy(original_terms)
        local new_terms::Vector{Float64} = zeros(4)

        local error::Float64 = 1.0

        # Newtonian loop
        while @fastmath (error >= heur_tolerance)
            local sigma_a_1::Float64 = @inbounds PhysicsFunctions.sigma_a(opacity_1, old_terms[3])
            local sigma_a_2::Float64 = @inbounds PhysicsFunctions.sigma_a(opacity_2, old_terms[4])
            local c_v_1::Float64 = @inbounds PhysicsFunctions.c_v(spec_heat_1, old_terms[3])
            local c_v_2::Float64 = @inbounds PhysicsFunctions.c_v(spec_heat_2, old_terms[4])
            local rho_1::Float64 = @inbounds PhysicsFunctions.rho(dens_1, old_terms[3])
            local rho_2::Float64 = @inbounds PhysicsFunctions.rho(dens_2, old_terms[4])

            #local jacobian::Array{Float64, 2} = PhysicsFunctions.make_jacobian_heuristic(old_terms...)
            local jacobian::Array{Float64, 2} = PhysicsFunctions.complex_step_jacobian_heuristic(old_terms...)
            local func_vector::Vector{Float64} = vec([
                @inbounds PhysicsFunctions.balance_f1(old_terms..., sigma_a_1, sigma_a_2, c_v_1, c_v_2, rho_1, rho_2, intensity_value_1)
                @inbounds PhysicsFunctions.balance_f2(old_terms..., sigma_a_1, sigma_a_2, c_v_1, c_v_2, rho_1, rho_2, intensity_value_2)
                @inbounds PhysicsFunctions.balance_g1(old_terms..., sigma_a_1, sigma_a_2, c_v_1, c_v_2, rho_1, rho_2, temp_value_1)
                @inbounds PhysicsFunctions.balance_g2(old_terms..., sigma_a_1, sigma_a_2, c_v_1, c_v_2, rho_1, rho_2, temp_value_2)
            ])

            local delta::Vector{Float64} = @inbounds @fastmath jacobian \ - func_vector
            @inbounds new_terms = @fastmath delta + old_terms
            old_terms = deepcopy(new_terms)

            error = PhysicsFunctions.relative_change(delta, original_terms)
        end
        intensity_value_1 = @inbounds new_terms[1]
        intensity_value_2 = @inbounds new_terms[2]
        temp_value_1 = @inbounds new_terms[3]
        temp_value_2 = @inbounds new_terms[4]
        @inbounds intensity_1[index] = intensity_value_1  # erg/cm^2-s
        @inbounds intensity_2[index] = intensity_value_2  # erg/cm^2-s
        @inbounds temperature_1[index] = temp_value_1  # eV
        @inbounds temperature_2[index] = temp_value_2  # eV
    end

    local tabular::DataFrame = DataFrame(time=times, intensity1=intensity_1, temperature1=temperature_1, intensity2=intensity_2, temperature2=temperature_2)

    CSV.write("out/nonlinear/data/heuristic.csv", tabular)

    return nothing
end

main()
