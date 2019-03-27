#!/usr/bin/env julia

include("GeometryGen.jl")
include("PhysicsFunctions.jl")
include("MeshMap.jl")
include("RunningStatistics.jl")
using Random
using LinearAlgebra
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

    # Arrays for online computation
    local stat_1_intensity::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:num_t]
    local stat_2_intensity::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:num_t]
    local stat_1_temp::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:num_t]
    local stat_2_temp::Vector{RunningStatistics.RunningStat} = [RunningStatistics.RunningStat() for i in 1:num_t]

    # Outer loop
    while (cont_calc_outer)
        local (t_delta::Vector{Float64}, t_arr::Vector{Float64}, materials::Vector{Int32}, num_cells::Int64) = GeometryGen.get_geometry(chord_1, chord_2, t_max, num_divs, rng=generator)

        local intensity::Vector{Float64} = zeros(num_cells)
        local temp::Vector{Float64} = zeros(num_cells)

        # First loop uses initial conditions
        local (intensity_value::Float64, temp_value::Float64) = (init_intensity, init_temp)

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

            local error::Float64 = 1.0

            # Newtonian loop
            while (error >= tolerance)
                local opacity::Float64 = PhysicsFunctions.sigma_a(opacity_term, old_terms[2])
                local spec_heat::Float64 = PhysicsFunctions.c_v(spec_heat_term, old_terms[2])

                local jacobian::Array{Float64, 2} = PhysicsFunctions.make_jacobian(old_terms[1], old_terms[2], delta_t_unstruct, opacity_term, dens, spec_heat_term)
                local func_vector::Vector{Float64} = [
                    PhysicsFunctions.balance_a(old_terms[1], old_terms[2], delta_t_unstruct, opacity, intensity_value),
                    PhysicsFunctions.balance_b(old_terms[1], old_terms[2], delta_t_unstruct, opacity, spec_heat, dens, temp_value)
                ]

                local delta::Vector{Float64} = @inbounds @fastmath jacobian \ - func_vector

                @inbounds new_terms = @fastmath delta + old_terms
                old_terms = new_terms

                error = PhysicsFunctions.relative_change(delta, original_terms)
            end
            intensity_value = new_terms[1]
            temp_value = new_terms[2]
            intensity[index] = intensity_value  # erg/cm^2-s
            temp[index] = temp_value  # eV
        end

        local material_intensity_array::Array{Float64, 2} = MeshMap.material_calc(intensity, t_delta, num_cells, materials, delta_t, num_t, convert(Int32, 2))
        local material_temp_array::Array{Float64, 2} = MeshMap.material_calc(temp, t_delta, num_cells, materials, delta_t, num_t, convert(Int32, 2))

        for k in 1:num_t
            if (material_intensity_array[k, 1] != 0.0)
                RunningStatistics.push(stat_1_intensity[k], material_intensity_array[k, 1])  # erg/cm^2-s
                RunningStatistics.push(stat_1_temp[k], material_temp_array[k, 1])  # eV
            else
                RunningStatistics.push(stat_2_intensity[k], material_intensity_array[k, 2])  # erg/cm^2-s
                RunningStatistics.push(stat_2_temp[k], material_temp_array[k, 2])  # eV
            end
        end

        iteration_number += 1

        if (iteration_number % num_say == 0)
            println(string("Realization Number ", iteration_number))
        end

        if (iteration_number > max_iterations)
            cont_calc_outer = false
        end
    end

    local material_1_intensity::Vector{Float64} = RunningStatistics.mean.(stat_1_intensity)
    local material_2_intensity::Vector{Float64} = RunningStatistics.mean.(stat_2_intensity)
    local material_1_temp::Vector{Float64} = RunningStatistics.mean.(stat_1_temp)
    local material_2_temp::Vector{Float64} = RunningStatistics.mean.(stat_2_temp)
    local variance_1_intensity::Vector{Float64} = RunningStatistics.variance.(stat_1_intensity)
    local variance_2_intensity::Vector{Float64} = RunningStatistics.variance.(stat_2_intensity)
    local variance_1_temp::Vector{Float64} = RunningStatistics.variance.(stat_1_temp)
    local variance_2_temp::Vector{Float64} = RunningStatistics.variance.(stat_2_temp)

    local tabular::DataFrame = DataFrame(time=times, intensity1=material_1_intensity, varintensity1=variance_1_intensity, temperature1=material_1_temp, vartemperature1=variance_1_temp, intensity2=material_2_intensity, varintensity2=variance_2_intensity, temperature2=material_2_temp, vartemperature2=variance_2_temp)

    CSV.write("out/nonlinear/data/nonlinear.csv", tabular)

    return nothing
end

main()
