#!/usr/bin/env julia

include("GeometryGen.jl")
include("PhysicsFunctions.jl")
include("MeshMap.jl")
include("RunningStatistics.jl")
using Random
using Future
using LinearAlgebra
using DataFrames
using CSV
using Base.Threads
include("Constants.jl")

function main()::Nothing
    set_zero_subnormals(true)
    # Iteration condition
    local gen_array::Array{MersenneTwister, 1} = let m::MersenneTwister = MersenneTwister(1234)
        @fastmath [m; accumulate(Future.randjump, fill(big(10)^20, nthreads() - 1), init=m)]
    end

    # Computational values
    local times::Vector{Float64} = @fastmath [(x * delta_t + t_init) * sol for x in 1:num_t]

    # Parallel arrays
    local stat_1_intensity::Array{RunningStatistics.RunningStat, 2} = RunningStatistics.threadarray(num_t)
    local stat_2_intensity::Array{RunningStatistics.RunningStat, 2} = RunningStatistics.threadarray(num_t)
    local stat_1_temp::Array{RunningStatistics.RunningStat, 2} = RunningStatistics.threadarray(num_t)
    local stat_2_temp::Array{RunningStatistics.RunningStat, 2} = RunningStatistics.threadarray(num_t)

    println("Proceeding with ", nthreads(), " computational threads...")

    local printlock::SpinLock = SpinLock()

    # Outer loop
    @threads for i = 1:max_iterations
        # Prevent random number clashing with discrete generators
        local (t_delta::Vector{Float64}, t_arr::Vector{Float64}, materials::Vector{Int32}, num_cells::Int64) = @inbounds GeometryGen.get_geometry(chord_1, chord_2, t_max, num_divs, rng=gen_array[threadid()])

        local intensity::Vector{Float64} = zeros(num_cells)
        local temp::Vector{Float64} = zeros(num_cells)

        # First loop uses initial conditions
        local (intensity_value::Float64, temp_value::Float64) = (init_intensity, init_temp)

        # Inner loop
        for (index, material) in enumerate(materials)
            local delta_t_unstruct::Float64 = @inbounds t_delta[index]
            local (opacity_term::Float64, spec_heat_term::Float64, dens_term::Float64) = @fastmath (material == 1) ? (opacity_1, spec_heat_1, dens_1) : (opacity_2, spec_heat_2, dens_2)

            local original_terms::Vector{Float64} = [
                intensity_value,
                temp_value
            ]
            local old_terms::Vector{Float64} = deepcopy(original_terms)
            local new_terms::Vector{Float64} = deepcopy(original_terms)

            local error::Float64 = 1.0

            # Newtonian loop
            while @fastmath(error >= tolerance)
                local opacity::Float64 = @inbounds PhysicsFunctions.sigma_a(opacity_term, old_terms[2])
                local spec_heat::Float64 = @inbounds PhysicsFunctions.c_v(spec_heat_term, old_terms[2])
                local dens::Float64 = @inbounds PhysicsFunctions.rho(dens_term, old_terms[2])

                local jacobian::Array{Float64, 2} = @inbounds PhysicsFunctions.make_jacobian(old_terms[1], old_terms[2], delta_t_unstruct, opacity_term, dens, spec_heat_term)
                local func_vector::Vector{Float64} = [
                    @inbounds PhysicsFunctions.balance_a(old_terms[1], old_terms[2], delta_t_unstruct, opacity, intensity_value)
                    @inbounds PhysicsFunctions.balance_b(old_terms[1], old_terms[2], delta_t_unstruct, opacity, spec_heat, dens, temp_value)
                ]

                local delta::Vector{Float64} = @inbounds @fastmath jacobian \ - func_vector

                new_terms = @fastmath delta + old_terms
                old_terms = new_terms

                error = PhysicsFunctions.relative_change(delta, original_terms)
            end
            intensity_value = @inbounds new_terms[1]
            temp_value = @inbounds new_terms[2]
            @inbounds intensity[index] = intensity_value  # erg/cm^2-s
            @inbounds temp[index] = temp_value  # eV
        end

        local material_intensity_array::Array{Float64, 2} = MeshMap.material_calc(intensity, t_delta, num_cells, materials, delta_t, num_t, convert(Int32, 2))
        local material_temp_array::Array{Float64, 2} = MeshMap.material_calc(temp, t_delta, num_cells, materials, delta_t, num_t, convert(Int32, 2))

        @simd for k in 1:num_t
            if @inbounds @fastmath(material_intensity_array[k, 1] != 0.0)
                @inbounds RunningStatistics.push(stat_1_intensity[k, threadid()], material_intensity_array[k, 1])  # erg/cm^2-s
                @inbounds RunningStatistics.push(stat_1_temp[k, threadid()], material_temp_array[k, 1])  # eV
            else
                @inbounds RunningStatistics.push(stat_2_intensity[k, threadid()], material_intensity_array[k, 2])  # erg/cm^2-s
                @inbounds RunningStatistics.push(stat_2_temp[k, threadid()], material_temp_array[k, 2])  # eV
            end
        end

        # Need to reference Core namespace for thread-safe printing
        if @fastmath(i % num_say == 0)
            lock(printlock) do
                Core.println("Realization Number ", i)
            end
        end
    end

    # Total tallies for each quantity
    local num_1::Vector{Float64} = RunningStatistics.total(stat_1_intensity)
    local num_2::Vector{Float64} = RunningStatistics.total(stat_2_intensity)

    # Mean values
    local material_1_intensity::Vector{Float64} = RunningStatistics.compute_mean(stat_1_intensity, num_1)
    local material_2_intensity::Vector{Float64} = RunningStatistics.compute_mean(stat_2_intensity, num_2)
    local material_1_temp::Vector{Float64} = RunningStatistics.compute_mean(stat_1_temp, num_1)
    local material_2_temp::Vector{Float64} = RunningStatistics.compute_mean(stat_2_temp, num_2)

    # Variance values
    local variance_1_intensity::Vector{Float64} = RunningStatistics.compute_variance(stat_1_intensity, num_1, material_1_intensity)
    local variance_2_intensity::Vector{Float64} = RunningStatistics.compute_variance(stat_2_intensity, num_2, material_2_intensity)
    local variance_1_temp::Vector{Float64} = RunningStatistics.compute_variance(stat_1_temp, num_1, material_1_temp)
    local variance_2_temp::Vector{Float64} = RunningStatistics.compute_variance(stat_2_temp, num_2, material_2_temp)

    # Maximum and minimum values
    local max_intensity_1::Vector{Float64} = RunningStatistics.compute_max(stat_1_intensity)
    local min_intensity_1::Vector{Float64} = RunningStatistics.compute_min(stat_1_intensity)
    local max_intensity_2::Vector{Float64} = RunningStatistics.compute_max(stat_2_intensity)
    local min_intensity_2::Vector{Float64} = RunningStatistics.compute_min(stat_2_intensity)
    local max_temp_1::Vector{Float64} = RunningStatistics.compute_max(stat_1_temp)
    local min_temp_1::Vector{Float64} = RunningStatistics.compute_min(stat_1_temp)
    local max_temp_2::Vector{Float64} = RunningStatistics.compute_max(stat_2_temp)
    local min_temp_2::Vector{Float64} = RunningStatistics.compute_min(stat_2_temp)

    local tabular::DataFrame = DataFrame(time=times, intensity1=material_1_intensity, varintensity1=variance_1_intensity, maxintensity1=max_intensity_1, minintensity1=min_intensity_1, temperature1=material_1_temp, vartemperature1=variance_1_temp, maxtemperature1=max_temp_1, mintemperature1=min_temp_1, intensity2=material_2_intensity, varintensity2=variance_2_intensity, maxintensity2=max_intensity_2, minintensity2=min_intensity_2, temperature2=material_2_temp, vartemperature2=variance_2_temp, maxtemperature2=max_temp_2, mintemperature2=min_temp_2)

    CSV.write("out/nonlinear/data/nonlinear.csv", tabular)

    return nothing
end

main()
