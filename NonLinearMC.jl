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
    local material_1_temp::Vector{Float64} = zeros(num_t)
    local material_2_temp::Vector{Float64} = zeros(num_t)
    local times::Vector{Float64} = [(x * delta_t + t_init) * sol for x in 1:num_t]

    # Probability for material sampling
    local prob_1::Float64 = chord_1 / (chord_1 + chord_2)

    function sigma_a(opacity_term::Float64, temp::Float64)::Float64
        return opacity_term / temp^3
    end

    function c_v(spec_heat_term::Float64, temp::Float64)::Float64
        return spec_heat_term
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
        local rand_num::Float64
        local material_num::Int32

        # First loop uses initial conditions
        local (intensity_value::Float64, temp_value::Float64) = (init_intensity, init_temp)

        # Sample initial starting material
        lock(random_lock) do
            rand_num = rand(generator, Float64)
        end
        material_num = (rand_num < prob_1) ? 1 : 2

        # Inner loop
        for (index, time) in enumerate(times)
            local (opacity_term::Float64, spec_heat_term::Float64, dens::Float64) = (material == 1) ? (opacity_1, spec_heat_1, dens_1) : (opacity_2, spec_heat_2, dens_2)
        end

        atomic_add!(iteration_number, 1)

        if (iteration_number[] % num_say == 0)
            println(string("Iteration Number ", iteration_number[]))
        end
    end

    return nothing
end

main()
