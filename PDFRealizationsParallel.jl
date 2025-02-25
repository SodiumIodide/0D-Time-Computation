#!/usr/bin/env julia

include("ExponentialHist.jl")
include("GeometryGen.jl")
include("MeshMap.jl")
include("PDFFunctions.jl")
include("PhysicsFunctions.jl")
using Random
using Future
using LinearAlgebra
using DataFrames
using CSV
using Base.Threads
include("Constants.jl")

function main()::Nothing
    set_zero_subnormals(true)

    # Read previous data to explicitly create binned elements
    # The Monte Carlo data is read, because it is the faster computation method
    # (in comparison to realization generation)
    local mc_data_location::String = "out/nonlinear/data/nonlinearmc.csv"

    # User Input for data existence
    println("\n")
    println("Did you remember to run the Monte Carlo computation for a two-pass binning?")
    println("Press Enter to continue (no additional input)")
    println("Press C and Enter to run the Monte Carlo computation using current values in Constants.jl")
    println("Press K and Enter to quit")
    local continue_calc::String = lowercase(chomp(readline()))
    if (occursin("k", continue_calc))
        println("Terminating execution")
        return nothing
    elseif (occursin("c", continue_calc))
        println("Running Monte Carlo computation now! Please wait...")
        local cmd::Cmd = Base.julia_cmd()
        run(`$cmd NonLinearMC.jl`)
        println("Finished")
    end

    # Datafile
    local old_data::DataFrame

    if (isfile(mc_data_location))
        old_data = CSV.File(mc_data_location) |> DataFrame
    else
        println("Cannot find results of a Monte Carlo run, terminating")
        return nothing
    end

    local intensity_1_ss::Int64 = @inbounds PDFFunctions.locate_steady_state(vec([old_data.intensity1...]))
    local intensity_2_ss::Int64 = @inbounds PDFFunctions.locate_steady_state(vec([old_data.intensity2...]))
    local temperature_1_ss::Int64 = @inbounds PDFFunctions.locate_steady_state(vec([old_data.temperature1...]))
    local temperature_2_ss::Int64 = @inbounds PDFFunctions.locate_steady_state(vec([old_data.temperature2...]))

    local steady_state_index::Int64 = @fastmath max(intensity_1_ss, intensity_2_ss, temperature_1_ss, temperature_2_ss)
    local steady_state_time::Float64 = old_data.time[steady_state_index]

    println("Steady state index found to occur at point ", steady_state_index, ", time=", steady_state_time, " ct")

    local chosenindex::Int64
    local chosentime::Float64
    println("Select data point to create histograms...")
    println("[1] Steady State -DEFAULT-")
    println("[2] Early Time-Step")
    local timestepselect::String = chomp(readline())
    if (occursin("2", timestepselect))
        chosenindex = hist_early_index
        chosentime = old_data.time[chosenindex]
    else
        chosenindex = steady_state_index
        chosentime = steady_state_time
    end

    # Time (constant array in csv)
    local time_array::Vector{Float64} = fill(chosentime, num_bins + 1)

    # Bounds for binning
    local intensity_1_min::Float64 = old_data.minintensity1[chosenindex]
    local intensity_1_max::Float64 = old_data.maxintensity1[chosenindex]
    local temperature_1_min::Float64 = old_data.mintemperature1[chosenindex]
    local temperature_1_max::Float64 = old_data.maxtemperature1[chosenindex]
    local intensity_2_min::Float64 = old_data.minintensity2[chosenindex]
    local intensity_2_max::Float64 = old_data.maxintensity2[chosenindex]
    local temperature_2_min::Float64 = old_data.mintemperature2[chosenindex]
    local temperature_2_max::Float64 = old_data.maxtemperature2[chosenindex]

    local opacity_1_op::Vector{Float64} = [
        PhysicsFunctions.sigma_a(opacity_1, temperature_1_min);
        PhysicsFunctions.sigma_a(opacity_1, temperature_1_max)
    ]
    local opacity_2_op::Vector{Float64} = [
        PhysicsFunctions.sigma_a(opacity_2, temperature_2_min);
        PhysicsFunctions.sigma_a(opacity_2, temperature_2_max)
    ]
    local opacity_1_min::Float64 = @fastmath minimum(opacity_1_op)
    local opacity_1_max::Float64 = @fastmath maximum(opacity_1_op)
    local opacity_2_min::Float64 = @fastmath minimum(opacity_2_op)
    local opacity_2_max::Float64 = @fastmath maximum(opacity_2_op)

    # Arrays for binning
    local intensity_1_bin::Vector{ExponentialHist.ExpHist} = [ExponentialHist.ExpHist(num_bins, intensity_1_min, intensity_1_max) for i in 1:nthreads()]
    local intensity_2_bin::Vector{ExponentialHist.ExpHist} = [ExponentialHist.ExpHist(num_bins, intensity_2_min, intensity_2_max) for i in 1:nthreads()]
    local temperature_1_bin::Vector{ExponentialHist.ExpHist} = [ExponentialHist.ExpHist(num_bins, temperature_1_min, temperature_1_max) for i in 1:nthreads()]
    local temperature_2_bin::Vector{ExponentialHist.ExpHist} = [ExponentialHist.ExpHist(num_bins, temperature_2_min, temperature_2_max) for i in 1:nthreads()]
    local opacity_1_bin::Vector{ExponentialHist.ExpHist} = [ExponentialHist.ExpHist(num_bins, opacity_1_min, opacity_1_max) for i in 1:nthreads()]
    local opacity_2_bin::Vector{ExponentialHist.ExpHist} = [ExponentialHist.ExpHist(num_bins, opacity_2_min, opacity_2_max) for i in 1:nthreads()]

    # BEGIN REALIZATIONS SECOND PASS

    # Iteration condition
    local gen_array::Array{MersenneTwister, 1} = let m::MersenneTwister = MersenneTwister(1234)
        @fastmath [m; accumulate(Future.randjump, fill(big(10)^20, nthreads() - 1), init=m)]
    end

    # Computational values
    local new_end_time::Float64 = @fastmath chosentime / sol
    local new_delta_t::Float64 = @fastmath new_end_time / num_t_hist
    local times::Vector{Float64} = @fastmath [(x * new_delta_t + t_init) * sol for x in 1:num_t_hist]

    println("Proceeding with ", nthreads(), " computational threads...")

    local printlock::SpinLock = SpinLock()

    # Outer loop
    @threads for i = 1:max_iterations_hist
        local (t_delta::Vector{Float64}, t_arr::Vector{Float64}, materials::Vector{Int32}, num_cells::Int64) = @inbounds GeometryGen.get_geometry(chord_1, chord_2, new_end_time, num_divs_hist, rng=gen_array[threadid()])

        local intensity::Vector{Float64} = zeros(num_cells)
        local temp::Vector{Float64} = zeros(num_cells)

        # First loop uses initial conditions
        local (intensity_value::Float64, temp_value::Float64) = (init_intensity, init_temp)
        local material_num::Int32 = @inbounds materials[1]

        # Inner loop
        for index = 1:num_cells
            material_num = @inbounds materials[index]
            local delta_t_unstruct::Float64 = @inbounds t_delta[index]
            local (opacity_term::Float64, spec_heat_term::Float64, dens_term::Float64) = @fastmath (material_num == 1) ? (opacity_1, spec_heat_1, dens_1) : (opacity_2, spec_heat_2, dens_2)

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
            intensity_value = @inbounds new_terms[1]  # erg/cm^2-s
            temp_value = @inbounds new_terms[2]  # eV

            if @fastmath(material_num == 1)
                @inbounds ExponentialHist.push(intensity_1_bin[threadid()], intensity_value)
                @inbounds ExponentialHist.push(temperature_1_bin[threadid()], temp_value)
                @inbounds ExponentialHist.push(opacity_1_bin[threadid()], PhysicsFunctions.sigma_a(opacity_1, temp_value))
            else
                @inbounds ExponentialHist.push(intensity_2_bin[threadid()], intensity_value)
                @inbounds ExponentialHist.push(temperature_2_bin[threadid()], temp_value)
                @inbounds ExponentialHist.push(opacity_2_bin[threadid()], PhysicsFunctions.sigma_a(opacity_2, temp_value))
            end
        end

        # Need to reference Core namespace for thread-safe printing
        if @fastmath(i % num_say == 0)
            lock(printlock) do
                Core.println("Iteration Number ", i)
            end
        end
    end

    local material_1_intensity_bin::Vector{Float64} = zeros(num_bins + 1)
    local material_2_intensity_bin::Vector{Float64} = zeros(num_bins + 1)
    local material_1_temperature_bin::Vector{Float64} = zeros(num_bins + 1)
    local material_2_temperature_bin::Vector{Float64} = zeros(num_bins + 1)
    local material_1_opacity_bin::Vector{Float64} = zeros(num_bins + 1)
    local material_2_opacity_bin::Vector{Float64} = zeros(num_bins + 1)
    @simd for i in 1:nthreads()
        @inbounds @fastmath material_1_intensity_bin += ExponentialHist.histogram(intensity_1_bin[i])
        @inbounds @fastmath material_2_intensity_bin += ExponentialHist.histogram(intensity_2_bin[i])
        @inbounds @fastmath material_1_temperature_bin += ExponentialHist.histogram(temperature_1_bin[i])
        @inbounds @fastmath material_2_temperature_bin += ExponentialHist.histogram(temperature_2_bin[i])
        @inbounds @fastmath material_1_opacity_bin += ExponentialHist.histogram(opacity_1_bin[i])
        @inbounds @fastmath material_2_opacity_bin += ExponentialHist.histogram(opacity_2_bin[i])
    end

    # Normalize histograms
    @fastmath material_1_intensity_bin ./= sum(material_1_intensity_bin)
    @fastmath material_2_intensity_bin ./= sum(material_2_intensity_bin)
    @fastmath material_1_temperature_bin ./= sum(material_1_temperature_bin)
    @fastmath material_2_temperature_bin ./= sum(material_2_temperature_bin)
    @fastmath material_1_opacity_bin ./= sum(material_1_opacity_bin)
    @fastmath material_2_opacity_bin ./= sum(material_2_opacity_bin)

    # The distribution can be computed from the first array point
    local material_1_intensity_array::Vector{Float64} = @inbounds ExponentialHist.distribution(intensity_1_bin[1])
    local material_2_intensity_array::Vector{Float64} = @inbounds ExponentialHist.distribution(intensity_2_bin[1])
    local material_1_temperature_array::Vector{Float64} = @inbounds ExponentialHist.distribution(temperature_1_bin[1])
    local material_2_temperature_array::Vector{Float64} = @inbounds ExponentialHist.distribution(temperature_2_bin[1])
    local material_1_opacity_array::Vector{Float64} = @inbounds ExponentialHist.distribution(opacity_1_bin[1])
    local material_2_opacity_array::Vector{Float64} = @inbounds ExponentialHist.distribution(opacity_2_bin[1])

    local tabular::DataFrame = DataFrame(time=time_array, intensity1arr=material_1_intensity_array, freqintensity1=material_1_intensity_bin, intensity2arr=material_2_intensity_array, freqintensity2=material_2_intensity_bin, temperature1arr=material_1_temperature_array, freqtemperature1=material_1_temperature_bin, temperature2arr=material_2_temperature_array, freqtemperature2=material_2_temperature_bin, opacity1arr=material_1_opacity_array, freqopacity1=material_1_opacity_bin, opacity2arr=material_2_opacity_array, freqopacity2=material_2_opacity_bin)

    CSV.write("out/nonlinear/pdf_data/realizations_pdf.csv", tabular)

    return nothing
end

main()
