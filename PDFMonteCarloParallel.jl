#!/usr/bin/env julia

include("ExponentialHist.jl")
include("PhysicsFunctions.jl")
using Random
using DataFrames
using CSV
using Base.Threads
using Future
include("Constants.jl")

function main()::Nothing
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

    local intensity_1_ss::Int64 = PDFFunctions.locate_steady_state(vec([old_data.intensity1...]))
    local intensity_2_ss::Int64 = PDFFunctions.locate_steady_state(vec([old_data.intensity2...]))
    local temperature_1_ss::Int64 = PDFFunctions.locate_steady_state(vec([old_data.temperature1...]))
    local temperature_2_ss::Int64 = PDFFunctions.locate_steady_state(vec([old_data.temperature2...]))

    local steady_state_index::Int64 = max(intensity_1_ss, intensity_2_ss, temperature_1_ss, temperature_2_ss)
    local steady_state_time::Float64 = old_data.time[steady_state_index]

    println("Steady state index found to occur at point ", steady_state_index, ", time=", steady_state_time, " ct")

    local chosenindex::Int64
    println("Select data point to create histograms...")
    println("[1] Steady State -DEFAULT-")
    println("[2] Early Time-Step")
    local timestepselect::String = chomp(readline())
    if (occursin("2", timestepselect))
        chosenindex = 3
    else
        chosenindex = steady_state_index
    end

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
    local opacity_1_min::Float64 = minimum(opacity_1_op)
    local opacity_1_max::Float64 = maximum(opacity_1_op)
    local opacity_2_min::Float64 = minimum(opacity_2_op)
    local opacity_2_max::Float64 = maximum(opacity_2_op)

    # Arrays for binning
    local intensity_1_bin::Vector{ExponentialHist.ExpHist} = [ExponentialHist.ExpHist(num_bins, intensity_1_min, intensity_1_max) for i in 1:nthreads()]
    local intensity_2_bin::Vector{ExponentialHist.ExpHist} = [ExponentialHist.ExpHist(num_bins, intensity_2_min, intensity_2_max) for i in 1:nthreads()]
    local temperature_1_bin::Vector{ExponentialHist.ExpHist} = [ExponentialHist.ExpHist(num_bins, temperature_1_min, temperature_1_max) for i in 1:nthreads()]
    local temperature_2_bin::Vector{ExponentialHist.ExpHist} = [ExponentialHist.ExpHist(num_bins, temperature_2_min, temperature_2_max) for i in 1:nthreads()]
    local opacity_1_bin::Vector{ExponentialHist.ExpHist} = [ExponentialHist.ExpHist(num_bins, opacity_1_min, opacity_1_max) for i in 1:nthreads()]
    local opacity_2_bin::Vector{ExponentialHist.ExpHist} = [ExponentialHist.ExpHist(num_bins, opacity_2_min, opacity_2_max) for i in 1:nthreads()]

    # BEGIN MONTE CARLO SECOND PASS

    # Iteration condition
    local gen_array::Array{MersenneTwister, 1} = let m::MersenneTwister = MersenneTwister(1234)
        [m; accumulate(Future.randjump, fill(big(10)^20, nthreads() - 1), init=m)]
    end

    # Computational values
    local new_delta_t::Float64 = (steady_state_time / sol) / num_t_hist
    local times::Vector{Float64} = [(x * new_delta_t + t_init) * sol for x in 1:num_t_hist]

    # Probability for material sampling
    local prob_1::Float64 = chord_1 / (chord_1 + chord_2)
    local change_prob_1::Float64 = 1.0 / chord_1 * delta_t
    local change_prob_2::Float64 = 1.0 / chord_2 * delta_t
    if ((change_prob_1 > 1.0) || (change_prob_2 > 1.0))
        println("The value for delta_t is too large for sampling")
        return nothing
    end

    println("Proceeding with ", nthreads(), " computational threads...")

    local printlock::SpinLock = SpinLock()

    # Outer loop
    @threads for i = 1:max_iterations_hist
        local rand_num::Float64

        # First loop uses initial conditions
        local (intensity_value::Float64, temp_value::Float64) = (init_intensity, init_temp)
        # Sample initial starting material
        rand_num = rand(gen_array[threadid()], Float64)
        local material_num::Int32 = (rand_num < prob_1) ? 1 : 2

        # Innter loop
        for (index, time) in enumerate(times)
            local (opacity_term::Float64, spec_heat_term::Float64, dens::Float64, change_prob::Float64) = (material_num == 1) ? (opacity_1, spec_heat_1, dens_1, change_prob_1) : (opacity_2, spec_heat_2, dens_2, change_prob_2)

            # Sample whether material changes
            rand_num = rand(gen_array[threadid()], Float64)

            if (rand_num > change_prob)
                local opacity::Float64 = PhysicsFunctions.sigma_a(opacity_term, temp_value)
                local spec_heat::Float64 = PhysicsFunctions.c_v(spec_heat_term, temp_value)

                local new_intensity_value::Float64 = PhysicsFunctions.balance_intensity(opacity, new_delta_t, intensity_value, temp_value)
                local new_temp_value::Float64 = PhysicsFunctions.balance_temp(opacity, spec_heat, dens, new_delta_t, intensity_value, temp_value)

                (intensity_value, temp_value) = (new_intensity_value, new_temp_value)

                local (intensity_bin_no::Int64, temp_bin_no::Int64)
            else
                material_num = (material_num == 1) ? 2 : 1
            end

            # Tally histogram data
            if (material_num == 1)
                ExponentialHist.push(intensity_1_bin[threadid()], intensity_value)
                ExponentialHist.push(temperature_1_bin[threadid()], temp_value)
                ExponentialHist.push(opacity_1_bin[threadid()], PhysicsFunctions.sigma_a(opacity_1, temp_value))
            else
                ExponentialHist.push(intensity_2_bin[threadid()], intensity_value)
                ExponentialHist.push(temperature_2_bin[threadid()], temp_value)
                ExponentialHist.push(opacity_2_bin[threadid()], PhysicsFunctions.sigma_a(opacity_2, temp_value))
            end
        end

        # Need to reference Core namespace for thread-safe printing
        if (i % num_say_hist == 0)
            lock(printlock) do
                Core.println("History Number ", i)
            end
        end
    end

    local material_1_intensity_bin::Vector{Float64} = zeros(num_bins + 1)
    local material_2_intensity_bin::Vector{Float64} = zeros(num_bins + 1)
    local material_1_temperature_bin::Vector{Float64} = zeros(num_bins + 1)
    local material_2_temperature_bin::Vector{Float64} = zeros(num_bins + 1)
    local material_1_opacity_bin::Vector{Float64} = zeros(num_bins + 1)
    local material_2_opacity_bin::Vector{Float64} = zeros(num_bins + 1)
    for i in 1:nthreads()
        material_1_intensity_bin += ExponentialHist.histogram(intensity_1_bin[i])
        material_2_intensity_bin += ExponentialHist.histogram(intensity_2_bin[i])
        material_1_temperature_bin += ExponentialHist.histogram(temperature_1_bin[i])
        material_2_temperature_bin += ExponentialHist.histogram(temperature_2_bin[i])
        material_1_opacity_bin += ExponentialHist.histogram(opacity_1_bin[i])
        material_2_opacity_bin += ExponentialHist.histogram(opacity_2_bin[i])
    end

    # Normalize histograms
    material_1_intensity_bin ./= sum(material_1_intensity_bin)
    material_2_intensity_bin ./= sum(material_2_intensity_bin)
    material_1_temperature_bin ./= sum(material_1_temperature_bin)
    material_2_temperature_bin ./= sum(material_2_temperature_bin)
    material_1_opacity_bin ./= sum(material_1_opacity_bin)
    material_2_opacity_bin ./= sum(material_2_opacity_bin)

    # The distribution can be computed from the first array point
    local material_1_intensity_array::Vector{Float64} = ExponentialHist.distribution(intensity_1_bin[1])
    local material_2_intensity_array::Vector{Float64} = ExponentialHist.distribution(intensity_2_bin[1])
    local material_1_temperature_array::Vector{Float64} = ExponentialHist.distribution(temperature_1_bin[1])
    local material_2_temperature_array::Vector{Float64} = ExponentialHist.distribution(temperature_2_bin[1])
    local material_1_opacity_array::Vector{Float64} = ExponentialHist.distribution(opacity_1_bin[1])
    local material_2_opacity_array::Vector{Float64} = ExponentialHist.distribution(opacity_2_bin[1])

    # Time (constant array in csv)
    local time_array::Vector{Float64} = fill(steady_state_time, num_bins + 1)

    local tabular::DataFrame = DataFrame(time=time_array, intensity1arr=material_1_intensity_array, freqintensity1=material_1_intensity_bin, intensity2arr=material_2_intensity_array, freqintensity2=material_2_intensity_bin, temperature1arr=material_1_temperature_array, freqtemperature1=material_1_temperature_bin, temperature2arr=material_2_temperature_array, freqtemperature2=material_2_temperature_bin, opacity1arr=material_1_opacity_array, freqopacity1=material_1_opacity_bin, opacity2arr=material_2_opacity_array, freqopacity2=material_2_opacity_bin)

    CSV.write("out/nonlinear/pdf_data/mc_pdf.csv", tabular)

    return nothing
end

main()
