#!/usr/bin/env julia

include("Histogram.jl")
include("PhysicsFunctions.jl")
include("PDFFunctions.jl")
using Random
using DataFrames
using CSV
using ProgressMeter
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

    local steady_state_index::Int64 = @fastmath maximum([intensity_1_ss, intensity_2_ss, temperature_1_ss, temperature_2_ss])

    #local steady_state_index::Int64 = PDFFunctions.steady_state_fix(vec([old_data.time...]))

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
    local opacity_1_min::Float64 = old_data.minopacity1[chosenindex]
    local opacity_1_max::Float64 = old_data.maxopacity1[chosenindex]
    local opacity_2_min::Float64 = old_data.minopacity2[chosenindex]
    local opacity_2_max::Float64 = old_data.maxopacity2[chosenindex]

    #=
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
    =#

    # Arrays for binning
    local intensity_1_bin::Histogram.Hist = Histogram.Hist(num_bins, intensity_1_min, intensity_1_max)
    local intensity_2_bin::Histogram.Hist = Histogram.Hist(num_bins, intensity_2_min, intensity_2_max)
    local temperature_1_bin::Histogram.Hist = Histogram.Hist(num_bins, temperature_1_min, temperature_1_max)
    local temperature_2_bin::Histogram.Hist = Histogram.Hist(num_bins, temperature_2_min, temperature_2_max)
    local opacity_1_bin::Histogram.Hist = Histogram.Hist(num_bins, opacity_1_min, opacity_1_max)
    local opacity_2_bin::Histogram.Hist = Histogram.Hist(num_bins, opacity_2_min, opacity_2_max)

    # BEGIN MONTE CARLO SECOND PASS

    # Iteration condition
    local generator::MersenneTwister = MersenneTwister(1234)

    # Computational values
    local new_delta_t::Float64 = @fastmath (chosentime / sol) / num_t_hist
    local times::Vector{Float64} = @fastmath [(x * new_delta_t + t_init) * sol for x in 1:num_t_hist]

    # Probability for material sampling
    local prob_1::Float64 = @fastmath chord_1 / (chord_1 + chord_2)
    local change_prob_1::Float64 = @fastmath 1.0 / chord_1 * new_delta_t
    local change_prob_2::Float64 = @fastmath 1.0 / chord_2 * new_delta_t
    if @fastmath((change_prob_1 > 1.0) || (change_prob_2 > 1.0))
        println("The value for delta_t is too large for sampling")
        return nothing
    end

    # Outer loop
    @showprogress 1 for iteration_number=1:max_iterations_hist
        local rand_num::Float64

        # First loop uses initial conditions
        local (intensity_value::Float64, temp_value::Float64) = (init_intensity, init_temp)
        # Sample initial starting material
        rand_num = @fastmath rand(generator, Float64)
        local material_num::Int32 = @fastmath (rand_num < prob_1) ? 1 : 2

        # Innter loop
        for (index, time) in enumerate(times)
            local (opacity_term::Float64, spec_heat_term::Float64, dens_term::Float64, change_prob::Float64) = @fastmath (material_num == 1) ? (opacity_1, spec_heat_1, dens_1, change_prob_1) : (opacity_2, spec_heat_2, dens_2, change_prob_2)

            # Sample whether material changes
            rand_num = @fastmath rand(generator, Float64)

            if @fastmath(rand_num > change_prob)
                local opacity::Float64 = PhysicsFunctions.sigma_a(opacity_term, temp_value)
                local spec_heat::Float64 = PhysicsFunctions.c_v(spec_heat_term, temp_value)
                local dens::Float64 = PhysicsFunctions.rho(dens_term, temp_value)

                local new_intensity_value::Float64 = PhysicsFunctions.balance_intensity(opacity, new_delta_t, intensity_value, temp_value)
                local new_temp_value::Float64 = PhysicsFunctions.balance_temp(opacity, spec_heat, dens, new_delta_t, intensity_value, temp_value)

                (intensity_value, temp_value) = (new_intensity_value, new_temp_value)
            else
                material_num = @fastmath (material_num == 1) ? 2 : 1
            end
        end

        # Tally histogram data
        if @fastmath(material_num == 1)
            Histogram.push(intensity_1_bin, intensity_value)
            Histogram.push(temperature_1_bin, temp_value)
            Histogram.push(opacity_1_bin, PhysicsFunctions.sigma_a(opacity_1, temp_value))
        else
            Histogram.push(intensity_2_bin, intensity_value)
            Histogram.push(temperature_2_bin, temp_value)
            Histogram.push(opacity_2_bin, PhysicsFunctions.sigma_a(opacity_2, temp_value))
        end
    end

    local material_1_intensity_bin::Vector{Float64} = Histogram.histogram(intensity_1_bin)
    local material_2_intensity_bin::Vector{Float64} = Histogram.histogram(intensity_2_bin)
    local material_1_temperature_bin::Vector{Float64} = Histogram.histogram(temperature_1_bin)
    local material_2_temperature_bin::Vector{Float64} = Histogram.histogram(temperature_2_bin)
    local material_1_opacity_bin::Vector{Float64} = Histogram.histogram(opacity_1_bin)
    local material_2_opacity_bin::Vector{Float64} = Histogram.histogram(opacity_2_bin)

    # Normalize histograms
    local tot_m1_intensity::Float64 = @fastmath sum(material_1_intensity_bin)
    local tot_m2_intensity::Float64 = @fastmath sum(material_2_intensity_bin)
    local tot_m1_temp::Float64 = @fastmath sum(material_1_temperature_bin)
    local tot_m2_temp::Float64 = @fastmath sum(material_2_temperature_bin)
    local tot_m1_opacity::Float64 = @fastmath sum(material_1_opacity_bin)
    local tot_m2_opacity::Float64 = @fastmath sum(material_2_opacity_bin)
    @fastmath material_1_intensity_bin ./= tot_m1_intensity
    @fastmath material_2_intensity_bin ./= tot_m2_intensity
    @fastmath material_1_temperature_bin ./= tot_m1_temp
    @fastmath material_2_temperature_bin ./= tot_m2_temp
    @fastmath material_1_opacity_bin ./= tot_m1_opacity
    @fastmath material_2_opacity_bin ./= tot_m2_opacity

    local material_1_intensity_array::Vector{Float64} = Histogram.distribution(intensity_1_bin)
    local material_2_intensity_array::Vector{Float64} = Histogram.distribution(intensity_2_bin)
    local material_1_temperature_array::Vector{Float64} = Histogram.distribution(temperature_1_bin)
    local material_2_temperature_array::Vector{Float64} = Histogram.distribution(temperature_2_bin)
    local material_1_opacity_array::Vector{Float64} = Histogram.distribution(opacity_1_bin)
    local material_2_opacity_array::Vector{Float64} = Histogram.distribution(opacity_2_bin)

    local tabular::DataFrame = DataFrame(time=time_array, intensity1arr=material_1_intensity_array, freqintensity1=material_1_intensity_bin, intensity2arr=material_2_intensity_array, freqintensity2=material_2_intensity_bin, temperature1arr=material_1_temperature_array, freqtemperature1=material_1_temperature_bin, temperature2arr=material_2_temperature_array, freqtemperature2=material_2_temperature_bin, opacity1arr=material_1_opacity_array, freqopacity1=material_1_opacity_bin, opacity2arr=material_2_opacity_array, freqopacity2=material_2_opacity_bin)

    CSV.write("out/nonlinear/pdf_data/mc_pdf.csv", tabular)

    return nothing
end

main()
