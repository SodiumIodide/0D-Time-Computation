#!/usr/bin/env julia

include("Histogram.jl")
include("PhysicsFunctions.jl")
include("GeometryGen.jl")
include("PDFFunctions.jl")
using Random
using LinearAlgebra
using DataFrames
using CSV
include("Constants.jl")

function main()::Nothing
    # Read previous data to explicitly create binned elements
    # The Monte Carlo data is read, because it is the faster computation method
    # (in comparison to realization generation)
    local mc_data_location::String = "out/nonlinear/data/nonlinearmc.csv"
    local minmax_data_location::String = "out/nonlinear/pdf_data/minmax.csv"

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
    local intensity_1_bin::Histogram.Hist = Histogram.Hist(num_bins, intensity_1_min, intensity_1_max)
    local intensity_2_bin::Histogram.Hist = Histogram.Hist(num_bins, intensity_2_min, intensity_2_max)
    local temperature_1_bin::Histogram.Hist = Histogram.Hist(num_bins, temperature_1_min, temperature_1_max)
    local temperature_2_bin::Histogram.Hist = Histogram.Hist(num_bins, temperature_2_min, temperature_2_max)
    local opacity_1_bin::Histogram.Hist = Histogram.Hist(num_bins, opacity_1_min, opacity_1_max)
    local opacity_2_bin::Histogram.Hist = Histogram.Hist(num_bins, opacity_2_min, opacity_2_max)

    # BEGIN REALIZATIONS SECOND PASS

    # Iteration condition
    local cont_calc_outer::Bool = true
    local generator::MersenneTwister = MersenneTwister(1234)

    # Computational values
    local iteration_number::Int64 = 0
    local new_delta_t::Float64 = (steady_state_time / sol) / num_t_hist
    local times::Vector{Float64} = [(x * new_delta_t + t_init) * sol for x in 1:num_t_hist]

    # Outer loop
    while (cont_calc_outer)
        local (t_delta::Vector{Float64}, t_arr::Vector{Float64}, materials::Vector{Int32}, num_cells::Int64) = GeometryGen.get_geometry(chord_1, chord_2, steady_state_time, num_divs_hist, rng=generator)

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

                local delta::Vector{Float64} = jacobian \ - func_vector

                new_terms = delta + old_terms
                old_terms = new_terms

                error = PhysicsFunctions.relative_change(delta, original_terms)
            end
            intensity_value = new_terms[1]  # erg/cm^2-s
            temp_value = new_terms[2]  # eV

            if (material == 1)
                Histogram.push(intensity_1_bin, intensity_value)
                Histogram.push(temperature_1_bin, temp_value)
                Histogram.push(opacity_1_bin, PhysicsFunctions.sigma_a(opacity_1, temp_value))
            else
                Histogram.push(intensity_2_bin, intensity_value)
                Histogram.push(temperature_2_bin, temp_value)
                Histogram.push(opacity_2_bin, PhysicsFunctions.sigma_a(opacity_2, temp_value))
            end
        end

        iteration_number += 1

        if (iteration_number % num_say_hist == 0)
            println("Realization Number ", iteration_number)
        end

        if (iteration_number > max_iterations_hist)
            cont_calc_outer = false
        end
    end

    local material_1_intensity_bin::Vector{Float64} = Histogram.histogram(intensity_1_bin)
    local material_2_intensity_bin::Vector{Float64} = Histogram.histogram(intensity_2_bin)
    local material_1_temperature_bin::Vector{Float64} = Histogram.histogram(temperature_1_bin)
    local material_2_temperature_bin::Vector{Float64} = Histogram.histogram(temperature_2_bin)
    local material_1_opacity_bin::Vector{Float64} = Histogram.histogram(opacity_1_bin)
    local material_2_opacity_bin::Vector{Float64} = Histogram.histogram(opacity_2_bin)

    # Normalize histograms
    material_1_intensity_bin ./= sum(material_1_intensity_bin)
    material_2_intensity_bin ./= sum(material_2_intensity_bin)
    material_1_temperature_bin ./= sum(material_1_temperature_bin)
    material_2_temperature_bin ./= sum(material_2_temperature_bin)
    material_1_opacity_bin ./= sum(material_1_opacity_bin)
    material_2_opacity_bin ./= sum(material_2_opacity_bin)

    local material_1_intensity_array::Vector{Float64} = Histogram.distribution(intensity_1_bin)
    local material_2_intensity_array::Vector{Float64} = Histogram.distribution(intensity_2_bin)
    local material_1_temperature_array::Vector{Float64} = Histogram.distribution(temperature_1_bin)
    local material_2_temperature_array::Vector{Float64} = Histogram.distribution(temperature_2_bin)
    local material_1_opacity_array::Vector{Float64} = Histogram.distribution(opacity_1_bin)
    local material_2_opacity_array::Vector{Float64} = Histogram.distribution(opacity_2_bin)

    # Time (constant array in csv)
    local time_array::Vector{Float64} = fill(steady_state_time, num_bins + 1)

    local tabular::DataFrame = DataFrame(time=time_array, intensity1arr=material_1_intensity_array, freqintensity1=material_1_intensity_bin, intensity2arr=material_2_intensity_array, freqintensity2=material_2_intensity_bin, temperature1arr=material_1_temperature_array, freqtemperature1=material_1_temperature_bin, temperature2arr=material_2_temperature_array, freqtemperature2=material_2_temperature_bin, opacity1arr=material_1_opacity_array, freqopacity1=material_1_opacity_bin, opacity2arr=material_2_opacity_array, freqopacity2=material_2_opacity_bin)

    CSV.write("out/nonlinear/pdf_data/realizations_pdf.csv", tabular)

    return nothing
end

main()
