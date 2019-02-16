#!/usr/bin/env julia

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
    local minmax_data_location::String = "out/nonlinear/pdf_data/minmax.csv"

    # User Input for data existence
    println("\n")
    println("Did you remember to run the Monte Carlo computation for a two-pass binning?")
    println("Press Enter to continue (no additional input)")
    println("Press C and Enter to run the Monte Carlo computation using current values in Constants.jl")
    println("Press K and Enter to quit")
    local continue_calc::String = lowercase(chomp(readline()))
    if (continue_calc == "k")
        println("Terminating exectuion")
        return nothing
    elseif (continue_calc == "c")
        println("Running Monte Carlo computation now! Please wait...")
        local cmd::Cmd = Base.julia_cmd()
        run(`$cmd NonLinearMC.jl`)
        println("Finished")
    end

    # Datafile
    local (old_data::DataFrame, minmax_data::DataFrame)

    if (isfile(mc_data_location) && isfile(minmax_data_location))
        old_data = CSV.File(mc_data_location) |> DataFrame
        minmax_data = CSV.File(minmax_data_location) |> DataFrame
    else
        println("Cannot find results of a Monte Carlo run, terminating")
        return nothing
    end

    # Bounds for binning
    local intensity_1_min::Float64 = minmax_data.minint1[]
    local intensity_1_max::Float64 = minmax_data.maxint1[]
    local temperature_1_min::Float64 = minmax_data.mintemp1[]
    local temperature_1_max::Float64 = minmax_data.maxtemp1[]
    local intensity_2_min::Float64 = minmax_data.minint2[]
    local intensity_2_max::Float64 = minmax_data.maxint2[]
    local temperature_2_min::Float64 = minmax_data.mintemp2[]
    local temperature_2_max::Float64 = minmax_data.maxtemp2[]

    # Delta values
    local delta_intensity_1::Float64 = (intensity_1_max - intensity_1_min) / num_bins
    local delta_temperature_1::Float64 = (temperature_1_max - temperature_1_min) / num_bins
    local delta_intensity_2::Float64 = (intensity_2_max - intensity_2_min) / num_bins
    local delta_temperature_2::Float64 = (temperature_2_max - temperature_2_min) / num_bins

    # Arrays for binning
    local intensity_1_bin::Array{Float64, 2} = zeros(nthreads(), num_bins + 1)
    local intensity_2_bin::Array{Float64, 2} = zeros(nthreads(), num_bins + 1)
    local temperature_1_bin::Array{Float64, 2} = zeros(nthreads(), num_bins + 1)
    local temperature_2_bin::Array{Float64, 2} = zeros(nthreads(), num_bins + 1)

    function locate_steady_state(data::Vector{Float64})::Int64
        local point::Float64 = 0.0
        local last_point::Float64 = 0.0
        local sec_last_point::Float64 = 0.0
        local index::Int64 = 0
        local last_index::Int64 = num_t
        local steady_state_search::Bool = true

        # Naive pass through data - three same points in a row = steady state
        while (steady_state_search)
            index += 1
            (point, last_point, sec_last_point) = (data[index], point, last_point)
            if ((point == last_point) && (last_point == sec_last_point))
                last_index = index
                steady_state_search = false
            end
        end

        return last_index
    end

    local intensity_1_ss::Int64 = locate_steady_state(vec([old_data.intensity1...]))
    local intensity_2_ss::Int64 = locate_steady_state(vec([old_data.intensity2...]))
    local temperature_1_ss::Int64 = locate_steady_state(vec([old_data.temperature1...]))
    local temperature_2_ss::Int64 = locate_steady_state(vec([old_data.temperature2...]))

    local steady_state_index::Int64 = max(intensity_1_ss, intensity_2_ss, temperature_1_ss, temperature_2_ss)
    local steady_state_time::Float64 = old_data.time[steady_state_index]  # s
    local new_delta_t::Float64 = (steady_state_time - t_init) / num_t  # s

    println(string("Steady state index found to occur at point ", steady_state_index))

    # BEGIN MONTE CARLO SECOND PASS

    # Iteration condition
    local gen_array::Array{MersenneTwister, 1} = let m::MersenneTwister = MersenneTwister(1234)
        [m; accumulate(Future.randjump, fill(big(10)^20, nthreads() - 1), init=m)]
    end

    # Computational values
    local times::Vector{Float64} = [(x * new_delta_t + t_init) * sol for x in 1:num_t]

    # Probability for material sampling
    local prob_1::Float64 = chord_1 / (chord_1 + chord_2)
    local change_prob_1::Float64 = 1.0 / chord_1 * delta_t
    local change_prob_2::Float64 = 1.0 / chord_2 * delta_t
    if ((change_prob_1 > 1.0) || (change_prob_2 > 1.0))
        println("The value for delta_t is too large for sampling")
        return nothing
    end

    function sigma_a(opacity_term::Float64, temp::Float64)::Float64
        return opacity_term / temp^3
    end

    function c_v(spec_heat_term::Float64, temp::Float64)::Float64
        return spec_heat_term
    end

    function balance_intensity(opacity::Float64, past_intensity::Float64, past_temp::Float64)::Float64
        local term_1::Float64 = delta_t * sol^2 * opacity * arad * past_temp^4
        local term_2::Float64 = (1.0 - delta_t * sol * opacity) * past_intensity

        return term_1 + term_2
    end

    function balance_temp(opacity::Float64, spec_heat::Float64, density::Float64, past_intensity::Float64, past_temp::Float64)::Float64
        local term_1::Float64 = delta_t / (density * spec_heat) * opacity * past_intensity
        local term_2::Float64 = - delta_t / (density * spec_heat) * sol * opacity * arad * past_temp^4
        local term_3::Float64 = past_temp

        return term_1 + term_2 + term_3
    end

    println(string("Proceeding with ", nthreads(), " computational threads..."))

    local printlock::SpinLock = SpinLock()

    function bin_no(value::Float64, delta::Float64, d_minimum::Float64)::Int64
        local number::Int64 = floor((value + delta / 2.0 - d_minimum) / delta) + 1

        # Fit into pre-determined number of bins: excess in either direction is assumed to be maximum or minimum
        number = (number < 1) ? 1 : number
        number = (number > num_bins + 1) ? num_bins + 1 : number

        return number
    end

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
                local opacity::Float64 = sigma_a(opacity_term, temp_value)
                local spec_heat::Float64 = c_v(spec_heat_term, temp_value)

                local new_intensity_value::Float64 = balance_intensity(opacity, intensity_value, temp_value)
                local new_temp_value::Float64 = balance_temp(opacity, spec_heat, dens, intensity_value, temp_value)

                (intensity_value, temp_value) = (new_intensity_value, new_temp_value)

                local (intensity_bin_no::Int64, temp_bin_no::Int64)

                if (material_num == 1)
                    intensity_bin_no = bin_no(intensity_value, delta_intensity_1, intensity_1_min)
                    temp_bin_no = bin_no(temp_value, delta_temperature_1, temperature_1_min)
                    intensity_1_bin[threadid(), intensity_bin_no] += 1.0
                    temperature_1_bin[threadid(), temp_bin_no] += 1.0
                else
                    intensity_bin_no = bin_no(intensity_value, delta_intensity_2, intensity_2_min)
                    temp_bin_no = bin_no(temp_value, delta_temperature_2, temperature_2_min)
                    intensity_2_bin[threadid(), intensity_bin_no] += 1.0
                    temperature_2_bin[threadid(), temp_bin_no] += 1.0
                end
            else
                material_num = (material_num == 1) ? 2 : 1
            end
        end

        # Need to reference Core namespace for thread-safe printing
        if (i % num_say == 0)
            lock(printlock) do
                Core.println(string("Iteration Number ", i))
            end
        end
    end

    local material_1_intensity_bin::Vector{Float64} = vec(sum(intensity_1_bin, dims=1))
    local material_2_intensity_bin::Vector{Float64} = vec(sum(intensity_2_bin, dims=1))
    local material_1_temperature_bin::Vector{Float64} = vec(sum(temperature_1_bin, dims=1))
    local material_2_temperature_bin::Vector{Float64} = vec(sum(temperature_2_bin, dims=1))

    material_1_intensity_bin /= sum(material_1_intensity_bin)
    material_2_intensity_bin /= sum(material_2_intensity_bin)
    material_1_temperature_bin /= sum(material_1_temperature_bin)
    material_2_temperature_bin /= sum(material_2_temperature_bin)

    local material_1_intensity_array::Vector{Float64} = [i * delta_intensity_1 + intensity_1_min - delta_intensity_1 / 2.0 for i in 1:(num_bins + 1)]
    local material_2_intensity_array::Vector{Float64} = [i * delta_intensity_2 + intensity_2_min - delta_intensity_2 / 2.0 for i in 1:(num_bins + 1)]
    local material_1_temperature_array::Vector{Float64} = [i * delta_temperature_1 + temperature_1_min - delta_temperature_1 / 2.0 for i in 1:(num_bins + 1)]
    local material_2_temperature_array::Vector{Float64} = [i * delta_temperature_2 + temperature_2_min - delta_temperature_2 / 2.0 for i in 1:(num_bins + 1)]

    tabular::DataFrame = DataFrame(intensity1arr=material_1_intensity_array, freqintensity1=material_1_intensity_bin, intensity2arr=material_2_intensity_array, freqintensity2=material_2_intensity_bin, temperature1arr=material_1_temperature_array, freqtemperature1=material_1_temperature_bin, temperature2arr=material_2_temperature_array, freqtemperature2=material_2_temperature_bin)

    CSV.write("out/nonlinear/pdf_data/mc_pdf.csv", tabular)

    return nothing
end

main()
