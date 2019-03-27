#!/usr/bin/env julia

include("ExponentialHist.jl")
include("GeometryGen.jl")
include("MeshMap.jl")
using Random
using Future
using LinearAlgebra
using DataFrames
using CSV
using Base.Threads
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

    # Arrays for binning
    local intensity_1_bin::Vector{ExponentialHist.ExpHist} = [ExponentialHist.ExpHist(num_bins, intensity_1_min, intensity_1_max) for i in 1:nthreads()]
    local intensity_2_bin::Vector{ExponentialHist.ExpHist} = [ExponentialHist.ExpHist(num_bins, intensity_2_min, intensity_2_max) for i in 1:nthreads()]
    local temperature_1_bin::Vector{ExponentialHist.ExpHist} = [ExponentialHist.ExpHist(num_bins, temperature_1_min, temperature_1_max) for i in 1:nthreads()]
    local temperature_2_bin::Vector{ExponentialHist.ExpHist} = [ExponentialHist.ExpHist(num_bins, temperature_2_min, temperature_2_max) for i in 1:nthreads()]

    function locate_steady_state(data::Vector{Float64})::Int64
        local point::Float64 = 0.0
        local last_point::Float64 = 0.0
        local sec_last_point::Float64 = 0.0
        local index::Int64 = 0
        local last_index::Int64 = num_t
        local steady_state_search::Bool = true

        # Naive pass through data - three same points in a row = steady_state
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
    local steady_state_time::Float64 = old_data.time[steady_state_index] / sol  # s
    local new_delta_t::Float64 = (steady_state_time - t_init) / num_t_hist  # s

    println("Steady state index found to occur at point ", steady_state_index)

    # BEGIN REALIZATIONS SECOND PASS

    # Iteration condition
    local gen_array::Array{MersenneTwister, 1} = let m::MersenneTwister = MersenneTwister(1234)
        [m; accumulate(Future.randjump, fill(big(10)^20, nthreads() - 1), init=m)]
    end

    # Computational values
    local times::Vector{Float64} = [(x * new_delta_t + t_init) * sol for x in 1:num_t_hist]

    function sigma_a(opacity_term::Float64, temp::Float64)::Float64
        return opacity_term / temp^3
    end

    function c_v(spec_heat_term::Float64, temp::Float64)::Float64
        return spec_heat_term
    end

    function balance_a(intensity::Float64, temp::Float64, delta_t::Float64, opacity::Float64, past_intensity::Float64)::Float64
        local term_1::Float64 = intensity
        local term_2::Float64 = - delta_t * sol^2 * opacity * arad * temp^4
        local term_3::Float64 = delta_t * sol * opacity * intensity
        local term_4::Float64 = - past_intensity

        return term_1 + term_2 + term_3 + term_4
    end

    function balance_b(intensity::Float64, temp::Float64, delta_t::Float64, opacity::Float64, spec_heat::Float64, density::Float64, past_temp::Float64)::Float64
        local term_1::Float64 = temp
        local term_2::Float64 = - delta_t / (density * spec_heat) * opacity * intensity
        local term_3::Float64 = delta_t / (density * spec_heat) * sol * opacity * arad * temp^4
        local term_4::Float64 = - past_temp

        return term_1 + term_2 + term_3 + term_4
    end

    function make_jacobian(intensity::Float64, temp::Float64, delta_t::Float64, opacity_term::Float64, density::Float64, spec_heat_term::Float64)::Array{Float64, 2}
        local term_1::Float64 = 1.0 + delta_t * sol * opacity_term / temp^3
        local term_2::Float64 = - delta_t * sol^2 * opacity_term * arad - 3.0 * delta_t * sol * opacity_term / temp^4 * intensity
        local term_3::Float64 = - delta_t / (density * spec_heat_term) * opacity_term / temp^3
        local term_4::Float64 = 1.0 + 3.0 * delta_t * opacity_term / (density * spec_heat_term * temp^4) * intensity + delta_t / (density * spec_heat_term) * sol * opacity_term * arad

        return [
            term_1 term_2;
            term_3 term_4
        ]
    end

    function relative_change(delta::Vector{Float64}, original_terms::Vector{Float64})::Float64
        local current_norm::Float64 = sqrt(sum([x^2 for x in delta]))
        local original_norm::Float64 = sqrt(sum([x^2 for x in original_terms]))

        return current_norm / original_norm
    end

    println("Proceeding with ", nthreads(), " computational threads...")

    local printlock::SpinLock = SpinLock()

    # Outer loop
    @threads for i = 1:max_iterations_hist
        local (t_delta::Vector{Float64}, t_arr::Vector{Float64}, materials::Vector{Int32}, num_cells::Int64) = GeometryGen.get_geometry(chord_1, chord_2, steady_state_time, num_divs_hist, rng=gen_array[threadid()])

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
                local opacity::Float64 = sigma_a(opacity_term, old_terms[2])
                local spec_heat::Float64 = c_v(spec_heat_term, old_terms[2])

                local jacobian::Array{Float64, 2} = make_jacobian(old_terms[1], old_terms[2], delta_t_unstruct, opacity_term, dens, spec_heat_term)
                local func_vector::Vector{Float64} = [
                    balance_a(old_terms[1], old_terms[2], delta_t_unstruct, opacity, intensity_value),
                    balance_b(old_terms[1], old_terms[2], delta_t_unstruct, opacity, spec_heat, dens, temp_value)
                ]

                local delta::Vector{Float64} = jacobian \ - func_vector

                new_terms = delta + old_terms
                old_terms = new_terms

                error = relative_change(delta, original_terms)
            end
            intensity_value = new_terms[1]  # erg/cm^2-s
            temp_value = new_terms[2]  # eV
            #intensity[index] = intensity_value  # erg/cm^2-s
            #temp[index] = temp_value  # eV

            if (material == 1)
                ExponentialHist.push(intensity_1_bin[threadid()], intensity_value)
                ExponentialHist.push(temperature_1_bin[threadid()], temp_value)
            else
                ExponentialHist.push(intensity_2_bin[threadid()], intensity_value)
                ExponentialHist.push(temperature_2_bin[threadid()], temp_value)
            end
        end

        #local material_intensity_array::Array{Float64, 2} = MeshMap.material_calc(intensity, t_delta, num_cells, materials, new_delta_t, num_t_hist, convert(Int32, 2))
        #local material_temp_array::Array{Float64, 2} = MeshMap.material_calc(temp, t_delta, num_cells, materials, new_delta_t, num_t_hist, convert(Int32, 2))

        #for k in 1:num_t_hist
        #    if (material_intensity_array[k, 1] != 0.0)
        #        ExponentialHist.push(intensity_1_bin[threadid()], material_intensity_array[k, 1])
        #        ExponentialHist.push(temperature_1_bin[threadid()], material_temp_array[k, 1])
        #    else
        #        ExponentialHist.push(intensity_2_bin[threadid()], material_intensity_array[k, 2])
        #        ExponentialHist.push(temperature_2_bin[threadid()], material_temp_array[k, 2])
        #    end
        #end

        # Need to reference Core namespace for thread-safe printing
        if (i % num_say == 0)
            lock(printlock) do
                Core.println("Iteration Number ", i)
            end
        end
    end

    local material_1_intensity_bin::Vector{Float64} = zeros(num_bins + 1)
    local material_2_intensity_bin::Vector{Float64} = zeros(num_bins + 1)
    local material_1_temperature_bin::Vector{Float64} = zeros(num_bins + 1)
    local material_2_temperature_bin::Vector{Float64} = zeros(num_bins + 1)
    for i in 1:nthreads()
        material_1_intensity_bin += ExponentialHist.histogram(intensity_1_bin[i])
        material_2_intensity_bin += ExponentialHist.histogram(intensity_2_bin[i])
        material_1_temperature_bin += ExponentialHist.histogram(temperature_1_bin[i])
        material_2_temperature_bin += ExponentialHist.histogram(temperature_2_bin[i])
    end

    # Normalize histograms
    material_1_intensity_bin ./= sum(material_1_intensity_bin)
    material_2_intensity_bin ./= sum(material_2_intensity_bin)
    material_1_temperature_bin ./= sum(material_1_temperature_bin)
    material_2_temperature_bin ./= sum(material_2_temperature_bin)

    # The distribution can be computed from the first array point
    local material_1_intensity_array::Vector{Float64} = ExponentialHist.distribution(intensity_1_bin[1])
    local material_2_intensity_array::Vector{Float64} = ExponentialHist.distribution(intensity_2_bin[1])
    local material_1_temperature_array::Vector{Float64} = ExponentialHist.distribution(temperature_1_bin[1])
    local material_2_temperature_array::Vector{Float64} = ExponentialHist.distribution(temperature_2_bin[1])

    # Compute opacities as a function of temperature
    local xs_1_array::Vector{Float64} = @. sigma_a(opacity_1, material_1_temperature_array)
    local xs_2_array::Vector{Float64} = @. sigma_a(opacity_2, material_2_temperature_array)

    local tabular::DataFrame = DataFrame(intensity1arr=material_1_intensity_array, freqintensity1=material_1_intensity_bin, intensity2arr=material_2_intensity_array, freqintensity2=material_2_intensity_bin, temperature1arr=material_1_temperature_array, freqtemperature1=material_1_temperature_bin, temperature2arr=material_2_temperature_array, freqtemperature2=material_2_temperature_bin, opacity1arr=xs_1_array, opacity2arr=xs_2_array)

    CSV.write("out/nonlinear/pdf_data/realizations_pdf.csv", tabular)

    return nothing
end

main()
