#!/usr/bin/env julia

include("GeometryGen.jl")
include("MeshMap.jl")
include("Constants.jl")
#using Plots
using LinearAlgebra
using DataFrames
using CSV

function main()::Nothing
    # Computational values
    local material_intensity::Vector{Float64} = zeros(num_t)
    local material_temp::Vector{Float64} = zeros(num_t)
    local times::Vector{Float64} = [(x * delta_t + t_init) * sol for x in 1:num_t]

    function sigma_a(opacity_term::Float64, temp::Float64)::Float64
        return opacity_term / temp^3
    end

    function c_v(spec_heat_term::Float64, temp::Float64)::Float64
        return spec_heat_term
    end

    function balance_intensity(delta_t::Float64, opacity::Float64, past_intensity::Float64, past_temp::Float64)::Float64
        local term_1::Float64 = delta_t * sol^2 * opacity * arad * past_temp^4
        local term_2::Float64 = (1.0 - delta_t * sol * opacity) * past_intensity

        return term_1 + term_2
    end

    function balance_temp(delta_t::Float64, opacity::Float64, spec_heat::Float64, density::Float64, past_intensity::Float64, past_temp::Float64)::Float64
        local term_1::Float64 = delta_t / (density * spec_heat) * opacity * past_intensity
        local term_2::Float64 = - delta_t / (density * spec_heat) * sol * opacity * arad * past_temp^4
        local term_3::Float64 = past_temp

        return term_1 + term_2 + term_3
    end

    # Loop through time discretizations
    local intensity_value::Float64 = init_intensity
    local temp_value::Float64 = init_temp
    local opacity_term::Float64 = volfrac_1 * opacity_1 + volfrac_2 * opacity_2
    local spec_heat_term::Float64 = volfrac_1 * spec_heat_1 + volfrac_2 * spec_heat_2
    local dens::Float64 = volfrac_1 * dens_1 + volfrac_2 * dens_2
    for index = 1:num_t
        local opacity::Float64 = sigma_a(opacity_term, temp_value)
        local spec_heat::Float64 = c_v(spec_heat_term, temp_value)

        local new_intensity_value::Float64 = balance_intensity(delta_t, opacity, intensity_value, temp_value)
        local new_temp_value::Float64 = balance_temp(delta_t, opacity, spec_heat, dens, intensity_value, temp_value)

        (intensity_value, temp_value) = (new_intensity_value, new_temp_value)

        material_intensity[index] = intensity_value  # erg/cm^2-s
        material_temp[index] = temp_value  # eV
    end

    function plot_limits_log10(series::Vector{Float64})::Tuple{Float64, Float64}
        local new_series::Vector{Float64} = log10.(series)
        local minimum_value::Float64 = floor(minimum(new_series))
        local maximum_value::Float64 = ceil(maximum(new_series))

        return (10.0^minimum_value, 10.0^maximum_value)
    end

    tabular::DataFrame = DataFrame(time=times, intensity=material_intensity, temperature=material_temp)

    CSV.write("out/nonlinear/data/homog_nonlinear_exp.csv", tabular)

#    plot(times, material_intensity, label="Intensity", lc=:magenta, ls=:solid, lw=3)
#    plot!(yscale=:log10, xscale=:log10,
#        title="Explicit Homogeneous Intensity Plot", xlabel="Time - ct (cm)", ylabel="Intensity (erg/cm^2-s)",
#        legend=:best, ylims=plot_limits_log10(material_intensity))
#    png("out/nonlinear/exp_homog_intensity")

#    plot(times, material_temp, label="Temperature", lc=:green, ls=:solid, lw=3)
#    plot!(yscale=:log10, xscale=:log10,
#        title="Explicit Homogeneous Temperature Plot", xlabel="Time - ct (cm)", ylabel="Temperature (eV)",
#        legend=:best, ylims=plot_limits_log10(material_temp))
#    png("out/nonlinear/exp_homog_temperature")

    return nothing
end

main()
