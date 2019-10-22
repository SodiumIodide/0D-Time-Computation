#!/usr/bin/env julia

include("PhysicsFunctions.jl")
using DataFrames
using CSV
include("Constants.jl")

# Explicit Solver
function main()::Nothing
    set_zero_subnormals(true)

    # Computational values
    local times::Vector{Float64} = @fastmath [(x * delta_t + t_init) * sol for x in 1:num_t]

    local dens_term::Float64 = volfrac_1 * dens_1 + volfrac_2 * dens_2  # g/cm^3
    local opacity_term::Float64 = volfrac_1 * opacity_1 + volfrac_2 * opacity_2  # cm^-1
    local spec_heat_term::Float64 = volfrac_1 * spec_heat_1 + volfrac_2 * spec_heat_2  # erg/g-eV

    # Parallel arrays
    local intensity_arr::Vector{Float64} = zeros(num_t)
    local temperature_arr::Vector{Float64} = zeros(num_t)
    local opacity_arr::Vector{Float64} = zeros(num_t)

    # First value uses initial conditions
    local (intensity_value::Float64, temp_value::Float64) = (init_intensity, init_temp)

    # Loop
    for index = 1:num_t
        local opacity::Float64 = PhysicsFunctions.sigma_a(opacity_term, temp_value)
        local spec_heat::Float64 = PhysicsFunctions.c_v(spec_heat_term, temp_value)
        local dens::Float64 = PhysicsFunctions.rho(dens_term, temp_value)

        local new_intensity_value::Float64 = PhysicsFunctions.balance_intensity(opacity, delta_t, intensity_value, temp_value)
        local new_temp_value::Float64 = PhysicsFunctions.balance_temp(opacity, spec_heat, dens, delta_t, intensity_value, temp_value)

        (intensity_value, temp_value) = (new_intensity_value, new_temp_value)

        intensity_arr[index] = new_intensity_value  # erg/cm^2-s
        temperature_arr[index] = new_temp_value  # eV
        opacity_arr[index] = opacity  # cm^-1
    end

    local tabular::DataFrame = DataFrame(time=times, intensityarr=intensity_arr, temperaturearr=temperature_arr, opacityarr=opacity_arr)

    CSV.write("out/nonlinear/data/atomicmix.csv", tabular)

    return nothing
end

main()
