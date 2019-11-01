#!/usr/bin/env julia

include("PhysicsFunctions.jl")
using DataFrames
using CSV
include("Constants.jl")

@inline function explicit_intensity1(intensity1_prev::Float64, intensity2_prev::Float64, temp1_prev::Float64, opacity::Float64)::Float64
    local emission::Float64 = @fastmath delta_t * sol^2 * opacity * arad * temp1_prev^4
    local absorption::Float64 = @fastmath delta_t * sol * destruction_factor * opacity * intensity1_prev
    local lp_term::Float64 = @fastmath delta_t * (volfrac_2 * intensity2_prev / chord_2 - volfrac_1 * intensity1_prev / chord_1)
    local return_value = @fastmath intensity1_prev + emission - absorption + lp_term

    if @fastmath (return_value < 0.0)
        return_value = 0.0
    end

    return return_value
end

@inline function explicit_intensity2(intensity2_prev::Float64, intensity1_prev::Float64, temp2_prev::Float64, opacity::Float64)::Float64
    local emission::Float64 = @fastmath delta_t * sol^2 * opacity * arad * temp2_prev^4
    local absorption::Float64 = @fastmath delta_t * sol * destruction_factor * opacity * intensity2_prev
    local lp_term::Float64 = @fastmath delta_t * (volfrac_1 * intensity1_prev / chord_1 - volfrac_2 * intensity2_prev / chord_2)
    local return_value = @fastmath intensity2_prev + emission - absorption + lp_term

    if @fastmath (return_value < 0.0)
        return_value = 0.0
    end

    return return_value
end

@inline function explicit_temp1(temp1_prev::Float64, temp2_prev::Float64, intensity1_prev::Float64, opacity::Float64, dens::Float64, spec_heat::Float64)::Float64
    local alpha::Float64 = @fastmath delta_t / (dens * spec_heat)
    local absorption::Float64 = @fastmath alpha * opacity * intensity1_prev
    local emission::Float64 = @fastmath alpha * sol * opacity * arad * temp1_prev^4
    local lp_term::Float64 = @fastmath delta_t * (volfrac_2 * temp2_prev / chord_2 - volfrac_1 * temp1_prev / chord_1)

    return @fastmath temp1_prev + absorption - emission + lp_term
end

@inline function explicit_temp2(temp2_prev::Float64, temp1_prev::Float64, intensity2_prev::Float64, opacity::Float64, dens::Float64, spec_heat::Float64)::Float64
    local alpha::Float64 = @fastmath delta_t / (dens * spec_heat)
    local absorption::Float64 = @fastmath alpha * opacity * intensity2_prev
    local emission::Float64 = @fastmath alpha * sol * opacity * arad * temp2_prev^4
    local lp_term::Float64 = @fastmath delta_t * (volfrac_1 * temp1_prev / chord_1 - volfrac_2 * temp2_prev / chord_2)

    return @fastmath temp2_prev + absorption - emission + lp_term
end

# Explicit solver
function main()::Nothing
    set_zero_subnormals(true)

    # Computational values
    local times::Vector{Float64} = @fastmath [(x * delta_t + t_init) * sol for x in 1:num_t]

    # Parallel arrays
    local intensity_arr::Array{Float64, 2} = zeros(num_t, 2)
    local temperature_arr::Array{Float64, 2} = zeros(num_t, 2)
    local opacity_arr::Array{Float64, 2} = zeros(num_t, 2)

    # First value uses initial conditions
    local (intensity_1::Float64, intensity_2::Float64, temp_1::Float64, temp_2::Float64) = (init_intensity, init_intensity, init_temp, init_temp)

    # Loop
    for index = 1:num_t
        local sigma_1::Float64 = PhysicsFunctions.sigma_a(opacity_1, temp_1)
        local sigma_2::Float64 = PhysicsFunctions.sigma_a(opacity_2, temp_2)
        local c_v_1 = PhysicsFunctions.c_v(spec_heat_1, temp_1)
        local c_v_2 = PhysicsFunctions.c_v(spec_heat_2, temp_2)
        local rho_1 = PhysicsFunctions.rho(dens_1, temp_1)
        local rho_2 = PhysicsFunctions.rho(dens_2, temp_2)

        local new_intensity_1::Float64 = explicit_intensity1(intensity_1, intensity_2, temp_1, sigma_1)
        local new_intensity_2::Float64 = explicit_intensity2(intensity_2, intensity_1, temp_2, sigma_2)
        local new_temp_1::Float64 = explicit_temp1(temp_1, temp_2, intensity_1, sigma_1, rho_1, c_v_1)
        local new_temp_2::Float64 = explicit_temp2(temp_2, temp_1, intensity_2, sigma_2, dens_2, c_v_2)

        (intensity_1, intensity_2) = (new_intensity_1, new_intensity_2)
        (temp_1, temp_2) = (new_temp_1, new_temp_2)

        intensity_arr[index, 1] = new_intensity_1  # erg/cm^2-s
        intensity_arr[index, 2] = new_intensity_2  # erg/cm^2-s
        temperature_arr[index, 1] = new_temp_1  # eV
        temperature_arr[index, 2] = new_temp_2  # eV
        opacity_arr[index, 1] = sigma_1  # cm^-1
        opacity_arr[index, 2] = sigma_2  # cm^-1
    end

    local avg_intensity::Vector{Float64} = @fastmath @. intensity_arr[:, 1] * volfrac_1 + intensity_arr[:, 2] * volfrac_2
    local avg_temperature::Vector{Float64} = @fastmath @. temperature_arr[:, 1] * volfrac_1 + temperature_arr[:, 2] * volfrac_2
    local avg_opacity::Vector{Float64} = @fastmath @. opacity_arr[:, 1] * volfrac_1 + opacity_arr[:, 2] * volfrac_2

    local tabular::DataFrame = DataFrame(time=times, intensity1=intensity_arr[:, 1], intensity2=intensity_arr[:, 2], intensityarr=avg_intensity, temperature1=temperature_arr[:, 1], temperature2=temperature_arr[:, 2], temperaturearr=avg_temperature, opacity1=opacity_arr[:, 1], opacity2=opacity_arr[:, 2], opacityarr=avg_opacity)

    CSV.write("out/nonlinear/data/levermorepomraning.csv", tabular)

    return nothing
end

main()
