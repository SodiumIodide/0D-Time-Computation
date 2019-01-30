#!/usr/bin/env julia

# Problem constants
global const sol = 29979245800.0  # cm/s
global const arad = 137.2017050419133  # erg/cm^2-eV^4
global const t_max = 1e-11  # s
global const t_init = 0.0  # s
global const num_t = convert(Int64, 1e4)
global const delta_t = (t_max - t_init) / num_t  # s
global const chord_1 = 1.0 / sol  # s
global const chord_2 = 5.0 / sol  # s
global const dens_1 = 1.0  # g/cm^3
global const dens_2 = 1.0  # g/cm^3
global const opacity_1 = 1.0  # cm^-1
global const opacity_2 = 5.0  # cm^-1
global const spec_heat_1 = 1.0  # erg/g-eV
global const spec_heat_2 = 1.0  # erg/g-eV
global const volfrac_1 = chord_1 / (chord_1 + chord_2)
global const volfrac_2 = 1.0 - volfrac_1
global const sb_const = 5.6704e-5 * 11604.525^4  # erg/cm^2-s-eV^4
global const num_divs = convert(Int64, 1e4)
global const max_iterations = convert(Int64, 1e5)
global const num_say = convert(Int64, 1e2)

# Linear factors
global const factor_1 = 4.0 * arad / (dens_1 * spec_heat_1)  # eV^-3
global const factor_2 = 4.0 * arad / (dens_2 * spec_heat_2)  # eV^-3

# Initial conditions
global const init_intensity = 1e0  # erg/cm^2-s
global const init_temp = 1e0  # eV

# Iteration condition
global const tolerance = convert(Float64, eps(Float32))
