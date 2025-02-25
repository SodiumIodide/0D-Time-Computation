#!/usr/bin/env julia

# Histogram and problem parameters
global const num_t = convert(Int64, 1e5)
global const num_t_hist = convert(Int64, 1e5)
global const num_divs = convert(Int64, 1e3)
global const num_divs_hist = convert(Int64, 1e3)
global const max_iterations = convert(Int64, 1e5)
global const max_iterations_hist = convert(Int64, 1e7)

global const ss_fix_time = 1e0  # ct = cm

# Initial conditions
global const init_intensity = 1e0  # erg/cm^2-s
global const init_temp = 1e0  # eV

# Iteration condition
#global const tolerance = convert(Float64, eps(Float32))
global const tolerance = 1e-3
global const heur_tolerance = 1e-6
global const num_say = convert(Int64, 1e3)
global const num_say_hist = convert(Int64, 1e4)

# Number of bins to produce in PDF
global const num_bins = convert(Int64, 3e1)

# Histogram early time index number
global const hist_early_index = 100

# Physics constants
global const sol = 29979245800.0  # cm/s
global const arad = 137.2017050419133  # erg/cm^3-eV^4
global const sb_const = 5.6704e-5 * 11604.525^4  # erg/cm^2-s-eV^4

# Problem physics parameters
global const t_max = 1e-11  # s
global const t_init = 0.0  # s
global const delta_t = (t_max - t_init) / num_t  # s
global const chord_1 = 1e-3 / sol  # s
global const chord_2 = 5e-3 / sol  # s
global const dens_1 = 1e0  # g/cm^3
global const dens_2 = 1e0  # g/cm^3
global const opacity_1 = 1.0  # cm^-1
global const opacity_2 = 5.0  # cm^-1
global const spec_heat_1 = 1e0  # erg/g-eV
global const spec_heat_2 = 1e0  # erg/g-eV
global const volfrac_1 = chord_1 / (chord_1 + chord_2)
global const volfrac_2 = 1.0 - volfrac_1
global const destruction_factor = 1.0
# Nonhomogeneous chord lengths
global const start_chord_1 = 101.0 / 20.0 / sol * 1e-3  # s
global const start_chord_2 = 101.0 / 20.0 / sol * 1e-3  # s
global const end_chord_1 = 99.0 / 10.0 / sol * 1e-3  # s
global const end_chord_2 = 11.0 / 10.0 / sol * 1e-3  # s
global const quad = false  # quadratic or linear models for nonhomogeneous system

# Linear factors
global const factor_1 = 4.0 * arad / (dens_1 * spec_heat_1)  # eV^-3
global const factor_2 = 4.0 * arad / (dens_2 * spec_heat_2)  # eV^-3
