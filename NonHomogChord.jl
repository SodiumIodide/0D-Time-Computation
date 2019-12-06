#!/usr/bin/env julia

module NonHomogChord
    export linear_chord, quad_chord, mid_value_1_func, mid_value_2_func, calc_c, calc_b, calc_a, get_quad_params

    @inline function linear_chord(chord_initial::Float64, chord_slope::Float64, distance::Float64)::Float64
        return @fastmath chord_initial + chord_slope * distance
    end

    @inline function quad_chord(param_a::Float64, param_b::Float64, param_c::Float64, distance::Float64)::Float64
        return @fastmath param_a * distance^2 + param_b * distance + param_c
    end

    @inline function mid_value_1_func(start_value_1::Float64, end_value_1::Float64)::Float64
        return @fastmath (start_value_1 > end_value_1) ? 2.0 * start_value_1 : 2.0 * end_value_1
    end

    @inline function mid_value_2_func(start_value_2::Float64, end_value_2::Float64)::Float64
        return @fastmath (start_value_2 > end_value_2) ? end_value_2 / 2.0 : start_value_2 / 2.0
    end

    @inline function calc_c(start_value::Float64)::Float64
        return start_value
    end

    @inline function calc_b(delta_m::Float64, delta_t::Float64, end_dist::Float64)::Float64
        return (2.0 * delta_m) / end_dist * (1.0 + sqrt(1.0 - delta_t / delta_m))
    end

    @inline function calc_a(delta_t::Float64, param_b::Float64, end_dist::Float64)::Float64
        return (delta_t - param_b * end_dist) / end_dist^2
    end

    function get_quad_params(start_value_1::Float64, end_value_1::Float64, start_value_2::Float64, end_value_2::Float64, end_dist::Float64)::Tuple{Float64, Float64, Float64, Float64, Float64, Float64}
        local mid_value_1::Float64 = mid_value_1_func(start_value_1, end_value_1)
        local mid_value_2::Float64 = mid_value_2_func(start_value_2, end_value_2)
        local delta_t_1::Float64 = @fastmath end_value_1 - start_value_1
        local delta_t_2::Float64 = @fastmath end_value_2 - start_value_2
        local delta_m_1::Float64 = @fastmath mid_value_1 - start_value_1
        local delta_m_2::Float64 = @fastmath mid_value_2 - start_value_2
        local param_c_1::Float64 = calc_c(start_value_1)
        local param_b_1::Float64 = calc_b(delta_m_1, delta_t_1, end_dist)
        local param_a_1::Float64 = calc_a(delta_t_1, param_b_1, end_dist)
        local param_c_2::Float64 = calc_c(start_value_2)
        local param_b_2::Float64 = calc_b(delta_m_2, delta_t_2, end_dist)
        local param_a_2::Float64 = calc_a(delta_t_2, param_b_2, end_dist)
        return (param_a_1, param_b_1, param_c_1, param_a_2, param_b_2, param_c_2)
    end
end
