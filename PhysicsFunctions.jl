module PhysicsFunctions
    include("Constants.jl")

    export sigma_a, c_v, balance_intensity, balance_temp, balance_a, balance_b, make_jacobian, balance_f1, balance_f2, balance_g1, balance_g1, make_jacobian_heuristic

    function sigma_a(opacity_term::Float64, temp::Float64)::Float64
        set_zero_subnormals(true)
        return @fastmath opacity_term / temp^3
    end

    function c_v(spec_heat_term::Float64, temp::Float64)::Float64
        set_zero_subnormals(true)
        return spec_heat_term
    end

    # Monte Carlo functions
    function balance_intensity(opacity::Float64, dt::Float64, past_intensity::Float64, past_temp::Float64)::Float64
        set_zero_subnormals(true)
        local term_1::Float64 = @fastmath dt * sol^2 * opacity * arad * past_temp^4
        local term_2::Float64 = @fastmath (1.0 - dt * sol * opacity) * past_intensity

        return @fastmath term_1 + term_2
    end

    function balance_temp(opacity::Float64, spec_heat::Float64, density::Float64, dt::Float64, past_intensity::Float64, past_temp::Float64)::Float64
        set_zero_subnormals(true)
        local term_1::Float64 = @fastmath dt / (density * spec_heat) * opacity * past_intensity
        local term_2::Float64 = @fastmath - dt / (density * spec_heat) * sol * opacity * arad * past_temp^4
        local term_3::Float64 = past_temp

        return @fastmath term_1 + term_2 + term_3
    end

    # Realization numerical functions
    function balance_a(intensity::Float64, temp::Float64, delta_t::Float64, opacity::Float64, past_intensity::Float64)::Float64
        set_zero_subnormals(true)
        local term_1::Float64 = intensity
        local term_2::Float64 = @fastmath - delta_t * sol^2 * opacity * arad * temp^4
        local term_3::Float64 = @fastmath delta_t * sol * opacity * intensity
        local term_4::Float64 = @fastmath - past_intensity

        return @fastmath term_1 + term_2 + term_3 + term_4
    end

    function balance_b(intensity::Float64, temp::Float64, delta_t::Float64, opacity::Float64, spec_heat::Float64, density::Float64, past_temp::Float64)::Float64
        set_zero_subnormals(true)
        local term_1::Float64 = temp
        local term_2::Float64 = @fastmath - delta_t / (density * spec_heat) * opacity * intensity
        local term_3::Float64 = @fastmath delta_t / (density * spec_heat) * sol * opacity * arad * temp^4
        local term_4::Float64 = @fastmath - past_temp

        return @fastmath term_1 + term_2 + term_3 + term_4
    end

    function make_jacobian(intensity::Float64, temp::Float64, delta_t::Float64, opacity_term::Float64, density::Float64, spec_heat_term::Float64)::Array{Float64, 2}
        set_zero_subnormals(true)
        local term_1::Float64 = @fastmath 1.0 + delta_t * sol * opacity_term / temp^3
        local term_2::Float64 = @fastmath - delta_t * sol^2 * opacity_term * arad - 3.0 * delta_t * sol * opacity_term / temp^4 * intensity
        local term_3::Float64 = @fastmath - delta_t / (density * spec_heat_term) * opacity_term / temp^3
        local term_4::Float64 = @fastmath 1.0 + 3.0 * delta_t * opacity_term / (density * spec_heat_term * temp^4) * intensity + delta_t / (density * spec_heat_term) * sol * opacity_term * arad

        return [
            term_1 term_2;
            term_3 term_4
        ]
    end

    # Heuristic model numerical functions
    function balance_f1(intensity_1::Float64, intensity_2::Float64, temp_1::Float64, temp_2::Float64, sigma_a_1::Float64, sigma_a_2::Float64, c_v_1::Float64, c_v_2::Float64, past_intensity_1::Float64)::Float64
        set_zero_subnormals(true)
        local term_1::Float64 = intensity_1
        local term_2::Float64 = @fastmath delta_t * sol * sigma_a_1 * intensity_1
        local term_3::Float64 = @fastmath - delta_t * sol^2 * sigma_a_1 * arad * temp_1^4
        local term_4::Float64 = @fastmath - delta_t * volfrac_2 / (volfrac_1 * chord_2) * intensity_2
        local term_5::Float64 = @fastmath delta_t / chord_1 * intensity_1
        local term_6::Float64 = - past_intensity_1

        return @fastmath term_1 + term_2 + term_3 + term_4 + term_5 + term_6
    end

    function balance_f2(intensity_1::Float64, intensity_2::Float64, temp_1::Float64, temp_2::Float64, sigma_a_1::Float64, sigma_a_2::Float64, c_v_1::Float64, c_v_2::Float64, past_intensity_2::Float64)::Float64
        set_zero_subnormals(true)
        local term_1::Float64 = intensity_2
        local term_2::Float64 = @fastmath delta_t * sol * sigma_a_2 * intensity_2
        local term_3::Float64 = @fastmath - delta_t * sol^2 * sigma_a_2 * arad * temp_2^4
        local term_4::Float64 = @fastmath - delta_t * volfrac_1 / (volfrac_2 * chord_1) * intensity_1
        local term_5::Float64 = @fastmath delta_t / chord_2 * intensity_2
        local term_6::Float64 = - past_intensity_2

        return @fastmath term_1 + term_2 + term_3 + term_4 + term_5 + term_6
    end

    function balance_g1(intensity_1::Float64, intensity_2::Float64, temp_1::Float64, temp_2::Float64, sigma_a_1::Float64, sigma_a_2::Float64, c_v_1::Float64, c_v_2::Float64, past_temp_1::Float64)::Float64
        set_zero_subnormals(true)
        local term_1::Float64 = temp_1
        local term_2::Float64 = @fastmath delta_t * sol / (dens_1 * c_v_1) * sigma_a_1 * arad * temp_1^4
        local term_3::Float64 = @fastmath - delta_t / (dens_1 * c_v_1) * sigma_a_1 * intensity_1
        local term_4::Float64 = @fastmath - delta_t * volfrac_2 / (volfrac_1 * chord_2) * temp_2
        local term_5::Float64 = @fastmath delta_t / chord_1 * temp_1
        local term_6::Float64 = - past_temp_1

        return @fastmath term_1 + term_2 + term_3 + term_4 + term_5 + term_6
    end

    function balance_g2(intensity_1::Float64, intensity_2::Float64, temp_1::Float64, temp_2::Float64, sigma_a_1::Float64, sigma_a_2::Float64, c_v_1::Float64, c_v_2::Float64, past_temp_2::Float64)::Float64
        set_zero_subnormals(true)
        local term_1::Float64 = temp_2
        local term_2::Float64 = @fastmath delta_t * sol / (dens_2 * c_v_2) * sigma_a_2 * arad * temp_2^4
        local term_3::Float64 = @fastmath - delta_t / (dens_2 * c_v_2) * sigma_a_2 * intensity_2
        local term_4::Float64 = @fastmath - delta_t * volfrac_1 / (volfrac_2 * chord_1) * temp_1
        local term_5::Float64 = @fastmath delta_t / chord_2 * temp_2
        local term_6::Float64 = - past_temp_2

        return @fastmath term_1 + term_2 + term_3 + term_4 + term_5 + term_6
    end

    function make_jacobian_heuristic(intensity_1::Float64, intensity_2::Float64, temp_1::Float64, temp_2::Float64)::Array{Float64, 2}
        set_zero_subnormals(true)
        local term_1::Float64 = @fastmath 1.0 + delta_t * sol * opacity_1 / temp_1^3 + delta_t / chord_1
        local term_2::Float64 = @fastmath - delta_t * volfrac_2 / (volfrac_1 * chord_2)
        local term_3::Float64 = @fastmath - delta_t * sol^2 * opacity_1 * arad - 3.0 * delta_t * sol * opacity_1 / temp_1^4 * intensity_1
        local term_4::Float64 = 0.0
        local term_5::Float64 = @fastmath - delta_t * volfrac_1 / (volfrac_2 * chord_1)
        local term_6::Float64 = @fastmath 1.0 + delta_t * sol * opacity_2 / temp_2^3 + delta_t / chord_2
        local term_7::Float64 = 0.0
        local term_8::Float64 = @fastmath - delta_t * sol^2 * opacity_2 * arad - 3.0 * delta_t * sol * opacity_2 / temp_2^4 * intensity_2
        local term_9::Float64 = @fastmath - delta_t * opacity_1 / (dens_1 * spec_heat_1 * temp_1^3)
        local term_10::Float64 = 0.0
        local term_11::Float64 = @fastmath 1.0 + 3.0 * delta_t * opacity_1 / (dens_1 * spec_heat_1 * temp_1^4) * intensity_1 + delta_t / (dens_1 * spec_heat_1) * sol * opacity_1 * arad + delta_t / chord_1
        local term_12::Float64 = @fastmath - delta_t * volfrac_2 / (volfrac_1 * chord_2)
        local term_13::Float64 = 0.0
        local term_14::Float64 = @fastmath - delta_t * opacity_2 / (dens_2 * spec_heat_2 * temp_2^3)
        local term_15::Float64 = @fastmath - delta_t * volfrac_1 / (volfrac_2 * chord_1)
        local term_16::Float64 = @fastmath 1.0 + 3.0 * delta_t * opacity_2 / (dens_2 * spec_heat_2 * temp_2^4) * intensity_2 + delta_t / (dens_2 * spec_heat_2) * sol * opacity_2 * arad + delta_t / chord_2

        return [
            term_1 term_2 term_3 term_4;
            term_5 term_6 term_7 term_8;
            term_9 term_10 term_11 term_12;
            term_13 term_14 term_15 term_16
        ]
    end

    function relative_change(delta::Vector{Float64}, original_terms::Vector{Float64})::Float64
        set_zero_subnormals(true)
        local current_norm::Float64 = @fastmath sqrt(sum([x^2 for x in delta]))
        local original_norm::Float64 = @fastmath sqrt(sum([x^2 for x in original_terms]))

        return @fastmath current_norm / original_norm
    end
end
