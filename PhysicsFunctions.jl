module PhysicsFunctions
    include("Constants.jl")

    export sigma_a, c_v, balance_intensity, balance_temp

    function sigma_a(opacity_term::Float64, temp::Float64)::Float64
        return @fastmath opacity_term / temp^3
    end

    function c_v(spec_heat_term::Float64, temp::Float64)::Float64
        return spec_heat_term
    end

    # Monte Carlo functions
    function balance_intensity(opacity::Float64, dt::Float64, past_intensity::Float64, past_temp::Float64)::Float64
        local term_1::Float64 = @fastmath dt * sol^2 * opacity * arad * past_temp^4
        local term_2::Float64 = @fastmath (1.0 - dt * sol * opacity) * past_intensity

        return @fastmath term_1 + term_2
    end

    function balance_temp(opacity::Float64, spec_heat::Float64, density::Float64, dt::Float64, past_intensity::Float64, past_temp::Float64)::Float64
        local term_1::Float64 = @fastmath dt / (density * spec_heat) * opacity * past_intensity
        local term_2::Float64 = @fastmath - dt / (density * spec_heat) * sol * opacity * arad * past_temp^4
        local term_3::Float64 = past_temp

        return @fastmath term_1 + term_2 + term_3
    end

    # Realization numerical functions
    function balance_a(intensity::Float64, temp::Float64, delta_t::Float64, opacity::Float64, past_intensity::Float64)::Float64
        local term_1::Float64 = intensity
        local term_2::Float64 = @fastmath - delta_t * sol^2 * opacity * arad * temp^4
        local term_3::Float64 = @fastmath delta_t * sol * opacity * intensity
        local term_4::Float64 = - past_intensity

        return @fastmath term_1 + term_2 + term_3 + term_4
    end

    function balance_b(intensity::Float64, temp::Float64, delta_t::Float64, opacity::Float64, spec_heat::Float64, density::Float64, past_temp::Float64)::Float64
        local term_1::Float64 = temp
        local term_2::Float64 = @fastmath - delta_t / (density * spec_heat) * opacity * intensity
        local term_3::Float64 = @fastmath delta_t / (density * spec_heat) * sol * opacity * arad * temp^4
        local term_4::Float64 = - past_temp

        return @fastmath term_1 + term_2 + term_3 + term_4
    end

    function make_jacobian(intensity::Float64, temp::Float64, delta_t::Float64, opacity_term::Float64, density::Float64, spec_heat_term::Float64)::Array{Float64, 2}
        local term_1::Float64 = @fastmath 1.0 + delta_t * sol * opacity_term / temp^3
        local term_2::Float64 = @fastmath - delta_t * sol^2 * opacity_term * arad - 3.0 * delta_t * sol * opacity_term / temp^4 * intensity
        local term_3::Float64 = @fastmath - delta_t / (density * spec_heat_term) * opacity_term / temp^3
        local term_4::Float64 = @fastmath 1.0 + 3.0 * delta_t * opacity_term / (density * spec_heat_term * temp^4) * intensity + delta_t / (density * spec_heat_term) * sol * opacity_term * arad

        return [
            term_1 term_2;
            term_3 term_4
        ]
    end

    function relative_change(delta::Vector{Float64}, original_terms::Vector{Float64})::Float64
        local current_norm::Float64 = @fastmath sqrt(sum([x^2 for x in delta]))
        local original_norm::Float64 = @fastmath sqrt(sum([x^2 for x in original_terms]))

        return @fastmath current_norm / original_norm
    end
end
