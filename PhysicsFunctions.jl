module PhysicsFunctions
    include("Constants.jl")

    export sigma_a, c_v, balance_intensity, balance_temp

    function sigma_a(opacity_term::Float64, temp::Float64)::Float64
        return opacity_term / temp^3
    end

    function c_v(spec_heat_term::Float64, temp::Float64)::Float64
        return spec_heat_term
    end

    function balance_intensity(opacity::Float64, dt::Float64, past_intensity::Float64, past_temp::Float64)::Float64
        local term_1::Float64 = dt * sol^2 * opacity * arad * past_temp^4
        local term_2::Float64 = (1.0 - dt * sol * opacity) * past_intensity

        return term_1 + term_2
    end

    function balance_temp(opacity::Float64, spec_heat::Float64, density::Float64, dt::Float64, past_intensity::Float64, past_temp::Float64)::Float64
        local term_1::Float64 = dt / (density * spec_heat) * opacity * past_intensity
        local term_2::Float64 = - dt / (density * spec_heat) * sol * opacity * arad * past_temp^4
        local term_3::Float64 = past_temp

        return term_1 + term_2 + term_3
    end
end
