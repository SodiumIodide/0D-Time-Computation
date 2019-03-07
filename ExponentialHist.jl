module ExponentialHist
    export ExpHist, histogram, number, push, distribution

    mutable struct ExpHist
        m_length::Int64
        m_min::Float64
        m_max::Float64
        m_bin::Vector{Float64}
        m_delta::Float64

        # Internal constructor
        ExpHist(len::Int64, min::Float64, max::Float64) = new(len, log10(min), log10(max), zeros(len + 1), (log10(max) - log10(min)) / convert(Float64, len))
    end

    function histogram(exh::ExpHist)::Vector{Float64}
        return exh.m_bin
    end

    function bin_no(value::Float64, delta::Float64, d_minimum::Float64, num_bins::Int64)::Int64
        local number::Int64 = convert(Int64, floor((value + delta / 2.0 - d_minimum) / delta)) + 1

        # Fit into pre-determined number of gins: excess in either direction is assumed to be maximum or minimum
        number = (number < 1) ? 1 : number
        number = (number > num_bins + 1) ? num_bins + 1 : number

        return number
    end

    function push(exh::ExpHist, x::Float64)::Nothing
        local point::Float64 = log10(x)
        local data_bin_no::Int64 = bin_no(point, exh.m_delta, exh.m_min, exh.m_length)

        exh.m_bin[data_bin_no] += 1.0

        return nothing
    end

    function distribution(exh::ExpHist)::Vector{Float64}
        return vec([10.0^(i * exh.m_delta + exh.m_min - exh.m_delta / 2.0) for i in 1:(exh.m_length + 1)])
    end
end
