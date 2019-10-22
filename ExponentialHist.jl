module ExponentialHist
    export ExpHist, histogram, push, distribution

    mutable struct ExpHist
        m_length::Int64
        m_min::Float64
        m_max::Float64
        m_bin::Vector{Float64}
        m_delta::Float64

        set_zero_subnormals(true)

        # Internal constructor
        ExpHist(len::Int64, min_in::Float64, max_in::Float64) = @fastmath new(len, log10(min_in), log10(max_in), zeros(len + 1), (log10(max_in) - log10(min_in)) / convert(Float64, len))
    end

    function histogram(exh::ExpHist)::Vector{Float64}
        return exh.m_bin
    end

    function bin_no(value::Float64, delta::Float64, d_minimum::Float64, num_bins::Int64)::Int64
        set_zero_subnormals(true)
        local number::Int64 = @fastmath convert(Int64, (floor((value + delta / 2.0 - d_minimum) / delta))) + 1

        # Fit into pre-determined number of bins: excess in either direction is assumed to be maximum or minimum
        number = @fastmath (number < 1) ? 0 : number
        number = @fastmath (number > (num_bins + 1)) ? 0 : number

        return number
    end

    function push(exh::ExpHist, x::Float64)::Nothing
        set_zero_subnormals(true)
        local point::Float64 = @fastmath log10(x)
        local data_bin_no::Int64 = bin_no(point, exh.m_delta, exh.m_min, exh.m_length)

        if @fastmath((data_bin_no != 0) && (exh.m_length >= data_bin_no))
            @inbounds @fastmath exh.m_bin[data_bin_no] += 1.0
        end

        return nothing
    end

    function distribution(exh::ExpHist)::Vector{Float64}
        set_zero_subnormals(true)
        return @inbounds @fastmath vec([10.0^(i * exh.m_delta + exh.m_min - exh.m_delta / 2.0) for i in 1:(exh.m_length + 1)])
    end
end
