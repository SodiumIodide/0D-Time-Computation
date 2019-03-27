module Histogram
    export Hist, histogram, number, push, distribution

    mutable struct Hist
        m_length::Int64
        m_min::Float64
        m_max::Float64
        m_bin::Vector{Float64}
        m_delta::Float64

        # Internal constructor
        Hist(len::Int64, min::Float64, max::Float64) = @fastmath new(len, min, max, zeros(len + 1), (max - min) / convert(Float64, len))
    end

    function histogram(his::Hist)::Vector{Float64}
        return his.m_bin
    end

    function bin_no(value::Float64, delta::Float64, d_minimum::Float64, num_bins::Int64)::Int64
        local number::Int64 = convert(Int64, @fastmath(floor((value + delta / 2.0 - d_minimum) / delta))) + 1

        # Fit into pre-determined number of gins: excess in either direction is assumed to be maximum or minimum
        number = (number < 1) ? 0 : number
        number = (number > @fastmath(num_bins + 1)) ? 0 : number

        return number
    end

    function push(his::Hist, x::Float64)::Nothing
        local data_bin_no::Int64 = bin_no(x, his.m_delta, his.m_min, his.m_length)

        if ((data_bin_no != 0) && (his.m_length >= data_bin_no))
            @inbounds @fastmath his.m_bin[data_bin_no] += 1.0
        end

        return nothing
    end

    function distribution(his::Hist)::Vector{Float64}
        return @inbounds @fastmath vec([i * his.m_delta + his.m_min - his.m_delta / 2.0 for i in 1:(his.m_length + 1)])
    end
end
