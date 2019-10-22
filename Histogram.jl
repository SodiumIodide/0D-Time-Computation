module Histogram
    export Hist, histogram, push, distribution

    mutable struct Hist
        m_length::Int64
        m_min::Float64
        m_max::Float64
        m_bin::Vector{Float64}
        m_delta::Float64

        set_zero_subnormals(true)

        # Internal constructor
        Hist(len::Int64, min_in::Float64, max_in::Float64) = @fastmath new(len, min_in, max_in, zeros(len + 1), (max_in - min_in) / convert(Float64, len))
    end

    function histogram(his::Hist)::Vector{Float64}
        return his.m_bin
    end

    function bin_no(value::Float64, delta::Float64, d_minimum::Float64, num_bins::Int64)::Int64
        set_zero_subnormals(true)
        local number::Int64 = @fastmath convert(Int64, (floor((value + delta / 2.0 - d_minimum) / delta))) + 1

        # Fit into pre-determined number of bins: excess in either direction is assumed to be null data
        number = @fastmath (number < 1) ? 0 : number
        number = @fastmath (number > (num_bins + 1)) ? 0 : number

        return number
    end

    function push(his::Hist, x::Float64)::Nothing
        set_zero_subnormals(true)
        local data_bin_no::Int64 = bin_no(x, his.m_delta, his.m_min, his.m_length)

        if @fastmath((data_bin_no != 0) && (his.m_length >= data_bin_no))
            @inbounds @fastmath his.m_bin[data_bin_no] += 1.0
        end

        return nothing
    end

    function distribution(his::Hist)::Vector{Float64}
        set_zero_subnormals(true)
        return @inbounds @fastmath vec([i * his.m_delta + his.m_min - his.m_delta / 2.0 for i in 1:(his.m_length + 1)])
    end
end
