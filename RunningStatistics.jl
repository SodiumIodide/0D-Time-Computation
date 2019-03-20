module RunningStatistics
    export RunningStat, clear, push, num, mean, variance, standard_deviation

    using Base.Threads

    mutable struct RunningStat
        m_n::Int64
        m_oldM::Float64
        m_newM::Float64
        m_oldS::Float64
        m_newS::Float64
        m_min::Float64
        m_max::Float64

        # Internal constructor
        RunningStat() = new(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end

    function clear(r::RunningStat)::Nothing
        r.m_n = 0
        r.m_oldM = r.m_newM = r.m_oldS = r.m_newS = r.m_min = r.m_max = 0.0

        return nothing
    end

    function push(r::RunningStat, x::Float64)::Nothing
        r.m_n += 1

        # See Knuth TAOCP vol 2, 3rd edition, page 232
        if (r.m_n == 1)
            r.m_oldM = r.m_newM = x
            r.m_oldS = 0.0
            r.m_min = r.m_max = x
        else
            r.m_newM = r.m_oldM + (x - r.m_oldM) / r.m_n
            r.m_newS = r.m_oldS + (x - r.m_oldM) * (x - r.m_newM)

            # Set up for next iteration
            r.m_oldM = r.m_newM
            r.m_oldS = r.m_newS

            # Min and max
            if (x < r.m_min)
                r.m_min = x
            elseif (x > r.m_max)
                r.m_max = x
            end
        end

        return nothing
    end

    function num(r::RunningStat)::Int64
        return r.m_n
    end

    function mean(r::RunningStat)::Float64
        return (r.m_n > 0) ? r.m_newM : 0.0
    end

    function variance(r::RunningStat)::Float64
        return (r.m_n > 1) ? r.m_newS / (r.m_n - 1) : 0.0
    end

    function standard_devation(r::RunningStat)::Float64
        return sqrt(variance(r))
    end

    function least(r::RunningStat)::Float64
        return r.m_min
    end

    function greatest(r::RunningStat)::Float64
        return r.m_max
    end

    # Helper functions for threaded arrays
    function threadarray(length::Int64)::Array{RunningStatistics.RunningStat, 2}
        return [RunningStatistics.RunningStat() for i in 1:length, j in 1:nthreads()]
    end

    function total(mat_r::Array{RunningStat, 2})::Vector{Float64}
        return vec(sum((@. convert(Float64, num(mat_r))), dims=2))
    end

    function compute_mean(mat_r::Array{RunningStat, 2}, num_vec::Vector{Float64})::Vector{Float64}
        return vec(sum((@. mean(mat_r) * convert(Float64, num(mat_r))), dims=2) ./ num_vec)
    end

    function compute_variance(mat_r::Array{RunningStat, 2}, num_vec::Vector{Float64}, mean_vec::Vector{Float64})::Vector{Float64}
        local prefix::Vector{Float64} = vec(@. (num_vec - 1.0)^(-1))
        local first_sum::Vector{Float64} = vec(sum((@. (convert(Float64, num(mat_r)) - 1.0) * variance(mat_r)), dims=2))
        local second_sum::Vector{Float64} = vec(sum((@. convert(Float64, num(mat_r)) * (mean(mat_r) - mean_vec)^2), dims=2))

        return @. prefix * (first_sum + second_sum)
    end

    function compute_min(mat_r::Array{RunningStat, 2})::Vector{Float64}
        return vec(minimum((@. least(mat_r)), dims=2))
    end

    function compute_max(mat_r::Array{RunningStat, 2})::Vector{Float64}
        return vec(maximum((@. greatest(mat_r)), dims=2))
    end
end
