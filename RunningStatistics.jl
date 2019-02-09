module RunningStatistics
    export RunningStat, clear, push, num, mean, variance, standard_deviation

    mutable struct RunningStat
        m_n::Int64
        m_oldM::Float64
        m_newM::Float64
        m_oldS::Float64
        m_newS::Float64
    end

    function clear(r::RunningStat)::Nothing
        r.m_n = 0

        return nothing
    end

    function push(r::RunningStat, x::Float64)::Nothing
        r.m_n += 1

        # See Knuth TAOCP vol 2, 3rd edition, page 232
        if (r.m_n == 1)
            r.m_oldM = r.m_newM = x
            r.m_oldS = 0.0
        else
            r.m_newM = r.m_oldM + (x - r.m_oldM) / r.m_n
            r.m_newS = r.m_oldS + (x - r.m_oldM) * (x - r.m_newM)

            # Set up for next iteration
            r.m_oldM = r.m_newM
            r.m_oldS = r.m_newS
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
end
