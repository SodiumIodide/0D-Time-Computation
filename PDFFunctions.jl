module PDFFunctions
    include("Constants.jl")

    export locate_steady_state, steady_state_fix

    function locate_steady_state(data::Vector{Float64})::Int64
        local point::Float64 = 0.0
        local last_point::Float64 = 0.0
        local sec_last_point::Float64 = 0.0
        local index::Int64 = 0
        local last_index::Int64 = num_t
        local steady_state_search::Bool = true

        function frac_diff(point::Float64, pointn::Float64)::Float64
            return @fastmath abs(pointn - point) / point
        end

        # Naive pass through data - three same points in a row = steady state
        while (steady_state_search)
            @fastmath index += 1
            @inbounds (point, last_point, sec_last_point) = (data[index], point, last_point)

            if @fastmath((frac_diff(sec_last_point, point) <= 0.001) && (frac_diff(last_point, point) <= 0.001))
                last_index = index
                steady_state_search = false
            end
        end

        return last_index
    end

    function steady_state_fix(times::Vector{Float64})::Int64
        local index::Int64 = 0
        local steady_state_search::Bool = true
        local max_index::Int64 = length(times)

        while (steady_state_search)
            @fastmath index += 1
            if @fastmath ((times[index] > ss_fix_time) || (index >= max_index))
                steady_state_search = false
            end
        end

        return index
    end
end
