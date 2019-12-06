module GeometryGenThinning
    export get_geometry_thinning

    include("NonHomogChord.jl")
    using Random
    include("Constants.jl")

    # Nonhomogeneous geometry generation
    function get_geometry_thinning(end_time::Float64, num_divs::Int64; rng::MersenneTwister=MersenneTwister(1234))::Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Float64}
        set_zero_subnormals(true)
        # Computational utilities
        local material_num::Int32
        local rand_num::Float64
        local chord::Float64

        # Updateable values
        local cons_time::Float64 = 0.0  # s
        local time::Float64 = 0.0  # s

        # Return values
        local x_delta::Vector{Float64} = Float64[]
        local materials::Vector{Int32} = Int32[]
        local x_arr::Vector{Float64} = Float64[]
        local num_cells::Int64 = 0

        # Sample size
        chord = @fastmath minimum([start_chord_1, end_chord_1, start_chord_2, end_chord_2])
        local sample_time::Float64 = @fastmath chord * (-log(0.5))
        local est_size::Int64 = convert(Int64, ceil(@fastmath (end_time / sample_time * num_divs)))

        sizehint!(x_delta, est_size)
        sizehint!(materials, est_size)
        sizehint!(x_arr, est_size)

        # Problem parameters
        # Linear model slope
        local (slope_1::Float64, slope_2::Float64) = (0.0, 0.0)
        # Quadratic model equation parameters: a*x^2 + b*x + c
        local (param_a_1::Float64, param_b_1::Float64, param_c_1::Float64) = (0.0, 0.0, 0.0)
        local (param_a_2::Float64, param_b_2::Float64, param_c_2::Float64) = (0.0, 0.0, 0.0)
        # The value that the chord possesses for the maximum value computed (for rejection purposes)
        # As exponential distribution is computed using the inverse of the chord-length, the minimum value is chosen
        local (limiting_value_chord_1::Float64, limiting_value_chord_2::Float64) = (0.0, 0.0)
        if (quad)
            (param_a_1, param_b_1, param_c_1, param_a_2, param_b_2, param_c_2) = NonHomogChord.get_quad_params(start_chord_1, end_chord_1, start_chord_2, end_chord_2, end_time)
            limiting_value_chord_1 = minimum([start_chord_1, end_chord_1, NonHomogChord.mid_value_1_func(start_chord_1, end_chord_1)])
            limiting_value_chord_2 = minimum([start_chord_2, end_chord_2, NonHomogChord.mid_value_2_func(start_chord_2, end_chord_2)])
        else
            slope_1 = @fastmath (end_chord_1 - start_chord_1) / end_time
            slope_2 = @fastmath (end_chord_2 - start_chord_2) / end_time
            limiting_value_chord_1 = min(start_chord_1, end_chord_1)
            limiting_value_chord_2 = min(start_chord_2, end_chord_2)
        end

        # Calculation variables
        local limiting_value_chord::Float64 = 0.0
        local buffer_chord::Float64 = 0.0
        local prob_accept::Float64 = 0.0
        # Linear values
        local (chord_start::Float64, chord_slope::Float64) = (0.0, 0.0)
        # Quadratic values
        local (param_a::Float64, param_b::Float64, param_c::Float64) = (0.0, 0.0, 0.0)

        # Determine first material to use
        # For nonhomogeneous chords, initial time is at 0.0 and so the probability is equivalent to the constant term ratio
        local prob_1::Float64 = @fastmath start_chord_1 / (start_chord_1 + start_chord_2)
        rand_num = @fastmath rand(rng, Float64)
        material_num = @fastmath (rand_num < prob_1) ? 1 : 2

        # Loop to build geometry
        while (cons_time < end_time)
            time = 0.0  # s
            # Assign a chord length based on material number
            if @fastmath (material_num == 1)
                if (quad)
                    param_a = param_a_1
                    param_b = param_b_1
                    param_c = param_c_1
                else
                    chord_start = start_chord_1
                    chord_slope = slope_1
                end
                limiting_value_chord = limiting_value_chord_1
            else
                if (quad)
                    param_a = param_a_2
                    param_b = param_b_2
                    param_c = param_c_2
                else
                    chord_start = start_chord_2
                    chord_slope = slope_2
                end
                limiting_value_chord = limiting_value_chord_2
            end

            # Loop for rejection sampling
            local accepted::Bool = false
            while (!accepted)
                # Generate a random number
                rand_num = @fastmath rand(rng, Float64)

                # Sample from a homogeneous distribution of intensity equal to maximum chord
                @fastmath time -= limiting_value_chord * log(rand_num)
                if @fastmath (time > end_time)
                    break
                end

                # Maximum value achieved with minimum length chord
                # Or, conversely, the inverse of the chord (therefore inverse division)
                if (quad)
                    buffer_chord = NonHomogChord.quad_chord(param_a, param_b, param_c, cons_time + time)
                else
                    buffer_chord = NonHomogChord.linear_chord(chord_start, chord_slope, cons_time + time)
                end
                prob_accept = @fastmath limiting_value_chord / buffer_chord

                rand_num = @fastmath rand(rng, Float64)
                if @fastmath (rand_num < prob_accept)
                    accepted = true
                end
            end

            @fastmath cons_time += time  # s

            # Check on time to not overshoot the boundary
            if @fastmath (cons_time > end_time)
                @fastmath time += end_time - cons_time  # s
                cons_time = end_time  # s
            end

            # Further discretize geometry
            for i::Int64 = 1:num_divs
                push!(x_delta, @fastmath(time / convert(Float64, num_divs)))
                push!(x_arr, @fastmath(cons_time - time + (time / convert(Float64, num_divs) * convert(Float64, i))))
                push!(materials, material_num)
                @fastmath num_cells += 1
            end

            # Update material number
            material_num = @fastmath (material_num == 1) ? 2 : 1
        end

        return (x_delta, x_arr, materials, num_cells)
    end
end
