module GeometryGen
    export get_geometry

    using Random

    function get_geometry(chord_a::Float64, chord_b::Float64, end_time::Float64, num_divs::Int64; rng::MersenneTwister=MersenneTwister(1234), kwargs...)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Float64}
        # Parallel lock
        local lockpassed::Bool = haskey(kwargs, :lock)

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
        chord = (chord_a < chord_b) ? chord_a : chord_b
        local sample_time::Float64 = chord * (-log(0.5))
        local est_size::Int64 = convert(Int64, ceil(end_time / sample_time * num_divs))

        sizehint!(x_delta, est_size)
        sizehint!(materials, est_size)
        sizehint!(x_arr, est_size)

        # Determine first material to use
        local prob_a::Float64 = chord_a / (chord_a + chord_b)
        if (lockpassed)
            lock(kwargs[:lock]) do
                rand_num = rand(rng, Float64)
            end
        else
            rand_num = rand(rng, Float64)
        end
        material_num = (rand_num < prob_a) ? 1 : 2

        # Loop to build geometry
        while (cons_time < end_time)
            # Generate a random number
            if (lockpassed)
                lock(kwargs[:lock]) do
                    rand_num = rand(rng, Float64)
                end
            else
                rand_num = rand(rng, Float64)
            end

            # Assign a chord length based on material number
            chord = (material_num == 1) ? chord_a : chord_b  # s

            # Calculate and append the material length
            time = chord * (-log(rand_num))  # s
            cons_time += time  # s

            # Check on thickness to not overshoot the boundary
            if (cons_time > end_time)
                time += end_time - cons_time  # s
                cons_time = end_time  # s
            end

            # Further discretize geometry
            for i::Int64 = 1:num_divs
                push!(x_delta, (time / convert(Float64, num_divs)))
                push!(x_arr, (cons_time - time + (time / convert(Float64, num_divs) * convert(Float64, i))))
                push!(materials, material_num)
                num_cells += 1
            end

            # Update material number
            material_num = (material_num == 1) ? 2 : 1
        end

        return (x_delta, x_arr, materials, num_cells)
    end
end
