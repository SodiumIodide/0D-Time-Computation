module GeometryGen
    export get_geometry

    using Random

    function get_geometry(chord_a::Float64, chord_b::Float64, end_time::Float64, num_divs::Int64, generator::MersenneTwister=MersenneTwister(1234))::Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Float64}
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

        # Determine first material to use
        local prob_a::Float64 = chord_a / (chord_a + chord_b)
        rand_num = rand(generator, Float64)
        if (rand_num < prob_a)
            material_num = 1
        else
            material_num = 2
        end

        # Loop to build geometry
        while (cons_time < end_time)
            # Generate a random number
            rand_num = rand(generator, Float64)

            # Assign a chord length based on material number
            if (material_num == 1)
                chord = chord_a  # s
            else
                chord = chord_b  # s
            end

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
            if (material_num == 1)
                material_num = 2
            else
                material_num = 1
            end
        end

        return (x_delta, x_arr, materials, num_cells)
    end
end
