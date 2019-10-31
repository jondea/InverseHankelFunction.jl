
function testvalues(t::Type{CT}, removezeros=true) where CT <: Complex

    # Some special values in a grid shape and points nearby
    grid = [zero(CT), one(CT), -one(CT), one(CT)*im, -one(CT)*im]
    ε = 1e-6
    grid_with_perturbations = reshape([point1+ε*point2 for point1 in grid, point2 in grid], :)

    # Some random values with different magnitudes
    N = 10
    randvals(stddev) = randn(CT, N) .* stddev
    random_values = vcat(randvals.([1,10,100])...)

    testvalue_vec = vcat(grid_with_perturbations, random_values)
    testvalue_vec = grid_with_perturbations

    if removezeros filter!(!iszero, testvalue_vec) end

    return testvalue_vec
end
