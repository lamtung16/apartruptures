using LinearAlgebra
using Statistics

const π_ = 3.14159265358979323846

# Computes the squared geodesic distance between two angular vectors.
function sq_geo_dist(vect_a::AbstractVector, vect_b::AbstractVector)
    diff = abs.(vect_a .- vect_b)
    min_diff = min.(diff, 2π_ .- diff)
    return sum(min_diff .^ 2)
end

# Computes the circular mean for each column of a segment matrix.
function circular_mean(segment::AbstractMatrix)
    ncols = size(segment, 2)
    mean_angles = zeros(ncols)
    for j in 1:ncols
        col = segment[:, j]
        sin_sum = mean(sin.(col))
        cos_sum = mean(cos.(col))
        θ = atan(sin_sum, cos_sum)
        if θ < 0
            θ += 2π_
        end
        mean_angles[j] = θ
    end
    return mean_angles
end

# Initializes centroids by computing circular means for equally divided segments of the signal matrix.
function init_centroids(signal::AbstractMatrix, nStates::Int)
    nrows, ncols = size(signal)
    if nStates <= 0 || nrows < nStates
        throw(ArgumentError("nStates must be positive and ≤ signal length"))
    end

    centroids = zeros(nStates, ncols)
    segment_size = div(nrows, nStates)
    start_idx = 1
    for i in 1:nStates
        end_idx = start_idx + segment_size - 1
        centroids[i, :] = circular_mean(signal[start_idx:end_idx, :])
        start_idx = end_idx + 1
    end
    return centroids
end

# Tracks back changepoints from a path vector.
function track_back(path_vec::Vector{Int})
    chpnts = [length(path_vec)]
    s = path_vec[end]
    while s >= 0
        pushfirst!(chpnts, s)
        s = path_vec[s]
    end
    return chpnts
end

# apart
function apart(signal::AbstractMatrix, pen::Float64, nStates::Int)
    T, D = size(signal)
    centroids = init_centroids(signal, nStates)

    V = zeros(T + 1, nStates)
    tau = fill(-1, T + 1, nStates)
    path_vec = fill(-1, T)

    best_prev = 0.0
    for t in 1:T
        best_k = 1
        for k in 1:nStates
            V[t + 1, k] = V[t, k]
            tau[t + 1, k] = tau[t, k]
            if best_prev + pen < V[t, k]
                V[t + 1, k] = best_prev + pen
                tau[t + 1, k] = t - 1
            end
            V[t + 1, k] += sq_geo_dist(view(centroids, k, :), view(signal, t, :))
            if V[t + 1, k] < V[t + 1, best_k]
                best_k = k
            end
        end
        best_prev = V[t + 1, best_k]
        path_vec[t] = tau[t + 1, best_k]
    end

    return track_back(path_vec)
end