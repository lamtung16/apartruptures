include("apart.jl")

signal = [
    1.0 0.1;
    1.2 6.2;
    0.8 6.1;
    3.0 2.0;
    3.2 2.1;
    2.8 2.0;
    3.1 1.9;
    1.0 3.1;
    1.1 3.0;
    0.9 3.2
]

penalty = 0.1
num_states = 5

bkps = apart(signal, penalty, num_states)
println("Changepoints: ", bkps)