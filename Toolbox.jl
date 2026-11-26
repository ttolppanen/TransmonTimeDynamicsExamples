import QuantumStates as QStates
import QuantumOperators as QOps
import QuantumTimeEvolution as QTime
import QuantumToolbox as QBox
import OrdinaryDiffEq as ODE
using Plots
using SparseArrays

# https://qutip.org/QuantumToolbox.jl/stable/
# QuantumToolbox.jl is a Julia package for simulating quantum systems. I've understood it is similar
# to QuTiP for python. Since it is an established package, you should prefer it over my implementations.
# The package can solve master equations which my solvers do not do. This is advantageous over the trajectory
# approach for small systems typically. So if you are dealing with small open systems start with the toolbox.
# In general, I would advice using the QuantumToolbox for your simulations, until you run into problems with
# system size, then you can try the krylov implementation, or with open systems the trajectory approach. Also,
# I'm not sure how to implement measurement the way we have studied it with the package, for that it might be 
# easier to use my implementations. However, even here, you can override the unitary solver part with ToolBox,
# similar how in my implementation the unitary can be done via exact diagonalization or krylov method.

# The main point of these examples is to show how the existing state and operator generators can be used
# also with the QuantumToolbox.jl. Also here will be an example of the tunable coupler.

function three_transmons()
    d = 2 # dimension of each transmon 
    L = 3 # number of transmons

    w = 1000 # frequency of each transmon, you should think these are in the units of J, w/J etc.
    U = 40 # anharmonicity of each transmon
    J = 1 # coupling between transmons

    t = 10 # total time of evolution
    t_points = 1000 # number of time steps
    t_list = LinRange(0, t, t_points)

    initial_state = QStates.productstate(d, [1,0,0]) # |100>
    H = QOps.bosehubbard(d, L; w=w, U=U, J=J) # Bose-Hubbard Hamiltonian for 3 transmons

    # Turning the state and Hamiltonian to something QuantumToolbox can deal with.
    initial_state = QBox.QuantumObject(initial_state)
    H = QBox.QuantumObject(H)
    
    n1 = QBox.QuantumObject(QOps.singlesite_n(d, L, 1)) # n * I * I
    ntot = QBox.QuantumObject(QOps.nall(d, L))
    # these are used to calculate expectation values
    e_ops = [n1, ntot]

    # Giving the solver algorithm Vern9 (from ODE.jl), the tolerances and maxiters is not necessary, but
    # I ran into some numerical error withouth those, so I give them here as reference.
    abstol = 1e-14
    reltol = 1e-8
    # mesolve refers to master equation solver, but if it is not necessary QuantumToolbox will not promote
    # the state vector to a density operator etc., it uses sesolve instead. 
    time_evolution = QBox.mesolve(H, initial_state, t_list, e_ops = e_ops; alg = ODE.Vern9(), abstol, reltol, maxiters = 10^7)

    pl = plot(t_list, real.(time_evolution.expect[1, :]), label="n1", xlabel="time", ylabel="expectation value")
    plot!(t_list, real.(time_evolution.expect[2, :]), label="ntot", xlabel="time", ylabel="expectation value")
    display(pl)
end

# Lets do the previous example, but in the total boson number subspace while still using QuantumToolBox
function three_transmons_subspace()
    d = 2 # dimension of each transmon 
    L = 3 # number of transmons

    indices = QOps.total_boson_number_subspace_indices(d, L)
    ranges, perm_mat = QOps.total_boson_number_subspace_tools(d, L)

    w = 1000 # frequency of each transmon, you should think these are in the units of J, w/J etc.
    U = 40 # anharmonicity of each transmon
    J = 1 # coupling between transmons

    t = 10 # total time of evolution
    t_points = 1000 # number of time steps
    t_list = LinRange(0, t, t_points)

    initial_state = QStates.productstate(d, [1,0,0]) # |100>
    H = QOps.bosehubbard(d, L; w=w, U=U, J=J) # Bose-Hubbard Hamiltonian for 3 transmons

    # Lets find the subpace index of the initial state and hence the simulation
    # and then take only that part and represented it for QuantumToolbox
    i_id = QOps.find_subspace(initial_state, indices)
    initial_state = QOps.subspace_split(initial_state, ranges, perm_mat)
    initial_state = initial_state[i_id]
    initial_state = QBox.QuantumObject(initial_state)

    # Can also make a function that does it
    function in_subspace(x)
        x = QOps.subspace_split(x, ranges, perm_mat)
        return QBox.QuantumObject(x[i_id])
    end

    H = in_subspace(H)
    
    n1 = QOps.singlesite_n(d, L, 1) # n * I * I
    ntot = QOps.nall(d, L)
    # these are used to calculate expectation values
    e_ops = [in_subspace(n1), in_subspace(ntot)]

    # Giving the solver algorithm Vern9 (from ODE.jl), the tolerances and maxiters is not necessary, but
    # I ran into some numerical error withouth those, so I give them here as reference.
    abstol = 1e-14
    reltol = 1e-8
    # mesolve refers to master equation solver, but if it is not necessary QuantumToolbox will not promote
    # the state vector to a density operator etc., it uses sesolve instead. 
    time_evolution = QBox.mesolve(H, initial_state, t_list, e_ops = e_ops; alg = ODE.Vern9(), abstol, reltol, maxiters = 10^7)

    pl = plot(t_list, real.(time_evolution.expect[1, :]), label="n1", xlabel="time", ylabel="expectation value")
    plot!(t_list, real.(time_evolution.expect[2, :]), label="ntot", xlabel="time", ylabel="expectation value")
    display(pl)
end

# Hamiltonian for the tunable coupler setup, where the middle transmon is the coupler.
# Note how the coupling -g[2] is negative, due to which the off-position is below in frequency.
# The syntax is e.g. w = [w1, w2, w3] = [wq1, wc, wq2]
# and g = [g1c, g2c, g_direct]
function H_tunable_coupler(d, w, U, g)
    return QOps.bosehubbard(d, 3, [(1, 2, g[1]), (2, 3, -g[2]), (1, 3, g[3])]; w, U)
end

function tunable_coupler()
    d = 2 # dimension of each transmon 
    L = 3 # number of transmons, two "qubits" and one coupler

    w = 1000 # frequency of qubit transmons
    wc = 901 # frequency of coupler transmon, CHANGE HERE to see different results!
    U = 40 # anharmonicity of each transmon
    g12 = 1 # direct coupling between transmons
    gc = 10 # coupling to coupler

    ws = [w, wc, w]
    gs = [gc, gc, g12] # note that here there is no minus sign for g2c, it is accounted for in the function!

    t = 10 # total time of evolution
    t_points = 1000 # number of time steps
    t_list = LinRange(0, t, t_points)

    # After running and checking the results, you can change d = 3, and here the state to [2, 0, 0] |200>
    initial_state = QStates.productstate(d, [1,0,0]) # |100>
    H = H_tunable_coupler(d, ws, U, gs) # Tunable coupler Hamiltonian

    # Turning the state and Hamiltonian to something QuantumToolbox can deal with.
    initial_state = QBox.QuantumObject(initial_state)
    H = QBox.QuantumObject(H)
    
    n1 = QBox.QuantumObject(QOps.singlesite_n(d, L, 1)) # n * I * I
    ntot = QBox.QuantumObject(QOps.nall(d, L))
    # these are used to calculate expectation values
    e_ops = [n1, ntot]

    # Giving the solver algorithm Vern9 (from ODE.jl), the tolerances and maxiters is not necessary, but
    # I ran into some numerical error withouth those, so I give them here as reference.
    abstol = 1e-14
    reltol = 1e-8
    # mesolve refers to master equation solver, but if it is not necessary QuantumToolbox will not promote
    # the state vector to a density operator etc., it uses sesolve instead. 
    time_evolution = QBox.mesolve(H, initial_state, t_list, e_ops = e_ops; alg = ODE.Vern9(), abstol, reltol, maxiters = 10^7)

    pl = plot(t_list, real.(time_evolution.expect[1, :]), label="n1", xlabel="time", ylabel="expectation value")
    plot!(t_list, real.(time_evolution.expect[2, :]), label="ntot", xlabel="time", ylabel="expectation value")
    display(pl)
end

# three_transmons()
# three_transmons_subspace()
tunable_coupler()