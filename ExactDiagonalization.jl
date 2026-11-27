import QuantumStates as QStates
import QuantumOperators as QOps
import QuantumTimeEvolution as QTime
import QuantumToolbox as QBox
using Plots
using SparseArrays

# Exact diagonalization solves the time-evolution by calculating the time-evolution operator U(dt) = exp(-i H dt) and applying it.
# It is suitable for small systems, since calculating the exp(-i H dt) involves diagonalizing the Hamiltonian, which scales poorly with system size (even though it is done only once in start).

function three_transmons()
    d = 2 # dimension of each transmon 
    L = 3 # number of transmons

    w = 1000 # frequency of each transmon, you should think these are in the units of J, w/J etc.
    U = 40 # anharmonicity of each transmon
    J = 1 # coupling between transmons

    t = 10 # total time of evolution
    dt = 0.1 # time step of evolution

    initial_state = QStates.productstate(d, [1,0,0]) # |100>
    H = QOps.bosehubbard(d, L; w=w, U=U, J=J) # Bose-Hubbard Hamiltonian for 3 transmons
    U_op = exp(-im * dt * Matrix(H))
    
    n1 = QOps.singlesite_n(d, L, 1) # n * I * I
    ntot = QOps.nall(d, L)
    observables = [state -> QOps.expval(state, n1), state -> QOps.expval(state, ntot)] # observables are calculates as functions, you can also add more functions to the list

    time_evolution = QTime.exactevolve(initial_state, U_op, dt, t, observables...) # exact time evolution
    # display(time_evolution) 
    # The output is 2x101 matrix, where [1, :] gives you the n1 values, and [2, :] gives you the ntot values at each time step. The time steps are given by 0:dt:t.

    pl = plot(0:dt:t, time_evolution[1, :], label="n1", xlabel="time", ylabel="expectation value")
    plot!(0:dt:t, time_evolution[2, :], label="ntot", xlabel="time", ylabel="expectation value")
    display(pl)
end

# The code is pretty much the same as above.
# You can run this many times to see different results.
function three_transmons_measured()
    d = 2
    L = 3

    w = 1000 
    U = 40
    J = 1

    t = 10
    dt = 0.1
    
    msr_prob = 0.1 # probability of measuring a transmon in a time-step. The amount of measurements then depends on dt.
    op_to_measure = QOps.nop(d) # we will do measurements of the excitation number.
    # measurement operators for each transmon. The structure is
    # msr_ops[site][measurement_outcome]
    # so for example msr_ops[2][1] is the operator I * |0><0| * I
    msr_ops = QOps.measurementoperators(op_to_measure, L)
    # this is a function where in each time-step a random transmon can be measured with probability msr_prob.
    meffect!(state) = QTime.random_measurement!(state, msr_ops, msr_prob)

    initial_state = QStates.productstate(d, [1,0,0]) # |100>
    H = QOps.bosehubbard(d, L; w=w, U=U, J=J)
    U_op = exp(-im * dt * Matrix(H))
    
    n1 = QOps.singlesite_n(d, L, 1)
    ntot = QOps.nall(d, L)
    observables = [state -> QOps.expval(state, n1), state -> QOps.expval(state, ntot)]

    # the measurements are passed as a function with one argument, the state. This is the effect! keyword argument
    # You can pass also your own function, as long as it takes a single argument, the state.
    # You can also pass many functions, or effects, if you wrap them in a single function,
    # i.e. function effect!(state) = func1(state); func2(state); etc. end
    time_evolution = QTime.exactevolve(initial_state, U_op, dt, t, observables...; effect! = meffect!) # exact time evolution
    # display(time_evolution) 
    # The output is 2x101 matrix, where [1, :] gives you the n1 values, and [2, :] gives you the ntot values at each time step. The time steps are given by 0:dt:t.

    pl = plot(0:dt:t, time_evolution[1, :], label="n1", xlabel="time", ylabel="expectation value")
    plot!(0:dt:t, time_evolution[2, :], label="ntot", xlabel="time", ylabel="expectation value")
    display(pl)
end

# for open systems we usually look at average values, so this function calculates the mean over trajectories. 
function traj_mean(result)
        out = zeros(size(result[1]))
        for traj in result
            out .+= traj
        end
    return out ./ length(result)
end

function three_transmons_open()
    d = 2
    L = 3

    w = 1000
    U = 40
    J = 1
    gamma = 1 # dissipation rate
    n_traj = 100 # number of trajectories to average over

    t = 10
    dt = 0.1

    initial_state = QStates.productstate(d, [1,0,0]) # |100>
    H = QOps.bosehubbard(d, L; w=w, U=U, J=J)
    # linblad operators for the system, we will have just dissipation.
    lb_ops = sqrt(gamma) * [QOps.singlesite_a(d, L, 1), QOps.singlesite_a(d, L, 2), QOps.singlesite_a(d, L, 3)]
    for op in lb_ops
        H .-= 1.0im / 2.0 .* (op' * op) # add the non-hermitian part to the Hamiltonian
    end
    time_evolve_op = exp(-im * dt * Matrix(H))
    
    n1 = QOps.singlesite_n(d, L, 1) # n * I * I
    ntot = QOps.nall(d, L)
    observables = [state -> QOps.expval(state, n1), state -> QOps.expval(state, ntot)]

    # you can also do more inside the function, for example you can add a random measurement effect to the trajectories.
    function solve_single_trajectory()
        QTime.exactevolve(initial_state, time_evolve_op, lb_ops, dt, t, observables...)
    end
    result = QTime.solvetrajectories(solve_single_trajectory, n_traj; paral = :threads) # solves n_traj trajectories in parallel using threads, can also use paral = :none
    result = traj_mean(result) # average over trajectories

    pl = plot(0:dt:t, result[1, :], label="n1", xlabel="time", ylabel="expectation value")
    plot!(0:dt:t, result[2, :], label="ntot", xlabel="time", ylabel="expectation value")
    display(pl)
end

function three_transmons_disordered()
    d = 2
    L = 3

    w = 1000
    dw = 10 # disorder strength, the disorder is generated as w_i = w + [-dw, dw]
    U = 40
    J = 1
    n_traj = 1000 # number of trajectories to average over

    t = 10
    dt = 0.1

    initial_state = QStates.productstate(d, [1,0,0]) # |100>
    
    n1 = QOps.singlesite_n(d, L, 1) # n * I * I
    ntot = QOps.nall(d, L)
    observables = [state -> QOps.expval(state, n1), state -> QOps.expval(state, ntot)]

    # in each trajectory we need to generate a new disordered Hamiltonian, so we do that inside the function.
    function solve_single_trajectory()
        w = (2 .* rand(L) .- 1) * dw .+ w # random disorder, w_i = w + [-dw, dw]
        H = QOps.bosehubbard(d, L; w=w, U=U, J=J) # note that now w is a vector, and the function can handle that.
        time_evolve_op = exp(-im * dt * Matrix(H))
        QTime.exactevolve(initial_state, time_evolve_op, dt, t, observables...)
    end
    result = QTime.solvetrajectories(solve_single_trajectory, n_traj; paral = :threads) # solves n_traj trajectories in parallel using threads, can also use paral = :none
    result = traj_mean(result) # average over trajectories

    pl = plot(0:dt:t, result[1, :], label="n1", xlabel="time", ylabel="expectation value")
    plot!(0:dt:t, result[2, :], label="ntot", xlabel="time", ylabel="expectation value")
    display(pl)
end


# These examples do the same as above, but solve the problem in total excitation number N_tot subspaces.
function three_transmons_subspace()
    d = 2
    L = 3

    w = 1000
    U = 40
    J = 1

    t = 10
    dt = 0.1

    # gives you the indices of basis states in different subspaces, e.g. indices[1] = [1]
    # since the first subspace is N_tot = 0, which has only one basis state |000>.
    # if n2_inds = indices[2], then H[n2_inds, n2_inds] gives you the Hamiltonian in the N_tot = 1 subspace.
    indices = QOps.total_boson_number_subspace_indices(d, L)
    # display(indices) uncomment here and check

    initial_state = QStates.productstate(d, [1,0,0]) # |100>
    initial_id = QOps.find_subspace(initial_state, indices) # finds the subspace id of the initial state
    H = QOps.bosehubbard(d, L; w=w, U=U, J=J) # Bose-Hubbard Hamiltonian for 3 transmons

    # note that with this approach I build the full Hamiltonian, then split it. It would be more efficient to build
    # the Hamiltonian directly in the subspaces, but in my simulations it was not the bottleneck, and with the 
    # krylov you avoid it anyway. There was an implementation of building it only in the subspace from a previous
    # group member, but there was some issue with it, I think it didn't scale well either when going to large systems.
    # Anyways, if building the Hamiltonian is a problem, consider looking into how to build it directly in the subspace.
    U_op = exp(-im * dt * Matrix(H))

    # initial_state represented in all the subpaces. initial_state[initial_id] gives the state in the correct subspace.
    # So it is a list of states. It uses more memory than necessary, but avoids reallocating memory if changing subspaces.
    initial_state = QOps.subspace_split(initial_state, indices)
    U_op = QOps.subspace_split(U_op, indices) # U_op is now a list of operators, one for each subspace. U_op[initial_id] gives the operator in the correct subspace.

    n1 = QOps.subspace_split(QOps.singlesite_n(d, L, 1), indices) # n * I * I
    ntot = QOps.subspace_split(QOps.nall(d, L), indices)
    observables = [(state, id) -> QOps.expval(state[id], n1[id]), (state, id) -> QOps.expval(state[id], ntot[id])] # note the arguments now contain the subspace id.

    time_evolution = QTime.exactevolve(initial_state, initial_id, U_op, dt, t, observables...) # note initial_id argument.

    pl = plot(0:dt:t, time_evolution[1, :], label="n1", xlabel="time", ylabel="expectation value")
    plot!(0:dt:t, time_evolution[2, :], label="ntot", xlabel="time", ylabel="expectation value")
    display(pl)
end

function three_transmons_measured_subspace()
    d = 2
    L = 3

    indices = QOps.total_boson_number_subspace_indices(d, L)
    ranges, perm_mat = QOps.total_boson_number_subspace_tools(d, L) # these give tools some functions need for the subspace representation.

    w = 1000 
    U = 40
    J = 1

    t = 10
    dt = 0.1
    
    msr_prob = 0.1
    op_to_measure = QOps.nop(d)
    msr_ops = QOps.measurementoperators(op_to_measure, L)

    msr_ops_sb = QOps.measurement_subspace(msr_ops, ranges, perm_mat) # measurement operators in the subspaces
    meffect!(state, id) = QTime.random_measurement!(state, id, msr_ops_sb, msr_prob) # note how now the effect has two arguments, state and subspace id.

    initial_state = QStates.productstate(d, [1,0,0]) # |100>
    initial_id = QOps.find_subspace(initial_state, indices) # finds the subspace of the initial state
    H = QOps.bosehubbard(d, L; w=w, U=U, J=J)
    U_op = exp(-im * dt * Matrix(H))


    initial_state = QOps.subspace_split(initial_state, indices)
    U_op = QOps.subspace_split(U_op, indices) # U_op is now a list of operators, one for each subspace. U_op[initial_id] gives the operator in the correct subspace.

    n1 = QOps.subspace_split(QOps.singlesite_n(d, L, 1), indices) # n * I * I
    ntot = QOps.subspace_split(QOps.nall(d, L), indices)
    observables = [(state, id) -> QOps.expval(state[id], n1[id]), (state, id) -> QOps.expval(state[id], ntot[id])] # note the arguments now contain the subspace id.

    time_evolution = QTime.exactevolve(initial_state, initial_id, U_op, dt, t, observables...; effect! = meffect!) # note initial_id argument.
    # display(time_evolution) 
    # The output is 2x101 matrix, where [1, :] gives you the n1 values, and [2, :] gives you the ntot values at each time step. The time steps are given by 0:dt:t.

    pl = plot(0:dt:t, time_evolution[1, :], label="n1", xlabel="time", ylabel="expectation value")
    plot!(0:dt:t, time_evolution[2, :], label="ntot", xlabel="time", ylabel="expectation value")
    display(pl)
end

function three_transmons_open_subspace()
    # I did not write a way to do the open systems in the subspaces with
    # the exact diagonalization. For krylov it is done. If you need it, 
    # I think you can take example from the krylov implementation in
    # QuantumTimeEvolution.jl -> src/KrylovSubspaceDissipationDecoherence.jl
end

function three_transmons_disordered_subspace()
    d = 2
    L = 3

    indices = QOps.total_boson_number_subspace_indices(d, L)
    ranges, perm_mat = QOps.total_boson_number_subspace_tools(d, L) 

    w = 1000
    dw = 10 # disorder strength, the disorder is generated as w_i = w + [-dw, dw]
    U = 40
    J = 1
    n_traj = 1000 # number of trajectories to average over

    t = 10
    dt = 0.1

    initial_state = QStates.productstate(d, [1,0,0]) # |100>
    initial_id = QOps.find_subspace(initial_state, indices) # finds the subspace of the initial state
    initial_state = QOps.subspace_split(initial_state, indices)
    
    n1 = QOps.subspace_split(QOps.singlesite_n(d, L, 1), indices) # n * I * I
    ntot = QOps.subspace_split(QOps.nall(d, L), indices)
    observables = [(state, id) -> QOps.expval(state[id], n1[id]), (state, id) -> QOps.expval(state[id], ntot[id])] # note the arguments now contain the subspace id.

    # in each trajectory we need to generate a new disordered Hamiltonian, so we do that inside the function.
    function solve_single_trajectory()
        w = (2 .* rand(L) .- 1) * dw .+ w # random disorder, w_i = w + [-dw, dw]
        H = QOps.bosehubbard(d, L; w=w, U=U, J=J) # note that now w is a vector, and the function can handle that.
        time_evolve_op = exp(-im * dt * Matrix(H))
        time_evolve_op = QOps.subspace_split(time_evolve_op, ranges, perm_mat) # You can use indices or ranges, perm_mat, but I think with ranges its more efficient.
        QTime.exactevolve(initial_state, initial_id, time_evolve_op, dt, t, observables...)
    end
    result = QTime.solvetrajectories(solve_single_trajectory, n_traj; paral = :threads) # solves n_traj trajectories in parallel using threads, can also use paral = :none
    result = traj_mean(result) # average over trajectories

    pl = plot(0:dt:t, result[1, :], label="n1", xlabel="time", ylabel="expectation value")
    plot!(0:dt:t, result[2, :], label="ntot", xlabel="time", ylabel="expectation value")
    display(pl)
end

three_transmons()
# three_transmons_measured()
# three_transmons_open()
# three_transmons_disordered()
# three_transmons_subspace()
# three_transmons_measured_subspace()
# three_transmons_open_subspace()
# three_transmons_disordered_subspace()