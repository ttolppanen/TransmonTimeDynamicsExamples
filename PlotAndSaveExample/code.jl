import QuantumStates as QStates
import QuantumOperators as QOps
import QuantumTimeEvolution as QTime
import QuantumToolbox as QBox
import PlotAndSave as PaS
import OrderedCollections as OrdCol
using SparseArrays


# Often you need to run a simulation, store the date, make a figure, and have the parameters visible
# in the figure. PlotAndSave does that for simple figures. It also has some functions you can look into
# for generating the parameter texts etc. if you need it for more complicated plots.

# RUN THIS FILE TO GENERATE THE PLOT AND DATA
function three_transmons()
    # we will save the parameters in a dictionary
    param = OrdCol.OrderedDict()

    param[:d] = 2 # dimension of each transmon 
    param[:L] = 3 # number of transmons

    # With this format you can have the value stored in param[:w].val, while having
    # a custom text for the parameter info on top of the plot
    param[:w] = PaS.ValueInfo(1000, "$(1000 * 5 * 1e-3)GHz") # frequency of each transmon, you should think these are in the units of J, w/J etc.
    param[:U] = PaS.ValueInfo(40, "$(40 * 5)MHz") # anharmonicity of each transmon
    param[:J] = PaS.ValueInfo(1, "$(1 * 5)MHz") # coupling between transmons

    param[:t] = 10 # total time of evolution
    param[:dt] = 0.1 # time step of evolution

    # This way it also wont save the something very large states to the hard drive
    param[:initial_state] = PaS.ValueInfo(QStates.productstate(param[:d], [1,0,0]), "|100>") # |100>
    param[:H] = "Bose-Hubbard" #Can also do it like this
    H = QOps.bosehubbard(param[:d], param[:L]; w=param[:w].val, U=param[:U].val, J=param[:J].val) # Bose-Hubbard Hamiltonian for 3 transmons
    U_op = exp(-im * param[:dt] * Matrix(H))
    
    # In practise, you might want to do d = param[:d] after setting the parameter if you
    # have parameters like d and L that are used very often...
    n1 = QOps.singlesite_n(param[:d], param[:L], 1) # n * I * I
    ntot = QOps.nall(param[:d], param[:L])
    observables = [state -> QOps.expval(state, n1), state -> QOps.expval(state, ntot)] # observables are calculates as functions, you can also add more functions to the list

    time_evolution = QTime.exactevolve(param[:initial_state].val, U_op, param[:dt], param[:t], observables...) # exact time evolution
    # display(time_evolution) 
    # The output is 2x101 matrix, where [1, :] gives you the n1 values, and [2, :] gives you the ntot values at each time step. The time steps are given by 0:dt:t.

    lines = []
    dt = param[:dt]
    t = param[:t]
    push!(lines, PaS.LineInfo(0:dt:t, time_evolution[1, :], 1, "n1")) # x-axis, y-axis, number of trajectories, label
    push!(lines, PaS.LineInfo(0:dt:t, time_evolution[2, :], 1, "ntot")) 

    path = joinpath(@__DIR__)
    PaS.makeplot(path, lines...; xlabel = "t", ylabel = "expecation values", copy_code = false, param...)
    # copy_code = false so it doesnt copy this code.jl file again here, you want to use copy_code = true if
    # saving the data to a different folder than where the code exists.
    # passing the param as here, param..., opens the dictionary into keyword argument pairs, which
    # then the function knows how to handle to generate the text.
end

using Plots
# This is an example how to access the data generate by PlotAndSave.
# Test this after running the above function.
function reading_data()
    plotinfo = PaS.loaddata(@__DIR__) # find the data.jld2 file from this folder
    println("\n Whats inside plotinfo")
    display(fieldnames(typeof(plotinfo)))
    println("\n How to access the lines of the plot")
    display(keys(plotinfo.lines))
    
    # println("\n What do the parameters contain")
    # display(keys(plotinfo.parameters))
    # display(plotinfo.parameters[:d])
    # display(plotinfo.parameters[:initial_state])

    n1_line = plotinfo.lines["n1"] # note that these match from what is printed from above
    ntot_line = plotinfo.lines["ntot"]

    pl = plot(n1_line.x, n1_line.y, label = n1_line.tag)
    plot!(ntot_line.x, ntot_line.y, label = ntot_line.tag)
    plot!(pl, xlabel = plotinfo.xlabel, ylabel = plotinfo.ylabel)
    display(pl)
end

three_transmons()
# reading_data()