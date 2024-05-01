using MPI
MPI.Init()

include("Sample.jl")
using .Sample: nParamsSys,  ParamsRun

include("Afqmc.jl")
using .Afqmc: afqmc


# using Sample: nParamsSys, ParamsRun
# using Afqmc: afqmc

using Random
using Statistics

comm = MPI.COMM_WORLD
root = 0
const rank = MPI.Comm_rank(comm)
const ncores = MPI.Comm_size(comm)

nsysparams = nParamsSys()
runparams = ParamsRun()
rng = MersenneTwister(31415 + 11*rank)

#@show rank
#@show ncores

const taskspercore = div(runparams.nwalk, ncores)
const walk_ini = 1 + taskspercore*rank
const walk_fin = taskspercore * (rank+1)


Emeas = afqmc(nsysparams, runparams, rng, comm, taskspercore, walk_ini, walk_fin)



if MPI.Comm_rank(comm) == root
    println(" ")
    # @show(taskspercore)
    println(" ")
    #@show("Emeas",Emeas)
    println(" ")

    Ereal = real(Emeas)
    Eave = mean(Ereal)
    Eerr = std(Ereal)/sqrt(runparams.nmeasblocks)
    println(" ")
    @show(Eave, Eerr)
    println(" ")
end
