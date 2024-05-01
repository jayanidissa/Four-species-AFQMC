module Ortho

include("H_K.jl")
include("Initialization.jl")
include("HalfK_V.jl")

using .Main: nParamsSys, ParamsRun
using .H_K: nh_k,ph_k
using .Initialization: n_ini_internal,p_ini_internal,ini_Phi_E,ini_pop,ini_aux
using .HalfK_V: halfK!,V!



using LinearAlgebra
using MPI

function ortho!(nsysparams::nParamsSys, runparams::ParamsRun, Phis, Os, comm, taskspercore, walk_ini, walk_fin)
    nA = nsysparams.nA
    nB = nsysparams.nB
    nC = nsysparams.nC
    nD = nsysparams.nD
    npar = nA + nB
    apar = nA + nB + nC + nD
    nsites = nsysparams.nx*nsysparams.ny*nsysparams.nz
    walk = runparams.nwalk

    local_Os = Os[walk_ini:walk_fin]

    # Treating the Phis differently, since I want to
    # avoid needless copying of millions of numbers
    local_Phis = similar(Phis[:, :, walk_ini:walk_fin])

    for iwa = 1:taskspercore
        Fup_A = qr(Phis[:, 1:nA, walk_ini + iwa-1])
        local_Phis[:, 1:nA, iwa] = Matrix(Fup_A.Q)

        Fdo_B = qr(Phis[:, nA+1:npar, walk_ini + iwa-1])
        local_Phis[:, nA + 1:npar, iwa] = Matrix(Fdo_B.Q)

        Fup_C = qr(Phis[:, npar + 1:npar + nC, walk_ini + iwa-1])
        local_Phis[:, npar + 1:npar + nC, iwa] = Matrix(Fup_C.Q)

        Fdo_D = qr(Phis[:, npar + nC + 1:apar, walk_ini + iwa-1])
        local_Phis[:, npar + nC + 1:apar, iwa] = Matrix(Fdo_D.Q)

        local_Os[iwa] = local_Os[iwa]/( (((det(Fup_A.R)) * (det(Fdo_B.R))) * (det(Fup_C.R) * det(Fdo_D.R))))
        #local_Os[iwa] /= (((det(Fup_A.R)) * (det(Fdo_B.R))) * (det(Fup_C.R) * det(Fdo_D.R)))
    end

    MPI.Allgather!(local_Os, Os, taskspercore, comm)
    MPI.Allgather!(local_Phis, Phis, nsites*apar*taskspercore, comm)
    MPI.Barrier(comm)

    return Phis, Os
end


#sysparams = ParamsSys()
#runparams = ParamsRun()
#HK, ProjK = ini_internal(sysparams, runparams)
#PhiT, ET = ini_Phi_E(sysparams, HK)
#Phis, ws, Os = ini_pop(sysparams, runparams, PhiT)
#facnorm, auxfields = ini_aux(sysparams, runparams, ET)
#
#Edum = 0.
#Wdum = 0.
#flagmeas = 0
#Phis, ws, Os, Edum, Wdum = step!(sysparams, runparams, HK, ProjK, PhiT, Phis, ws, Os, facnorm, auxfields, Edum, Wdum, flagmeas)
#Phis, ws, Os, Edum, Wdum = step!(sysparams, runparams, HK, ProjK, PhiT, Phis, ws, Os, facnorm, auxfields, Edum, Wdum, flagmeas)
#println(" ")
#println(" ")
#display(Phis)
#println(" ")
#println(" ")
#display(Os)
#
#Phis, Os = ortho!(sysparams, runparams, Phis, Os)
#println(" ")
#println(" ")
#display(Phis)
#println(" ")
#println(" ")
#display(Os)
#println(" ")

end
