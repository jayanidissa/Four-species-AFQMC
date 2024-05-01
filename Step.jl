module Step

include("H_K.jl")
include("Initialization.jl")
include("HalfK_V.jl")
include("Measure.jl")
include("Population.jl")


using .Main: nParamsSys, ParamsRun
using .H_K: nh_k,ph_k
using .Initialization: n_ini_internal,p_ini_internal,ini_Phi_E,ini_pop,ini_aux

using .HalfK_V: halfK!,V!
using .Measure: measure
using .Population: population!

using LinearAlgebra
using MPI

function step!(nsysparams::nParamsSys, runparams::ParamsRun, nHK_up,nHK_do,pHK_up,pHK_do, nProjK_up,nProjK_do,pProjK_up,pProjK_do, PhiT, Phis, ws, Os, facnorm, AB_auxfields, AC_auxfields,AD_auxfields, BC_auxfields, BD_auxfields, CD_auxfields, Etot, Wtot, flagmeas, rng, comm, taskspercore, walk_ini, walk_fin)
    nA = nsysparams.nA
    nB = nsysparams.nB
    nC = nsysparams.nC
    nD = nsysparams.nD
    apar = nA + nB + nC + nD



    nsites = nsysparams.nx*nsysparams.ny*nsysparams.nz
    walk = runparams.nwalk


    es = zeros(Complex{Float64}, walk, 1)

    local_es = es[walk_ini:walk_fin]
    local_ws = ws[walk_ini:walk_fin]
    local_Os = Os[walk_ini:walk_fin]

    # Treating the Phis differently, since I want to
    # avoid needless copying of millions of numbers
    local_Phis = similar(Phis[:, :, walk_ini:walk_fin])

    for iwa = 1:taskspercore
        onephi = Phis[:,:,walk_ini + iwa-1]
        if real(local_ws[iwa]) > 0.
            local_ws[iwa] *= exp(facnorm )

            onephi, local_ws[iwa], local_Os[iwa],n_invO_A, n_invO_B,p_invO_C, p_invO_D = halfK!(nsysparams, onephi,local_ws[iwa], local_Os[iwa], nProjK_up,nProjK_do,pProjK_up,pProjK_do,PhiT)
            if real(local_ws[iwa]) > 0.
                for jsite = 1:nsites
                    if real(local_ws[iwa]) > 0.
                        onephi[jsite,:], local_ws[iwa], local_Os[iwa], n_invO_A, n_invO_B,p_invO_C, p_invO_D = V!(nsysparams, onephi[jsite,:], local_ws[iwa],  local_Os[iwa], n_invO_A, n_invO_B,p_invO_C, p_invO_D, PhiT[jsite,:], AB_auxfields,AC_auxfields, AD_auxfields,BC_auxfields,BD_auxfields,CD_auxfields,rng)
                    end
                end
            end

            if real(local_ws[iwa]) > 0.
                onephi, local_ws[iwa], local_Os[iwa], n_invO_A, n_invO_B,p_invO_C, p_invO_D = halfK!(nsysparams, onephi, local_ws[iwa], local_Os[iwa], nProjK_up,nProjK_do,pProjK_up,pProjK_do,PhiT)
                if real(local_ws[iwa]) > 0.
                    if flagmeas == 1
                        local_es[iwa] = measure(nsysparams::nParamsSys, nHK_up,nHK_do, pHK_up,pHK_do,onephi, n_invO_A, n_invO_B,p_invO_C, p_invO_D,PhiT)
                    end
                end
            end
        end
        local_Phis[:,:,iwa] = onephi
    end

    MPI.Allgather!(local_es, es, taskspercore, comm)
    MPI.Allgather!(local_ws, ws, taskspercore, comm)
    MPI.Allgather!(local_Os, Os, taskspercore, comm)
    MPI.Allgather!(local_Phis, Phis, nsites*apar*taskspercore, comm)
    MPI.Barrier(comm)

    if flagmeas == 1
        #println("e values")
        #display(es)

        for iwa = 1:walk
            if real(ws[iwa])> 0.
                Etot += es[iwa]*ws[iwa]
                Wtot += ws[iwa]
            end
        end
    end

    return Phis, ws, Os, Etot, Wtot
end

#sysparams = ParamsSys()
#runparams = ParamsRun()
#HK_up,HK_do, ProjK_up,ProjK_do = ini_internal(sysparams, runparams)
#PhiT, ET = ini_Phi_E(sysparams, HK_up,HK_do)
#Phis, ws, Os = ini_pop(sysparams, runparams, PhiT)
#facnorm, auxfields = ini_aux(sysparams, runparams, ET)

#display(HK_up,HK_do)
#println(" ")
#println(" ")
#println(" ")
#println(" ")
#display(Phis)

#Edum = 0.
#Wdum = 0.
#flagmeas = 0
#Phis, ws, Os, Edum, Wdum = step!(sysparams, runparams, HK_up,HK_do, ProjK_up,ProjK_do, PhiT, Phis, ws, Os, facnorm, auxfields, Edum, Wdum, flagmeas)
#flagmeas = 1
#Phis, ws, Os, Edum, Wdum = step!(sysparams, runparams, HK_up,HK_do, ProjK_up,ProjK_do, PhiT, Phis, ws, Os, facnorm, auxfields, Edum, Wdum, flagmeas)
##println(" ")
##println(" ")
#display(Phis)
#println(" ")
#println(" ")
#display(Edum)
#println(" ")
#println(" ")
#
## Changing weights by hand, in order to test Population module
#ws = [4 500 2 100.]
#display(ws)
#println(" ")
#println(" ")
#display(Os)
#println(" ")
#
#Phis, ws, Os = population!(sysparams, runparams, Phis, ws, Os)
#println(" ")
#println(" ")
#display(Phis)
#println(" ")
#println(" ")
#display(ws)
#println(" ")
#println(" ")
#display(Os)
#println(" ")

end
