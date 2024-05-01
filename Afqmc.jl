module Afqmc



include("H_K.jl")
include("Initialization.jl")
include("Step.jl")
include("Ortho.jl")
include("Population.jl")
include("Measure.jl")


using .Main: nParamsSys, ParamsRun
using .H_K: nh_k,ph_k
using .Initialization: n_ini_internal,p_ini_internal, ini_Phi_E,ini_pop,ini_aux
using .Step: step!
using .Ortho: ortho!
using .Population: population!
using .Measure: measure


function afqmc(nsysparams::nParamsSys, runparams::ParamsRun, rng, comm, taskspercore, walk_ini, walk_fin)
    nA = nsysparams.nA
    nB = nsysparams.nB
    nC = nsysparams.nC
    nD = nsysparams.nD
    dtau = runparams.dtau
    # write it in matrix format

    species = zeros(Int8, 1, 4)
    species[1] = nA
    species[2] = nB
    species[3] = nC
    species[4] = nD
    mat_species = species'* species

    if mat_species[1,2] > 0
        U_AB = nsysparams.U_AB
    else
        U_AB = 0
    end
    if mat_species[1,3] > 0
        U_AC = nsysparams.U_AC
    else
        U_AC = 0
    end
    if mat_species[1,4] > 0
        U_AD = nsysparams.U_AD
    else
        U_AD = 0
    end
    if mat_species[2,3] > 0
        U_BC = nsysparams.U_BC
    else
        U_BC = 0
    end
    if mat_species[2,4] > 0
        U_BD = nsysparams.U_BD
    else
        U_BD = 0
    end
    if mat_species[3,4] > 0
        U_CD = nsysparams.U_CD
    else
        U_CD = 0
    end


    nsites = nsysparams.nx*nsysparams.ny*nsysparams.nz
    dtau = runparams.dtau
    walk = runparams.nwalk
    nstepsperblock = runparams.nstepsperblock
    nequilblocks = runparams.nequilblocks
    nmeasblocks = runparams.nmeasblocks
    itv_ortho = runparams.itv_ortho
    itv_popu = runparams.itv_popu
    itv_meas = runparams.itv_meas


    nHK_up,nHK_do, nProjK_up, nProjK_do = n_ini_internal(nsysparams, runparams)
    pHK_up,pHK_do, pProjK_up, pProjK_do = p_ini_internal(nsysparams, runparams)
    PhiT, ET = ini_Phi_E(nsysparams, nHK_up,nHK_do,pHK_up,pHK_do)
    Phis, ws, Os = ini_pop(nsysparams, runparams, PhiT)
    facnorm,AB_auxfields,AC_auxfields,AD_auxfields,BC_auxfields,BD_auxfields,CD_auxfields = ini_aux(nsysparams, runparams, ET)


    # Equilibration stage
    flagmeas = 0
    Edum = 0.
    Wdum = 0.

    for iequil = 1:nequilblocks
        for jstep = 1:nstepsperblock
            Phis, ws, Os, Edum, Wdum = step!(nsysparams, runparams, nHK_up,nHK_do,pHK_up,pHK_do, nProjK_up,nProjK_do,pProjK_up,pProjK_do, PhiT, Phis, ws, Os, facnorm, AB_auxfields, AC_auxfields, AD_auxfields, BC_auxfields, BD_auxfields, CD_auxfields, Edum, Wdum, flagmeas, rng, comm, taskspercore, walk_ini, walk_fin)

            if (jstep % itv_ortho) == 0
                Phis, Os = ortho!(nsysparams, runparams, Phis, Os, comm, taskspercore, walk_ini, walk_fin)
            end
            if (jstep % itv_popu) == 0
                Phis, ws, Os = population!(nsysparams, runparams, Phis, ws, Os, rng)
            end
        end
    end


    # Measurement stage
    Emeas = zeros(Complex{Float64}, nmeasblocks,1)
    Wmeas = zeros(Complex{Float64}, nmeasblocks,1)

#    println("n_Wmeas")

    for imeas = 1:nmeasblocks
        for jstep = 1:nstepsperblock
            #println("imeas, jstep", imeas, " ", jstep)
            if (jstep % itv_meas) == 0
                flagmeas = 1
            else
                flagmeas = 0
                #println("check the ratio of steps_per_block/ itv_meas")
            end

            Phis, ws, Os, Emeas[imeas], Wmeas[imeas] = step!(nsysparams, runparams, nHK_up,nHK_do,pHK_up,pHK_do, nProjK_up,nProjK_do,pProjK_up,pProjK_do, PhiT, Phis, ws, Os, facnorm, AB_auxfields, AC_auxfields, AD_auxfields, BC_auxfields, BD_auxfields, CD_auxfields, Emeas[imeas], Wmeas[imeas], flagmeas, rng, comm, taskspercore, walk_ini, walk_fin)
            if (jstep % itv_ortho) == 0
                Phis, Os = ortho!(nsysparams, runparams, Phis, Os, comm, taskspercore, walk_ini, walk_fin)
                #else
                #println("check the ratio of steps_per_block/ itv_ortho")

            end
            if (jstep % itv_popu) == 0
                Phis, ws, Os = population!(nsysparams, runparams, Phis, ws, Os, rng)
                #else
                #println("check the ratio of steps_per_block/ itv_popu")

            end
            if (jstep % itv_meas) == 0

                # for repulsive interactions
                #facnorm = dtau * (real(Emeas[imeas]/Wmeas[imeas]) - 0.5*(U_AB*(nA + nB) + U_BC*(nB + nC) + U_AC*(nA + nC) + U_AD*(nA + nD) + U_BD*(nB + nD) + U_CD*(nC + nD)))

                # for attractive intercations
                facnorm = dtau * (real(Emeas[imeas]/Wmeas[imeas]))
                #else
                #println("check the ratio of steps_per_block/ itv_meas")


            end
        end

    #println(imeas, n_Wmeas[imeas])
    @show(imeas, Wmeas[imeas])

    Emeas[imeas] /= Wmeas[imeas]
    #@show(imeas,Emeas[imeas])
    end

    return Emeas
    @show Emeas
end

#sysparams = ParamsSys()
#runparams = ParamsRun()
#Emeas = afqmc(sysparams, runparams)
#println(" ")
#println(" ")
#display(Emeas)
#println(" ")
#
#Ereal = real(Emeas)
#Eave = mean(Ereal)
#Eerr = std(Ereal)/sqrt(runparams.nmeasblocks)
#println(" ")
#println(Eave, " ", Eerr)
#println(" ")

end
