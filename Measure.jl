module Measure

using .Main: nParamsSys, ParamsRun


using LinearAlgebra

function measure(nsysparams::nParamsSys, nHK_up,nHK_do, pHK_up,pHK_do,onephi, n_invO_A, n_invO_B,p_invO_C, p_invO_D,PhiT)
    nA = nsysparams.nA
    nB = nsysparams.nB
    nC = nsysparams.nC
    nD = nsysparams.nD
    npar = nA + nB
    apar = npar + nC + nD
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


    n_As = onephi[:, 1:nA]*n_invO_A
    n_Bs = onephi[:, nA+1:npar]*n_invO_B
    n_Cs = onephi[:,npar+1:npar + nC]*p_invO_C
    n_Ds = onephi[:,npar + nC + 1:apar]*p_invO_D
    n_GA = n_As*PhiT[:, 1:nA]'
    n_GB = n_Bs*PhiT[:, nA+1:npar]'
    n_GC = n_Cs*PhiT[:, npar + 1:npar + nC]'
    n_GD = n_Ds*PhiT[:, npar + nC + 1:apar]'

    AB_potfac = sum(diag(n_GA) .* diag(n_GB))
    AB_epot = U_AB*AB_potfac

    AC_potfac = sum(diag(n_GA) .* diag(n_GC))
    AC_epot = U_AC*AC_potfac

    AD_potfac = sum(diag(n_GA) .* diag(n_GD))
    AD_epot = U_AD*AD_potfac

    BC_potfac = sum(diag(n_GB) .* diag(n_GC))
    BC_epot = U_BC*BC_potfac

    BD_potfac = sum(diag(n_GB) .* diag(n_GD))
    BD_epot = U_BD*BD_potfac

    CD_potfac = sum(diag(n_GC) .* diag(n_GD))
    CD_epot = U_CD*CD_potfac

    epot = AB_epot + AC_epot + AD_epot + BC_epot + BD_epot + CD_epot

    interm_A = transpose(nHK_up) .* (n_GA)
    interm_B = transpose(nHK_do) .* (n_GB)
    interm_C = transpose(pHK_up) .* (n_GC)
    interm_D = transpose(pHK_do) .* (n_GD)
    ekin = sum(interm_A) + sum(interm_B) + sum(interm_C) + sum(interm_D)
    e = ekin + epot

    #e = ekin_up + ekin_do + epot
    #@show ekin,epot,e
    #@show ekin_up,ekin_do,epot,e
    return e
end


end

