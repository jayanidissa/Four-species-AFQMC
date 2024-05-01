module Initialization

include("H_K.jl")

using .Main: nParamsSys, ParamsRun
using .H_K: nh_k,ph_k

using LinearAlgebra
using Random

function n_ini_internal(nsysparams::nParamsSys, runparams::ParamsRun)
    nHK_up,nHK_do = nh_k(nsysparams)
    nProjK_up = exp(-0.5*runparams.dtau*(nHK_up))
    nProjK_do = exp(-0.5*runparams.dtau*(nHK_do))
    return nHK_up,nHK_do, nProjK_up, nProjK_do
end

function p_ini_internal(nsysparams::nParamsSys, runparams::ParamsRun)
    pHK_up,pHK_do = ph_k(nsysparams)
    pProjK_up = exp(-0.5*runparams.dtau*(pHK_up))
    pProjK_do = exp(-0.5*runparams.dtau*(pHK_do))
    return pHK_up,pHK_do, pProjK_up, pProjK_do
end



function ini_Phi_E(nsysparams::nParamsSys, nHK_up,nHK_do,pHK_up,pHK_do)

    nA = nsysparams.nA
    nB = nsysparams.nB
    nC = nsysparams.nC
    nD = nsysparams.nD

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
        #@show U_AB
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

    npar = nA + nB
    ppar = nC + nD

    n_E_NI_up = eigvals(nHK_up)
    n_E_NI_do = eigvals(nHK_do)

    n_psi_NI_up = eigvecs(nHK_up)
    n_psi_NI_do = eigvecs(nHK_do)

    n_PhiT = hcat(n_psi_NI_up[:,1:nA], n_psi_NI_do[:,1:nB])

    n_intup = diag(n_PhiT[:,1:nA]*(n_PhiT[:,1:nA])')
    n_intdo = diag(n_PhiT[:,nA+1:npar]*(n_PhiT[:,nA+1:npar])')


    p_E_NI_up = eigvals(pHK_up)
    p_E_NI_do = eigvals(pHK_do)


    AB_EK = sum(n_E_NI_up[1:nA]) + sum(n_E_NI_do[1:nB])
    AC_EK = sum(n_E_NI_up[1:nA]) + sum(p_E_NI_do[1:nC])
    AD_EK = sum(n_E_NI_up[1:nA]) + sum(p_E_NI_do[1:nD])
    CD_EK = sum(p_E_NI_up[1:nC]) + sum(p_E_NI_do[1:nD])
    BC_EK = sum(n_E_NI_do[1:nB]) + sum(p_E_NI_up[1:nC])
    BD_EK = sum(n_E_NI_do[1:nB]) + sum(p_E_NI_do[1:nD])

    p_psi_NI_up = eigvecs(pHK_up)
    p_psi_NI_do = eigvecs(pHK_do)

    p_PhiT = hcat(p_psi_NI_up[:,1:nC], p_psi_NI_do[:,1:nD])

    p_intup = diag(p_PhiT[:,1:nC]*(p_PhiT[:,1:nC])')
    p_intdo = diag(p_PhiT[:,nC+1:ppar]*(p_PhiT[:,nC+1:ppar])')

    PhiT = hcat(n_PhiT,p_PhiT)
    # @show size(PhiT)
    AB_EV = nsysparams.U_AB*(n_intup')*n_intdo
    AC_EV = nsysparams.U_AC*(n_intup')*p_intup
    AD_EV = nsysparams.U_AD*(n_intup')*p_intdo
    BD_EV = nsysparams.U_BD*(n_intdo')*p_intdo
    BC_EV = nsysparams.U_BC*(n_intdo')*p_intup
    CD_EV = nsysparams.U_CD*(p_intup')*p_intdo

    ET = AB_EK + CD_EK + AB_EV + BC_EV + CD_EV + AC_EV + AD_EV + BD_EV

    return PhiT, ET
end

function ini_pop(nsysparams::nParamsSys, runparams::ParamsRun, PhiT)
    nsites = nsysparams.nx*nsysparams.ny*nsysparams.nz
    nA = nsysparams.nA
    nB = nsysparams.nB
    nC = nsysparams.nC
    nD = nsysparams.nD
    apar = nA + nB + nC + nD
    nwalk = runparams.nwalk

    Phis = zeros(Complex{Float64}, nsites, apar, nwalk)
    for i = 1:nwalk
        Phis[:,:,i] = PhiT
    end

    ws = ones(Complex{Float64}, nwalk,1)
    Os = ones(Complex{Float64}, nwalk,1)

    return Phis, ws, Os
end


########## neutron up down interaction
function ini_aux(nsysparams::nParamsSys, runparams::ParamsRun, ET)
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
    #for repulsive attractions
    #facnorm = (real(ET) - 0.5*(U_AB*(nA + nB) + U_BC*(nB + nC) + U_AC*(nA + nC) + U_AD*(nA + nD) + U_BD*(nB + nD) + U_CD*(nC + nD)))*dtau
    #for attractive attraections
    facnorm = (real(ET))*dtau
    #@show facnorm
    if U_AB > 0

        AB_gamma = acosh(exp(0.5*dtau*abs(U_AB)))

        AB_auxfields = zeros(2,2)
        for i = 1:2
            for j = 1:2
                AB_auxfields[i,j] = exp(AB_gamma*(-1)^(i+j))
            end
        end

    elseif U_AB <= 0
    #facnorm = (real(ET) - 0.5*U*(npar-1))*dtau - this will blow the weights exponentially
        # AB_facnorm = (real(ET))*dtau
        AB_gamma = acosh(exp(0.5*dtau*abs(U_AB)))
        AB_auxfields = zeros(2,2)
        for i = 1:2
            for j = 1:2
                AB_auxfields[i,j] = exp((-1)*dtau*U_AB*0.5)*exp(AB_gamma*(-1)^(1+j))
            end
        end

    end
    @show(AB_auxfields)
    #@show facnorm

    if U_AC > 0
        # AC_facnorm = (- 0.5*U_AC*(nA + nC))*dtau
        AC_gamma = acosh(exp(0.5*dtau*abs(U_AC)))

        AC_auxfields = zeros(2,2)
        for i = 1:2
            for j = 1:2
                AC_auxfields[i,j] = exp(AC_gamma*(-1)^(i+j))
            end
        end


    elseif U_AC <= 0

    #facnorm = (real(ET) - 0.5*U*(npar-1))*dtau - this will blow the weights exponentially
        # AC_facnorm = (real(ET))*dtau
        AC_gamma = acosh(exp(0.5*dtau*abs(U_AC)))

        AC_auxfields = zeros(2,2)
        for i = 1:2
            for j = 1:2
                AC_auxfields[i,j] = exp((-1)*dtau*U_AC*0.5)*exp(AC_gamma*(-1)^(1+j))
            end
        end
    end
    @show(AC_auxfields)
    # @show AC_facnorm

    if U_AD > 0
        # AD_facnorm = (- 0.5*U_AD*(nA + nD))*dtau
        AD_gamma = acosh(exp(0.5*dtau*abs(U_AD)))

        AD_auxfields = zeros(2,2)
        for i = 1:2
            for j = 1:2
                AD_auxfields[i,j] = exp(AD_gamma*(-1)^(i+j))
            end
        end


    elseif U_AD <= 0

    #facnorm = (real(ET) - 0.5*U*(npar-1))*dtau - this will blow the weights exponentially
        # AD_facnorm = (real(ET))*dtau
        AD_gamma = acosh(exp(0.5*dtau*abs(U_AD)))

        AD_auxfields = zeros(2,2)
        for i = 1:2
            for j = 1:2
                AD_auxfields[i,j] = exp((-1)*dtau*U_AD*0.5)*exp(AD_gamma*(-1)^(1+j))
            end
        end
    end
    @show(AD_auxfields)
    # @show AD_facnorm


    if U_BC > 0
        # BC_facnorm = (- 0.5*U_BC*(nB + nC))*dtau
        BC_gamma = acosh(exp(0.5*dtau*abs(U_BC)))

        BC_auxfields = zeros(2,2)
        for i = 1:2
            for j = 1:2
                BC_auxfields[i,j] = exp(BC_gamma*(-1)^(i+j))
            end
        end


    elseif U_BC <= 0

    #facnorm = (real(ET) - 0.5*U*(npar-1))*dtau - this will blow the weights exponentially
        # BC_facnorm = (real(ET))*dtau
        BC_gamma = acosh(exp(0.5*dtau*abs(U_BC)))

        BC_auxfields = zeros(2,2)
        for i = 1:2
            for j = 1:2
                BC_auxfields[i,j] = exp((-1)*dtau*U_BC*0.5)*exp(BC_gamma*(-1)^(1+j))
            end
        end
    end
    @show(BC_auxfields)
    # @show BC_facnorm

    if U_BD > 0
        # BD_facnorm = (- 0.5*U_BD*(nB + nD))*dtau
        BD_gamma = acosh(exp(0.5*dtau*abs(U_BD)))

        BD_auxfields = zeros(2,2)
        for i = 1:2
            for j = 1:2
                BD_auxfields[i,j] = exp(BD_gamma*(-1)^(i+j))
            end
        end


    elseif U_BD <= 0

    #facnorm = (real(ET) - 0.5*U*(npar-1))*dtau - this will blow the weights exponentially
        # BD_facnorm = (real(ET))*dtau
        BD_gamma = acosh(exp(0.5*dtau*abs(U_BD)))

        BD_auxfields = zeros(2,2)
        for i = 1:2
            for j = 1:2
                BD_auxfields[i,j] = exp((-1)*dtau*U_BD*0.5)*exp(BD_gamma*(-1)^(1+j))
            end
        end
    end
    @show(BD_auxfields)
    # @show BD_facnorm

    if U_CD > 0
        # CD_facnorm = (- 0.5*U_CD*(nC + nD))*dtau
        CD_gamma = acosh(exp(0.5*dtau*abs(U_CD)))

        CD_auxfields = zeros(2,2)
        for i = 1:2
            for j = 1:2
                CD_auxfields[i,j] = exp(CD_gamma*(-1)^(i+j))
            end
        end


    elseif U_CD <= 0

    #facnorm = (real(ET) - 0.5*U*(npar-1))*dtau - this will blow the weights exponentially
        # CD_facnorm = (real(ET))*dtau
        CD_gamma = acosh(exp(0.5*dtau*abs(U_CD)))

        CD_auxfields = zeros(2,2)
        for i = 1:2
            for j = 1:2
                CD_auxfields[i,j] = exp((-1)*dtau*U_CD*0.5)*exp(CD_gamma*(-1)^(1+j))
            end
        end
    end
    @show(CD_auxfields)
    # @show CD_facnorm


    return facnorm, AB_auxfields,AC_auxfields,AD_auxfields,BC_auxfields,BD_auxfields,CD_auxfields

end


# use spaces when calling for functions
nsysparams = nParamsSys()

runparams = ParamsRun()

nHK_up,nHK_do, nProjK_up, nProjK_do = n_ini_internal(nsysparams, runparams)
pHK_up,pHK_do, pProjK_up, pProjK_do = p_ini_internal(nsysparams, runparams)

PhiT, ET = ini_Phi_E(nsysparams, nHK_up,nHK_do,pHK_up,pHK_do)

Phis, ws, Os = ini_pop(nsysparams, runparams, PhiT)


facnorm, AB_auxfields,AC_auxfields,AD_auxfields,BC_auxfields,BD_auxfields,CD_auxfields = ini_aux(nsysparams, runparams, ET)


end
