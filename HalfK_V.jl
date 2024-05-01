module HalfK_V

include("H_K.jl")
include("Initialization.jl")


using .Main: nParamsSys
using .H_K: nh_k,ph_k
using .Initialization: n_ini_internal,p_ini_internal,ini_Phi_E, ini_pop,ini_aux
using LinearAlgebra
using Random

function halfK!(nsysparams::nParamsSys, onephi, onew, oneO, nProjK_up,nProjK_do,pProjK_up,pProjK_do, PhiT)
    nup = nsysparams.nA
    ndo = nsysparams.nB
    pup = nsysparams.nC
    pdo = nsysparams.nD
    npar = nup + ndo
    apar = nup + ndo + pup + pdo


    onephi[:,1:nup] = nProjK_up*onephi[:,1:nup]
    onephi[:,nup + 1:npar] = nProjK_do*onephi[:,nup + 1:npar]

    onephi[:,npar + 1:npar + pup] = pProjK_up*onephi[:,npar + 1:npar + pup]
    onephi[:,npar + pup + 1:apar] = pProjK_do*onephi[:,npar + pup + 1:apar]

    n_invO_A = inv((PhiT[:,1:nup])' * onephi[:,1:nup])
    n_invO_B = inv((PhiT[:,nup + 1:npar])' * onephi[:,nup + 1:npar])
    #@show n_invO_A
    p_invO_C = inv((PhiT[:,(npar + 1):npar + pup])' * onephi[:,(npar + 1):npar + pup])
    p_invO_D = inv((PhiT[:,npar + pup + 1:apar])' * onephi[:,npar + pup + 1:apar])

    newO = 1/((det(n_invO_A)*det(n_invO_B)*det(p_invO_C)*det(p_invO_D)))

    ratioO = newO/oneO
    #@show ratioO
    if real(ratioO) > 0
        oneO = newO
        onew *= real(ratioO)
        #println("onew at HalfK: ",onew)
    else
        println("kinetic death")
        onew = 0
    end


    return onephi, onew, oneO, n_invO_A, n_invO_B,p_invO_C, p_invO_D
end


# Note that the phi and PhiT in this function are for a single lattice site,
# so I  use different variables for them, in order to avoid confusion

function V!(nsysparams::nParamsSys, onlyphi, onew, oneO, n_invO_A, n_invO_B,p_invO_C, p_invO_D, onlyPhiT, AB_auxfields,AC_auxfields, AD_auxfields,BC_auxfields,BD_auxfields,CD_auxfields,rng)
    nup = nsysparams.nA
    ndo = nsysparams.nB
    pup = nsysparams.nC
    pdo = nsysparams.nD
    npar = nup + ndo
    apar = nup + ndo + pup + pdo

    species = zeros(Int8, 1, 4)
    species[1] = nup
    species[2] = ndo
    species[3] = pup
    species[4] = pdo
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



 # AB overlap

    AB_Gii = zeros(Complex{Float64},2,1)
    AB_RR = ones(Complex{Float64},2,2)
    AB_matone = copy(AB_RR)

    AB_A1 = transpose(onlyphi[1:nup]) * n_invO_A
    AB_A2 = transpose(onlyphi[nup + 1:npar]) * n_invO_B
    AB_B1 = n_invO_A * conj(onlyPhiT[1:nup])
    AB_B2 = n_invO_B * conj(onlyPhiT[nup + 1:npar])
    AB_Gii[1] = AB_A1 * conj(onlyPhiT[1:nup])
    AB_Gii[2] = AB_A2 * conj(onlyPhiT[nup + 1:npar])
    AB_RR = (AB_auxfields - AB_matone) .* hcat(AB_Gii, AB_Gii) + AB_matone

    AB_ratioO_now = AB_RR[1,:] .* AB_RR[2,:]
    # probability for sampling the auxfields at each lattice point
    if U_AB >= 0.
        ratioO_re1 = max(real(AB_ratioO_now[1]), 0.)
        ratioO_re2 = max(real(AB_ratioO_now[2]), 0.)
    else
        AB_ratioO_now = AB_RR[1,:] .* AB_RR[2,:] ./AB_auxfields[1,:]
        # below two lines also should work, and don't update oneO with auxfiled selection if you use below two lines.
        #ratioO_re1 = max(real(AB_ratioO_now[1] ./AB_auxfields[1,1]), 0.)
        #ratioO_re2 = max(real(AB_ratioO_now[2] ./AB_auxfields[2,2]), 0.)
        ratioO_re1 = max(real(AB_ratioO_now[1]) , 0.)
        ratioO_re2 = max(real(AB_ratioO_now[2]) , 0.)

    end
    #normalization factor
    ratioO_sum = ratioO_re1 + ratioO_re2

    if ratioO_sum <= 0.
        onew = 0.
    end


    if real(onew)> 0.

        onew *= 0.5*ratioO_sum
        #println("onew updated after AB overlap: ", onew)
        if ratioO_re1/ratioO_sum >= rand(rng)
            AB_xspin = 1
        else
            AB_xspin = 2
        end
    end

    onlyphi[1:nup] *= (AB_auxfields[1,AB_xspin])
    onlyphi[nup+1:npar] *= (AB_auxfields[2,AB_xspin] )

    oneO *= AB_ratioO_now[AB_xspin]
    #This loop below is linked to the auxfield dependent AB_ratioO_now.
    if U_AB < 0
        oneO *= AB_auxfields[1,AB_xspin]
    end

    n_invO_A += (((1 - AB_auxfields[1,AB_xspin])/AB_RR[1,AB_xspin])* AB_B1 * AB_A1)
    n_invO_B += (((1 - AB_auxfields[2,AB_xspin])/AB_RR[2,AB_xspin]) * AB_B2 * AB_A2)

    #AC overlap

    AC_Gii = zeros(Complex{Float64},2,1)
    AC_RR = ones(Complex{Float64},2,2)
    AC_matone = copy(AC_RR)

    AC_A1 = transpose(onlyphi[1:nup]) * n_invO_A
    AC_A2 = transpose(onlyphi[npar + 1:npar + pup]) * p_invO_C
    AC_B1 = n_invO_A * conj(onlyPhiT[1:nup])
    AC_B2 = p_invO_C * conj(onlyPhiT[npar + 1:npar + pup])
    AC_Gii[1] = AC_A1 * conj(onlyPhiT[1:nup])
    AC_Gii[2] = AC_A2 * conj(onlyPhiT[npar + 1:npar + pup])
    AC_RR = (AC_auxfields - AC_matone) .* hcat(AC_Gii, AC_Gii) + AC_matone

    #@show AC_RR

    AC_ratioO_now = AC_RR[1,:] .* AC_RR[2,:]

    if U_AC >= 0.
        ratioO_re1 = max(real(AC_ratioO_now[1]), 0.)
        ratioO_re2 = max(real(AC_ratioO_now[2]), 0.)
    else
        AC_ratioO_now = AC_RR[1,:] .* AC_RR[2,:] ./ AC_auxfields[1,:]

        #ratioO_re1 = max(real(AC_ratioO_now[1] ./AC_auxfields[1,1]), 0.)
        #ratioO_re2 = max(real(AC_ratioO_now[2] ./AC_auxfields[2,2]), 0.)
        ratioO_re1 = max(real(AC_ratioO_now[1]) , 0.)
        ratioO_re2 = max(real(AC_ratioO_now[2]) , 0.)
    end

    ratioO_sum = (ratioO_re1 + ratioO_re2)

    if ratioO_sum <= 0.
        onew = 0.
    end


    if real(onew)> 0.

        onew *= 0.5*ratioO_sum
        #println("onew updated after AC overlap: ", onew)

        if ratioO_re1/ratioO_sum >= rand(rng)
            AC_xspin = 1
        else
            AC_xspin = 2
        end
    end

    onlyphi[1:nup] *= (AC_auxfields[1,AC_xspin])
    onlyphi[npar+1:npar+pup] *= (AC_auxfields[2,AC_xspin])

    oneO *= AC_ratioO_now[AC_xspin]
    if U_AC < 0
        oneO *= AC_auxfields[1,AC_xspin]
    end

    n_invO_A += (((1 - AC_auxfields[1,AC_xspin])/AC_RR[1,AC_xspin])* AC_B1 * AC_A1)
    p_invO_C += (((1 - AC_auxfields[2,AC_xspin])/AC_RR[2,AC_xspin]) * AC_B2 * AC_A2)


    # AD overlap

    AD_Gii = zeros(Complex{Float64},2,1)
    AD_RR = ones(Complex{Float64},2,2)
    AD_matone = copy(AD_RR)

    AD_A1 = transpose(onlyphi[1:nup]) * n_invO_A
    AD_A2 = transpose(onlyphi[npar + pup + 1:apar]) * p_invO_D
    AD_B1 = n_invO_A * conj(onlyPhiT[1:nup])
    AD_B2 = p_invO_D * conj(onlyPhiT[npar + pup + 1:apar])
    AD_Gii[1] = AD_A1 * conj(onlyPhiT[1:nup])
    AD_Gii[2] = AD_A2 * conj(onlyPhiT[npar + pup + 1:apar])
    AD_RR = (AD_auxfields - AD_matone) .* hcat(AD_Gii, AD_Gii) + AD_matone

    #@show AD_RR

    AD_ratioO_now = AD_RR[1,:] .* AD_RR[2,:]

    if U_AD >= 0.
        ratioO_re1 = max(real(AD_ratioO_now[1]), 0.)
        ratioO_re2 = max(real(AD_ratioO_now[2]), 0.)
    else
        AD_ratioO_now = AD_RR[1,:] .* AD_RR[2,:] ./ AD_auxfields[1,:]

        ratioO_re1 = max(real(AD_ratioO_now[1]) , 0.)
        ratioO_re2 = max(real(AD_ratioO_now[2]) , 0.)
        #ratioO_re1 = max(real(AD_ratioO_now[1] ./AD_auxfields[1,1]), 0.)
        #ratioO_re2 = max(real(AD_ratioO_now[2] ./AD_auxfields[2,2]), 0.)
    end

    ratioO_sum = (ratioO_re1 + ratioO_re2)

    if ratioO_sum <= 0.
        onew = 0.
    end


    if real(onew)> 0.

        onew *= 0.5*ratioO_sum
        #println("onew updated after AD overlap: ", onew)

        if ratioO_re1/ratioO_sum >= rand(rng)
            AD_xspin = 1
        else
            AD_xspin = 2
        end
    end

    onlyphi[1:nup] *= (AD_auxfields[1,AD_xspin])
    onlyphi[npar+nup+1:apar] *= (AD_auxfields[2,AD_xspin])

    oneO *= AD_ratioO_now[AD_xspin]
    if U_AD < 0
        oneO *= AD_auxfields[1,AD_xspin]
    end

    n_invO_A += (((1 - AD_auxfields[1,AD_xspin])/AD_RR[1,AD_xspin])* AD_B1 * AD_A1)
    p_invO_D += (((1 - AD_auxfields[2,AD_xspin])/AD_RR[2,AD_xspin])*AD_B2 * AD_A2)

    # BC overlap

    BC_Gii = zeros(Complex{Float64},2,1)
    BC_RR = ones(Complex{Float64},2,2)
    BC_matone = copy(BC_RR)

    BC_A1 = transpose(onlyphi[nup + 1:npar]) * n_invO_B
    BC_A2 = transpose(onlyphi[npar + 1:npar + pup]) * p_invO_C
    BC_B1 = n_invO_B * conj(onlyPhiT[nup + 1:npar])
    BC_B2 = p_invO_C * conj(onlyPhiT[npar + 1:npar + pup])
    BC_Gii[1] = BC_A1 * conj(onlyPhiT[nup + 1:npar])
    BC_Gii[2] = BC_A2 * conj(onlyPhiT[npar + 1:npar + pup])
    BC_RR = (BC_auxfields - BC_matone) .* hcat(BC_Gii, BC_Gii) + BC_matone

    #@show BC_RR


    BC_ratioO_now = BC_RR[1,:] .* BC_RR[2,:]

    if U_BC >= 0.
        ratioO_re1 = max(real(BC_ratioO_now[1]), 0.)
        ratioO_re2 = max(real(BC_ratioO_now[2]), 0.)
    else
        BC_ratioO_now = BC_RR[1,:] .* BC_RR[2,:] ./ BC_auxfields[1,:]

        #ratioO_re1 = max(real(BC_ratioO_now[1] ./BC_auxfields[1,1]), 0.)
        #ratioO_re2 = max(real(BC_ratioO_now[2] ./BC_auxfields[2,2]), 0.)
        ratioO_re1 = max(real(BC_ratioO_now[1]) , 0.)
        ratioO_re2 = max(real(BC_ratioO_now[2]) , 0.)
    end

    ratioO_sum = (ratioO_re1 + ratioO_re2)

    if ratioO_sum <= 0.
        onew = 0.
    end


    if real(onew)> 0.

        onew *= 0.5*ratioO_sum
        #println("onew updated after BC overlap: ", onew)

        if ratioO_re1/ratioO_sum >= rand(rng)
            BC_xspin = 1
        else
            BC_xspin = 2
        end
    end
    onlyphi[nup+1:npar] *= (BC_auxfields[1,BC_xspin] )
    onlyphi[npar+1:npar+pup] *= (BC_auxfields[2,BC_xspin])

    oneO *= BC_ratioO_now[BC_xspin]
    if U_BC < 0

        oneO *= BC_auxfields[1,BC_xspin]
    end

    n_invO_B += (((1 - BC_auxfields[1,BC_xspin])/BC_RR[1,BC_xspin]) * BC_B1 * BC_A1)
    p_invO_C += (((1 - BC_auxfields[2,BC_xspin])/BC_RR[2,BC_xspin]) * BC_B2 * BC_A2)

    # BD overlap

    BD_Gii = zeros(Complex{Float64},2,1)
    BD_RR = ones(Complex{Float64},2,2)
    BD_matone = copy(BD_RR)

    BD_A1 = transpose(onlyphi[nup + 1:npar]) * n_invO_B
    BD_A2 = transpose(onlyphi[npar + pup + 1:apar]) * p_invO_D
    BD_B1 = n_invO_B * conj(onlyPhiT[nup + 1:npar])
    BD_B2 = p_invO_D * conj(onlyPhiT[npar + pup + 1:apar])
    BD_Gii[1] = BD_A1 * conj(onlyPhiT[nup + 1:npar])
    BD_Gii[2] = BD_A2 * conj(onlyPhiT[npar + pup + 1:apar])
    BD_RR = (BD_auxfields - BD_matone) .* hcat(BD_Gii, BD_Gii) + BD_matone

    #@show BD_RR


    BD_ratioO_now = BD_RR[1,:] .* BD_RR[2,:]

    if U_BD >= 0.
        ratioO_re1 = max(real(BD_ratioO_now[1]), 0.)
        ratioO_re2 = max(real(BD_ratioO_now[2]), 0.)
    else
        BD_ratioO_now = BD_RR[1,:] .* BD_RR[2,:] ./ BD_auxfields[1,:]

        #ratioO_re1 = max(real(BD_ratioO_now[1] ./BD_auxfields[1,1]), 0.)
        #ratioO_re2 = max(real(BD_ratioO_now[2] ./BD_auxfields[2,2]), 0.)
        ratioO_re1 = max(real(BD_ratioO_now[1]) , 0.)
        ratioO_re2 = max(real(BD_ratioO_now[2]) , 0.)
    end

    ratioO_sum = (ratioO_re1 + ratioO_re2)

    if ratioO_sum <= 0.
        onew = 0.
    end


    if real(onew)> 0.

        onew *= 0.5*ratioO_sum
        #println("onew updated after BD overlap: ", onew)

        if ratioO_re1/ratioO_sum >= rand(rng)
            BD_xspin = 1
        else
            BD_xspin = 2
        end
    end

    onlyphi[nup+1:npar] *= (BD_auxfields[1,BD_xspin] )
    onlyphi[npar+nup+1:apar] *= (BD_auxfields[2,BD_xspin])

    oneO *= BD_ratioO_now[BD_xspin]
    if U_BD < 0
        oneO *= BD_auxfields[1,BD_xspin]
    end

    n_invO_B += (((1 - BD_auxfields[1,BD_xspin])/BD_RR[1,BD_xspin]) * BD_B1 * BD_A1)
    p_invO_D += (((1 - BD_auxfields[2,BD_xspin])/BD_RR[2,BD_xspin])*BD_B2 * BD_A2)

    # CD overlap

    CD_Gii = zeros(Complex{Float64},2,1)
    CD_RR = ones(Complex{Float64},2,2)
    CD_matone = copy(CD_RR)

    CD_A1 = transpose(onlyphi[(npar + 1): (npar + pup)]) * p_invO_C
    CD_A2 = transpose(onlyphi[(npar + pup + 1):apar]) * p_invO_D
    CD_B1 = p_invO_C * conj(onlyPhiT[(npar + 1): (npar+pup)])
    CD_B2 = p_invO_D * conj(onlyPhiT[(npar + pup +1):apar])
    CD_Gii[1] = CD_A1 * conj(onlyPhiT[(npar+1): (npar + pup)])
    CD_Gii[2] = CD_A2 * conj(onlyPhiT[(npar + pup + 1):apar])
    CD_RR = (CD_auxfields - CD_matone) .* hcat(CD_Gii, CD_Gii) + CD_matone

    #@show CD_RR


    CD_ratioO_now = CD_RR[1,:] .* CD_RR[2,:]

    if U_CD >= 0.
        ratioO_re1 = max(real(CD_ratioO_now[1]), 0.)
        ratioO_re2 = max(real(CD_ratioO_now[2]), 0.)
    else
        CD_ratioO_now = CD_RR[1,:] .* CD_RR[2,:] ./ CD_auxfields[1,:]

        #ratioO_re1 = max(real(CD_ratioO_now[1] ./CD_auxfields[1,1]), 0.)
        #ratioO_re2 = max(real(CD_ratioO_now[2] ./CD_auxfields[2,2]), 0.)
        ratioO_re1 = max(real(CD_ratioO_now[1]) , 0.)
        ratioO_re2 = max(real(CD_ratioO_now[2]) , 0.)
    end


    ratioO_sum = (ratioO_re1 + ratioO_re2)

    if ratioO_sum <= 0.
        onew = 0.
    end


    if real(onew)> 0.

        onew *= 0.5*ratioO_sum
        #println("onew updated after CD overlap: ", onew)

        if ratioO_re1/ratioO_sum >= rand(rng)
            CD_xspin = 1
        else
            CD_xspin = 2
        end
    end

    onlyphi[npar+1:npar+pup] *= (CD_auxfields[1,CD_xspin])
    onlyphi[npar+nup+1:apar] *= (CD_auxfields[2,CD_xspin])

    oneO *= CD_ratioO_now[CD_xspin]
    if U_CD < 0
        oneO *= CD_auxfields[1,CD_xspin]
    end

    p_invO_C += (((1 - CD_auxfields[1,CD_xspin])/CD_RR[1,CD_xspin]) * CD_B1 * CD_A1)
    p_invO_D += (((1 - CD_auxfields[2,CD_xspin])/CD_RR[2,CD_xspin])*CD_B2 * CD_A2)



    return onlyphi, onew, oneO, n_invO_A, n_invO_B,p_invO_C, p_invO_D
end


#sysparams = ParamsSys()
#runparams = ParamsRun()
#HK, ProjK = ini_internal(sysparams, runparams)
#PhiT, ET = ini_Phi_E(sysparams, HK)
#Phis, ws, Os = ini_pop(sysparams, runparams, PhiT)
#facnorm, auxfields = ini_aux(sysparams, runparams, ET)
#
#onephi = reshape(collect(1:24.) .+ collect(24:-1:1)*im, (8,3))
#
#println("BEFORE")
#display(onephi)
#display(0.7)
#display(0.8)
#println(" ")

#onephi, onew, oneO, invO_up, invO_do = halfK!(sysparams, onephi, 0.7, 0.8, ProjK, PhiT)
#println("AFTER")
#display(onephi)
#println(" ")
#display(onew)
#println(" ")
#display(oneO)
#println(" ")
#display(invO_up)
#println(" ")
#display(invO_do)
#println(" ")


#onephi = reshape(collect(11:34.) .- collect(34:-1:11)*im, (8,3))
#onlyphi = onephi[2,:]
#onlyPhiT = PhiT[2,:]
#println("BEFORE")
#display(onlyphi)
#display(1.7)
#display(1.8)
#println(" ")

#onlyphi, onew, oneO, invO_up, invO_do = V!(1, sysparams, onlyphi, 1.7, 1.8, invO_up, invO_do, onlyPhiT, auxfields)
#println("AFTER")
#display(onlyphi)
#println(" ")
#display(onew)
#println(" ")
#display(oneO)
#println(" ")
#display(invO_up)
#println(" ")
#display(invO_do)
#println(" ")
#onlyphi, onew, oneO, invO_up, invO_do = V!(2, sysparams, onlyphi, onew, oneO, invO_up, invO_do, onlyPhiT, auxfields)
#println("AFTER")
#display(onlyphi)
#println(" ")
#display(onew)
#println(" ")
#display(oneO)
#println(" ")
#display(invO_up)
#println(" ")
#display(invO_do)
#println(" ")

end
