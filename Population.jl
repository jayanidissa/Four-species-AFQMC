module Population

using .Main: nParamsSys, ParamsRun


using LinearAlgebra
using Random

function population!(nsysparams::nParamsSys, runparams::ParamsRun, Phis, ws, Os, rng)
    nA = nsysparams.nA
    nB = nsysparams.nB
    nC = nsysparams.nC
    nD = nsysparams.nD
    apar = nA + nB + nC + nD
    nsites = nsysparams.nx*nsysparams.ny*nsysparams.nz
    walk = runparams.nwalk

    newPhis = zeros(Complex{Float64}, nsites, apar, walk)
    newOs = zeros(Complex{Float64}, walk,1)

    sca = walk/sum(real(ws))
    #@show sca
    wsum = -rand(rng)
    mint = 0

    for iwa = 1:walk
        wsum += sca*ws[iwa]
        #@show wsum
        n = ceil(Integer, real(wsum))


        for j = (mint+1):n
            newPhis[:, :, j] = Phis[:, :, iwa]
            newOs[j] = Os[iwa]
        end
        mint = n
    end

    Phis = newPhis
    ws = ones(Complex{Float64}, walk,1)
    Os = newOs
    return Phis, ws, Os
end

end
