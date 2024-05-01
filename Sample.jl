module Sample
using LinearAlgebra

struct nParamsSys
    nx::Int64
    ny::Int64
    nz::Int64
    alpha_n::Float64
    nA::Int64
    nB::Int64
    nC::Int64
    nD::Int64

    mn::Float64


    ntwistx::Float64
    ntwisty::Float64
    ntwistz::Float64

    U_AB::Float64
    U_AC::Float64
    U_AD::Float64
    U_BD::Float64
    U_BC::Float64
    U_CD::Float64

    t::Float64

#Important!!! Assigning order matters
#Introducing neutrons and protons seperately. For the moment keep the same lattice for both protons and neutrons.
#Simple case of deuteron. My argument is deuteron is not a closed shell as we have neutron =1 and proton =1 . So we will need to account that for the calculation.
#Maybe we will need twists.
    #function nParamsSys(nx=3, ny=3, nz=3, alpha_n =1.0, nA = 1, nB = 1,nC = 1, nD = 1, mn=1.0, ntwistx= 0., ntwisty = 0., ntwistz=0.,U_AB=-11.839,U_AC= -11.839, U_AD= -11.839,U_BD= -11.839, U_BC= -11.839,U_CD=-11.839, t=1)
    #function nParamsSys(nx=3, ny=3, nz=3, alpha_n =1.0, nA = 1, nB = 1,nC = 1, nD = 1, mn=1.0, ntwistx= 0., ntwisty = 0., ntwistz=0.,U_AB=-5.9195,U_AC= -5.9195, U_AD= -5.9195,U_BD= -5.9195, U_BC= -5.9195,U_CD=-5.9195, t=1)
    function nParamsSys(nx=3, ny=3, nz=3, alpha_n =1.0, nA = 1, nB = 1,nC = 1, nD = 1, mn=1.0, ntwistx= 0., ntwisty = 0., ntwistz=0.,U_AB=-13.24,U_AC= -11.839, U_AD= -13.47,U_BD= -11.839, U_BC= -13.47,U_CD=-13.551, t=1)
    # function nParamsSys(nx=3, ny=3, nz=3, alpha_n =1.0, nA =1, nB = 1,nC = 1, nD = 1, mn=1.0, ntwistx= 0., ntwisty = 0., ntwistz=0.,U_AB=-10.695,U_AC= -8.571, U_AD= -10.886,U_BD= -8.571, U_BC= -10.886,U_CD=-11.024, t=1)
    #@show(nx, ny, nz,alpha_n, nA , nB ,nC , nD , mn, ntwistx, ntwisty, ntwistz,U_AB, U_AC, U_AD,U_BD, U_BC, U_CD, t)

    new(nx, ny, nz,alpha_n, nA, nB ,nC , nD , mn, ntwistx, ntwisty, ntwistz,U_AB, U_AC, U_AD,U_BD, U_BC, U_CD, t)
    end

end



struct ParamsRun
    dtau::Float64

    nwalk::Int64
    nstepsperblock::Int64
    nequilblocks::Int64
    nmeasblocks::Int64

    itv_ortho::Int64
    itv_popu::Int64
    itv_meas::Int64
    # Always keep ortho less than the number of steps per block. Keep itv_meas equals to nmeasblocks. Keep  nsteps per block/itv_popu =0.
    function ParamsRun(dtau=0.01, nwalk=50, nstepsperblock=10, nequilblocks=0, nmeasblocks=10, itv_ortho= 5, itv_popu=nstepsperblock, itv_meas=nmeasblocks)
    # @show(dtau, nwalk, nstepsperblock, nequilblocks, nmeasblocks, itv_ortho, itv_popu, itv_meas)

    new(dtau, nwalk, nstepsperblock, nequilblocks, nmeasblocks, itv_ortho, itv_popu, itv_meas)
    end

end

# nsysparams = nParamsSys()
# runparams = ParamsRun()

end
