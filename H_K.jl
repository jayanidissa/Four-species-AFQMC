module H_K

using .Main: nParamsSys
using LinearAlgebra
using Statistics

function Adisp(nsysparams::nParamsSys,kn)

    alpha_n = nsysparams.alpha_n
    h_bar= 1
    mn= nsysparams.mn
    nA = nsysparams.nA
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
        A_AB = 11.279
        beta_AB = 0.917
    else
        A_AB = 0
        beta_AB = 0
    end
    if mat_species[1,3] > 0
        A_AC = 5.567
        beta_AC = 0.5861
    else
        A_AC = 0
        beta_AC = 0
    end
    if mat_species[1,4] > 0
        A_AD = 11.465
        beta_AD = 0.924
    else
        A_AD = 0
        beta_AD = 0
    end

    A_hk_all = (A_AB + A_AC + A_AD)
    B_hk_all = -4.00


    beta_sq_A = (beta_AB^2 + beta_AC^2 + beta_AD^2)

    #beta_unitary = 0.16137
    beta_unitary = 0.5861

    #neutron dispersions

    #magic/ ek(n2)
    # ekh_A = (3 -cos(kn[1]*alpha_n) -cos(kn[2]*alpha_n) -cos(kn[3]*alpha_n))*h_bar^2/(mn*alpha_n^2)
    # #for non Unitary for 3D replace diagonal term 3 with 0 to map with the original CPMC and multiply  t_jm by 2. For 2D the diagonal term is 1, and for 1D the diagonal is 2.
    # X_A = ekh_A*mn*alpha_n^2/6
    # ek_A = ekh_A*(1 + A_hk_all*X_A + B_hk_all*X_A^2)
    #ekh
    # ek_A = (3 -cos(kn[1]*alpha_n) -cos(kn[2]*alpha_n) -cos(kn[3]*alpha_n))*h_bar^2/(mn*alpha_n^2)
    #ek(2)
    # ek_A = 0.5*h_bar^2(kn[1]^2 + kn[2]^2 + kn[3]^2)/mn
    # @show ek_A

    #ek(4) for unitary fermions : ( r_e = 0)
    ek_A = 0.5*h_bar^2*(kn[1]^2 + kn[2]^2 + kn[3]^2)*(1 + beta_unitary^2*(kn[1]^2 + kn[2]^2 + kn[3]^2)*alpha_n^2)/mn

    #ek(4)/ek(n1) for finite effective range
    #ek_A = 0.5*h_bar^2*(kn[1]^2 + kn[2]^2 + kn[3]^2)*(1 + beta_sq_A*(kn[1]^2 + kn[2]^2 + kn[3]^2)*alpha_n^2)/mn


    return ek_A
end

function Bdisp(nsysparams::nParamsSys,kn)

    alpha_n = nsysparams.alpha_n
    h_bar= 1
    mn= nsysparams.mn

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
        A_AB = 11.279
        beta_AB = 0.917
    else
        A_AB = 0
        beta_AB = 0
    end

    if mat_species[2,3] > 0
        A_BC = 11.465
        beta_BC = 0.924
    else
        A_BC = 0
        beta_BC = 0
    end
    if mat_species[2,4] > 0
        A_BD = 5.567
        beta_BD = 0.5861
    else
        A_BD = 0
        beta_BD = 0
    end


    A_hk_all = (A_AB + A_BD + A_BC)
    B_hk_all = -4.00


    beta_sq_B = (beta_AB^2 + beta_BC^2 + beta_BD^2)

    #beta_unitary = 0.16137
    beta_unitary = 0.5861

    #neutron dispersions

    #magic/ ek(n2)
        # ekh_B = (3 -cos(kn[1]*alpha_n) -cos(kn[2]*alpha_n) -cos(kn[3]*alpha_n))*h_bar^2/(mn*alpha_n^2)
        # #for non Unitary for 3D replace diagonal term 3 with 0 to map with the original CPMC and multiply  t_jm by 2. For 2D the diagonal term is 1, and for 1D the diagonal is 2.
        # X_B = ekh_B*mn*alpha_n^2/6
        # ek_B = ekh_B*(1 + A_hk_all*X_B + B_hk_all*X_B^2)
    #ekh
        # ek_B = (3 -cos(kn[1]*alpha_n) -cos(kn[2]*alpha_n) -cos(kn[3]*alpha_n))*h_bar^2/(mn*alpha_n^2)
    #ek(2)
        # ek_B = 0.5*h_bar^2(kn[1]^2 + kn[2]^2 + kn[3]^2)/mn
        # @show ek_B

    #ek(4) for unitary fermions : ( r_e = 0)
    ek_B = 0.5*h_bar^2*(kn[1]^2 + kn[2]^2 + kn[3]^2)*(1 + beta_unitary^2*(kn[1]^2 + kn[2]^2 + kn[3]^2)*alpha_n^2)/mn

    #ek(4)/ek(n1) for finite effective range
    #ek_B = 0.5*h_bar^2*(kn[1]^2 + kn[2]^2 + kn[3]^2)*(1 + beta_sq_B*(kn[1]^2 + kn[2]^2 + kn[3]^2)*alpha_n^2)/mn


    return ek_B
end

function Cdisp(nsysparams::nParamsSys,kn)

    alpha_n = nsysparams.alpha_n
    h_bar= 1
    mn= nsysparams.mn
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


    if mat_species[1,3] > 0
        A_AC = 5.567
        beta_AC = 0.5861
    else
        A_AC = 0
        beta_AC = 0
    end

    if mat_species[2,3] > 0
        A_BC = 11.465
        beta_BC = 0.924
    else
        A_BC = 0
        beta_BC = 0
    end

    if mat_species[3,4] > 0
        A_CD = 12.061
        beta_CD = 0.951
    else
        A_CD = 0
        beta_CD = 0
    end


    A_hk_all = (A_AC + A_BC + A_CD)
    B_hk_all = -4.00

    beta_sq_C = (beta_AC^2 + beta_BC^2 + beta_CD^2)

    #beta_unitary = 0.16137
    beta_unitary = 0.5861

    #proton dispersions

    #magic/ ek(n2)
        # ekh_C = (3 -cos(kn[1]*alpha_n) -cos(kn[2]*alpha_n) -cos(kn[3]*alpha_n))*h_bar^2/(mn*alpha_n^2)
        # #for non Unitary for 3D replace diagonal term 3 with 0 to map with the original CPMC and multiply  t_jm by 2. For 2D the diagonal term is 1, and for 1D the diagonal is 2.
        # X_C = ekh_C*mn*alpha_n^2/6
        # ek_C = ekh_C*(1 + A_hk_all*X_C + B_hk_all*X_C^2)
    #ekh
        # ek_C = (3 -cos(kn[1]*alpha_n) -cos(kn[2]*alpha_n) -cos(kn[3]*alpha_n))*h_bar^2/(mn*alpha_n^2)
    #ek(2)
        # ek_C = 0.5*h_bar^2(kn[1]^2 + kn[2]^2 + kn[3]^2)/mn
        # @show ek_C

    #ek(4) for unitary fermions : ( r_e = 0)
    ek_C = 0.5*h_bar^2*(kn[1]^2 + kn[2]^2 + kn[3]^2)*(1 + beta_unitary^2*(kn[1]^2 + kn[2]^2 + kn[3]^2)*alpha_n^2)/mn

    #ek(4)/ek(n1) for finite effective range
    #ek_C = 0.5*h_bar^2*(kn[1]^2 + kn[2]^2 + kn[3]^2)*(1 + beta_sq_C*(kn[1]^2 + kn[2]^2 + kn[3]^2)*alpha_n^2)/mn


    return ek_C
end

function Ddisp(nsysparams::nParamsSys,kn)

    alpha_n = nsysparams.alpha_n
    h_bar= 1
    mn= nsysparams.mn

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


    if mat_species[1,4] > 0
        A_AD = 11.465
        beta_AD = 0.924
    else
        A_AD = 0
        beta_AD = 0
    end
    if mat_species[2,4] > 0
        A_BD = 5.567
        beta_BD = 0.5861
    else
        A_BD = 0
        beta_BD = 0
    end
    if mat_species[3,4] > 0
        A_CD = 12.061
        beta_CD = 0.951
    else
        A_CD = 0
        beta_CD = 0
    end


    A_hk_all = (A_AD + A_BD + A_CD)
    B_hk_all = -4.00

    beta_sq_D = (beta_AD^2 + beta_BD^2 + beta_CD^2)

    #beta_unitary = 0.16137
    beta_unitary = 0.5861

    #proton dispersions

    #magic/ ek(n2)
        # ekh_D = (3 -cos(kn[1]*alpha_n) -cos(kn[2]*alpha_n) -cos(kn[3]*alpha_n))*h_bar^2/(mn*alpha_n^2)
        # #for non Unitary for 3D replace diagonal term 3 with 0 to map with the original CPMC and multiply  t_jm by 2. For 2D the diagonal term is 1, and for 1D the diagonal is 2.
        # X_D = ekh_D*mn*alpha_n^2/6
        # ek_D = ekh_D*(1 + A_hk_all*X_D + B_hk_all*X_D^2)
    #ekh
        # ek_D = (3 -cos(kn[1]*alpha_n) -cos(kn[2]*alpha_n) -cos(kn[3]*alpha_n))*h_bar^2/(mn*alpha_n^2)
    #ek(2)
        # ek_D = 0.5*h_bar^2(kn[1]^2 + kn[2]^2 + kn[3]^2)/mn
        # @show ek_D

    #ek(4) for unitary fermions : ( r_e = 0)
    ek_D = 0.5*h_bar^2*(kn[1]^2 + kn[2]^2 + kn[3]^2)*(1 + beta_unitary^2*(kn[1]^2 + kn[2]^2 + kn[3]^2)*alpha_n^2)/mn

    #ek(4)/ek(n1) for finite effective range
    #ek_D = 0.5*h_bar^2*(kn[1]^2 + kn[2]^2 + kn[3]^2)*(1 + beta_sq_D*(kn[1]^2 + kn[2]^2 + kn[3]^2)*alpha_n^2)/mn
    # @show ek_D

    return ek_D
end



function nh_k(nsysparams::nParamsSys)
    #neutrons
    nA = nsysparams.nA
    nB = nsysparams.nB
    nx = nsysparams.nx
    ny = nsysparams.ny
    nz = nsysparams.nz

    alpha_n = nsysparams.alpha_n

    ntx =nsysparams.ntwistx*pi
    nty = nsysparams.ntwisty*pi
    ntz = nsysparams.ntwistz*pi

    nsites = nx*ny*nz

    nLx = nx*alpha_n
    nLy = ny*alpha_n
    nLz = nz*alpha_n

    #Here I assume 1D,2D,3D starts in the order of nx,ny,nz.
    #E.g. if you need a 4x4 lattice put nx=4,ny=4,nz=1 in the Sample.jl
    if nx %2 == 0 && ny %2 == 1 && nz %2 == 1
    #this one is for nx even ny nz odd lattices. Works for 1D and 2D as well.(2x1x1, 4x3x1)

        nxs = [i*nLx/nx for i in ((-(nx/2) +1):((nx/2)))]
        nys = [i*nLy/ny for i in (-(ny-1)/2:(ny-1)/2)]
        nzs = [i*nLz/nz for i in (-(nz-1)/2:(nz-1)/2)]

        nRm =[((i),(j),(l)) for i in (nxs) for j in (nys) for l in (nzs)]

        nRj = nRm
    end
    if nx %2 == 1 && ny %2 == 0 && nz %2 == 1
    #this one is for nx odd ny even lattice. Works for 1D,2D as well. (1x4x1,3x4x1)
        nxs = [i*nLx/nx for i in (-(nx-1)/2:(nx-1)/2)]
        nys = [i*nLy/ny for i in ((-(ny/2) +1):((ny/2)))]
        nzs = [i*nLz/nz for i in (-(nz-1)/2:(nz-1)/2)]

        nRm =[((i),(j),(l)) for i in (nxs) for j in (nys) for l in (nzs)]

        nRj = nRm
    end

    if nx %2 == 0 && ny %2 == 0 && nz %2 == 1
        #this one is for nx,ny even nz odd lattice. Works for 2D as well.(4x4x1)

            nxs = [i*nLx/nx for i in ((-(nx/2) +1):((nx/2)))]
            nys = [i*nLy/ny for i in ((-(ny/2) +1):((ny/2)))]
            nzs = [i*nLz/nz for i in (-(nz-1)/2:(nz-1)/2)]

            nRm =[((i),(j),(l)) for i in (nxs) for j in (nys) for l in (nzs)]

            nRj = nRm
    end
    if nx %2 == 0 && ny %2 == 0 && nz %2 == 0
        #3D even lattice

        # nxs = [i*nLx/nx for i in ((-(nx/2) +1):((nx/2)))]
        # nys = [i*nLy/ny for i in ((-(ny/2) +1):((ny/2)))]
        # nzs = [i*nLz/nz for i in ((-(nz/2) +1):((nz/2)))]

        nxs = [i*nLx/nx for i in ((-(nx/2) + 0.5):((nx/2) - 0.5))]
        nys = [i*nLy/ny for i in ((-(ny/2) + 0.5):((ny/2) - 0.5))]
        nzs = [i*nLz/nz for i in ((-(nz/2) + 0.5):((nz/2) - 0.5))]

        nRm =[((i),(j),(l)) for i in (nxs) for j in (nys) for l in (nzs)]

        nRj = nRm
        # @show nRj
    end

    if nx %2 == 1 && ny %2 == 1 && nz %2 == 1
        #1D 2D or 3D odd lattices

        nxs = [i*nLx/nx for i in (-(nx-1)/2:(nx-1)/2)]
        nys = [i*nLy/ny for i in (-(ny-1)/2:(ny-1)/2)]
        nzs = [i*nLz/nz for i in (-(nz-1)/2:(nz-1)/2)]

        nRm =[((i),(j),(l)) for i in (nxs) for j in (nys) for l in (nzs)]

        nRj = nRm
    end
    #@show size(Rj)
    nkxs_up = [((2*pi)*i + ntx)/nLx for i in (nxs)/alpha_n]
    nkys_up = [((2*pi)*i + nty)/nLy for i in (nys)/alpha_n]
    nkzs_up = [((2*pi)*i + ntz)/nLz for i in (nzs)/alpha_n]

    # @show nkxs_up
    if nA == 0
        nks_up = [(0,0,0) for i in (nkxs_up) for j in (nkys_up) for l in (nkzs_up)]
    else
        nks_up = [((i),(j),(l)) for i in (nkxs_up) for j in (nkys_up) for l in (nkzs_up)]
    end

    nkxs_do = [((2*pi)*i - ntx)/nLx for i in (nxs)/alpha_n]
    nkys_do = [((2*pi)*i - nty)/nLy for i in (nys)/alpha_n]
    nkzs_do = [((2*pi)*i - ntz)/nLz for i in (nzs)/alpha_n]

    if nB == 0
        nks_do = [(0,0,0) for i in (nkxs_do) for j in (nkys_do) for l in (nkzs_do)]
    else
        nks_do = [((i),(j),(l)) for i in (nkxs_do) for j in (nkys_do) for l in (nkzs_do)]
    end

    nH_up =zeros(Complex{Float64}, nsites, nsites)
    nH_do =zeros(Complex{Float64}, nsites, nsites)

    for (jj,rj) in enumerate(nRj)

        for (mm,rm) in enumerate(nRm)

            t_jm_up = 0
            t_jm_do = 0

            for kn in (nks_up)
                # @show kn
                product = dot(kn,rj) - dot(kn,rm)
                t_jm_up += Adisp(nsysparams,kn)*exp(product*im)
            end

            t_jm_up = t_jm_up/nsites
            nH_up[jj,mm] =t_jm_up

            for kn in (nks_do)
                product = dot(kn,rj) - dot(kn,rm)
                t_jm_do += Bdisp(nsysparams,kn)*exp(product*im)
            end
            t_jm_do = t_jm_do/nsites
            nH_do[jj,mm] = t_jm_do


        end


    end


    return nH_up,nH_do



end

function ph_k(nsysparams::nParamsSys)

    #protons
    nC = nsysparams.nC
    nD = nsysparams.nD
    px = nsysparams.nx
    py = nsysparams.ny
    pz = nsysparams.nz

    alpha_n = nsysparams.alpha_n

    ptx =nsysparams.ntwistx*pi
    pty = nsysparams.ntwisty*pi
    ptz = nsysparams.ntwistz*pi

    psites = px*py*pz

    pLx = px*alpha_n
    pLy = py*alpha_n
    pLz = pz*alpha_n

    #Here I assume 1D,2D,3D starts in the order of px,py,pz.
    #E.g. if you need a 4x4 proton lattice put px=4,py=4,pz=1 in the Sample.jl
    if px %2 == 0 && py %2 == 1 && pz %2 == 1
    #this one is for px even py pz odd lattices. Works for 1D and 2D as well.(2x1x1, 4x3x1)

        pxs = [i*pLx/px for i in ((-(px/2) +1):((px/2)))]
        pys = [i*pLy/py for i in (-(py-1)/2:(py-1)/2)]
        pzs = [i*pLz/pz for i in (-(pz-1)/2:(pz-1)/2)]

        pRm =[((i),(j),(l)) for i in (pxs) for j in (pys) for l in (pzs)]

        pRj = pRm
    end
    if px %2 == 1 && py %2 == 0 && pz %2 == 1
    #this one is for px odd py even lattice. Works for 1D,2D as well. (1x4x1,3x4x1)
        pxs = [i*pLx/px for i in (-(px-1)/2:(px-1)/2)]
        pys = [i*pLy/py for i in ((-(py/2) +1):((py/2)))]
        pzs = [i*pLz/pz for i in (-(pz-1)/2:(pz-1)/2)]

        pRm =[((i),(j),(l)) for i in (pxs) for j in (pys) for l in (pzs)]

        pRj = pRm
    end

    if px %2 == 0 && py %2 == 0 && pz %2 == 1
        #this one is for px,py even pz odd lattice. Works for 2D as well.(4x4x1)

            pxs = [i*pLx/px for i in ((-(px/2) +1):((px/2)))]
            pys = [i*pLy/py for i in ((-(py/2) +1):((py/2)))]
            pzs = [i*pLz/pz for i in (-(pz-1)/2:(pz-1)/2)]

            pRm =[((i),(j),(l)) for i in (pxs) for j in (pys) for l in (pzs)]

            pRj = pRm
    end
    if px %2 == 0 && py %2 == 0 && pz %2 == 0
        #3D even lattice

        pxs = [i*pLx/px for i in ((-(px/2) +1):((px/2)))]
        pys = [i*pLy/py for i in ((-(py/2) +1):((py/2)))]
        pzs = [i*pLz/pz for i in ((-(pz/2) +1):((pz/2)))]

        pRm =[((i),(j),(l)) for i in (pxs) for j in (pys) for l in (pzs)]

        pRj = pRm
    end

    if px %2 == 1 && py %2 == 1 && pz %2 == 1
        #1D 2D or 3D odd lattices

        pxs = [i*pLx/px for i in (-(px-1)/2:(px-1)/2)]
        pys = [i*pLy/py for i in (-(py-1)/2:(py-1)/2)]
        pzs = [i*pLz/pz for i in (-(pz-1)/2:(pz-1)/2)]

        pRm =[((i),(j),(l)) for i in (pxs) for j in (pys) for l in (pzs)]

        pRj = pRm
    end
    #@show size(pRj)
    pkxs_up = [((2*pi)*i + ptx)/pLx for i in (pxs)/alpha_n]
    pkys_up = [((2*pi)*i + pty)/pLy for i in (pys)/alpha_n]
    pkzs_up = [((2*pi)*i + ptz)/pLz for i in (pzs)/alpha_n]
    # @show pkxs_up
    if nC == 0
        pks_up = [(0,0,0) for i in (pkxs_up) for j in (pkys_up) for l in (pkzs_up)]
    else
        pks_up = [((i),(j),(l)) for i in (pkxs_up) for j in (pkys_up) for l in (pkzs_up)]
    end
    pkxs_do = [((2*pi)*i - ptx)/pLx for i in (pxs)/alpha_n]
    pkys_do = [((2*pi)*i - pty)/pLy for i in (pys)/alpha_n]
    pkzs_do = [((2*pi)*i - ptz)/pLz for i in (pzs)/alpha_n]

    if nD == 0
        pks_do = [(0,0,0) for i in (pkxs_do) for j in (pkys_do) for l in (pkzs_do)]
    else
        pks_do = [((i),(j),(l)) for i in (pkxs_do) for j in (pkys_do) for l in (pkzs_do)]
    end


    pH_up =zeros(Complex{Float64}, psites, psites)
    pH_do =zeros(Complex{Float64}, psites, psites)

    for (jj,rj) in enumerate(pRj)

        for (mm,rm) in enumerate(pRm)

            t_jm_up = 0
            t_jm_do = 0

            for kp in (pks_up)
                product = dot(kp,rj) - dot(kp,rm)
                t_jm_up += Cdisp(nsysparams,kp)*exp(product*im)
            end

            t_jm_up = t_jm_up/psites
            pH_up[jj,mm] =t_jm_up

            for kp in (pks_do)
                product = dot(kp,rj) - dot(kp,rm)
                t_jm_do += Ddisp(nsysparams,kp)*exp(product*im)
            end
            t_jm_do = t_jm_do/psites
            pH_do[jj,mm] = t_jm_do


        end


    end

    return pH_up,pH_do
end



# runparams = ParamsRun()
# nsysparams = nParamsSys()

# ek_A = Adisp(nsysparams::nParamsSys,kn)
# ek_B = Bdisp(nsysparams::nParamsSys,kn)
# ek_C = Cdisp(nsysparams::nParamsSys,kn)
# ek_D = Ddisp(nsysparams::nParamsSys,kn)
# nH_up,nH_do = nh_k(nsysparams)
# pH_up,pH_do = ph_k(nsysparams)
# @show size(pH_up)
# @show eigvals(pH_do)
# @show eigvals(nH_do)
# @show eigvals(pH_up)
# @show eigvals(nH_up)

end
