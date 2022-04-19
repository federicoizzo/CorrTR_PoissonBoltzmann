function PB_gen_shape_WENO(N;
    w_ker::Function = Kinf, # averaging kernel
    surfTargetXYZ::Array{Array{Float64,1},1}=[[1;0.;0]], # surface targets: identified via 3D point, then projected
    outflag::Bool=true,  # interior or exterior problem: true for exterior
    shift::Array{Float64,1}=[0.;0.;0.], # shift of grid w.r.t. surface center (if any)
    epslI::Real=1.0, epslE::Real=1.0, kappa_val::Real=1.0,
    fε::Function=(t->2t), # function to get espilon from h
    plotting_surface::Bool=false, count::Int64=1
    )

    R=1.0
    R1 = 0.7; R2 = 0.2; 
    Rcil = 0.4; Zcil = 0.6;
    ZCIL = [0;0;Zcil]

    pl=4;
    h = 1.1*(2R+0.2)/(N)
        println(1.1*(2R+0.2)/(N))
    Xlim = [(-2.3)*1.1-pl*h; (1.5)*1.1+pl*h]
    Ylim = [(-1.5)*1.1-pl*h; (1.5)*1.1+pl*h]
    Zlim = [(-1.5)*1.1-pl*h; (1.5)*1.1+pl*h]
    
    ε = fε(h)
    # println("epsl=",ε)
    
    Rx(a) = [1 0 0;0 cos(a) sin(a);0 -sin(a) cos(a)]
    Ry(a) = [cos(a) 0 sin(a);0 1 0;-sin(a) 0 cos(a)]
    Rz(a) = [cos(a) -sin(a) 0;sin(a) cos(a) 0;0 0 1]
    xshift = zeros(3)

    # a,b,c = 0.24402412255508432, 0.7454097947651017, 2.219760487439292
    # a,b,c = 1.99487, 2.54097947651017, 4.219760487439292
    a,b,c = 0., 0., 0.
    
    # a,b,c = 0., 0., 0.
    Qrot = Rz(c)*Ry(b)*Rx(a);
    Rinv = Rx(-a)*Ry(-b)*Rz(-c);
    x0 = xshift; #zeros(3); 

    if false # sphere surface functions
        function Pgammafun(z)
            return (norm(z-x0)>1e-8)*((z-x0)/norm(z-x0)*R+x0) + (norm(z-x0)<=1e-8)*(x0+R*[0.;0;1.]);
        end
        function normalz(z)
            z = Pgammafun(z);
            return (norm(z-x0)>1e-8)*((z-x0)/norm(z-x0)) + (norm(z-x0)<=1e-8)*([0;0;1.]);
        end
        function insidepoint(z)
            z = Pgammafun(z)
            return norm(z-x0)<R
        end
        function far_insidepoint(z)
            z = Pgammafun(z)
            return norm(z-x0)<R-ε
        end
        function far_outsidepoint(z)
            z = Pgammafun(z)
            return norm(z-x0)>R+ε
        end
    end
    if false # torus surface functions
        torus(t,p) = [(R2*cos(t)+R1)*cos(p);(R2*cos(t)+R1)*sin(p);R2*sin(t)];
        normalt(t,p) = [cos(t)*cos(p);cos(t)*sin(p);sin(t)];
        #torus
        function Pgammafun(z)
            zn = Rinv*(z.-xshift);
            cosp, sinp = [zn[1];zn[2]]./sqrt(zn[1]^2+zn[2]^2)
            # phi = atan(zn[2],zn[1]);
            ctmp = R1*[cosp;sinp;0];
            return Qrot*((zn.-ctmp)/norm(zn.-ctmp)*R2.+ctmp).+xshift
        end
        function normalz(z)
            z = Pgammafun(z)
            zn = Rinv*(z.-xshift);
            cosp, sinp = [zn[1];zn[2]]./sqrt(zn[1]^2+zn[2]^2)
            ctmp = R1*[cosp;sinp;0];
            sint = R2/zn[3];
            cost = (abs(cosp)>1e-1)*(zn[1]/R2/cosp-R1/R2)+(abs(cosp)<=1e-1)*(zn[2]/R2/sinp-R1/R2)
            sint = (1-2*(sint<0))*sqrt(1-round(cost,digits=12)^2)
            v = (Qrot*[cost*cosp;cost*sinp;sint])
            return -v/norm(v)
        end
        function insidepoint(z)
            zn = Rinv*(z.-xshift);
            cosp, sinp = [zn[1];zn[2]]./sqrt(zn[1]^2+zn[2]^2)
            ctmp = R1*[cosp;sinp;0];
            return norm(ctmp-zn)<R2
        end
        function far_insidepoint(z)
            zn = Rinv*(z.-xshift);
            cosp, sinp = [zn[1];zn[2]]./sqrt(zn[1]^2+zn[2]^2)
            ctmp = R1*[cosp;sinp;0];
            return norm(ctmp-zn)<R2-ε
        end
        function far_outsidepoint(z)
            zn = Rinv*(z.-xshift);
            cosp, sinp = [zn[1];zn[2]]./sqrt(zn[1]^2+zn[2]^2)
            ctmp = R1*[cosp;sinp;0];
            return norm(ctmp-zn)>R2+ε
        end
        function Pgammafun_param(z) # for the torus, parameters theta and phi
            z=Pgammafun(z)
            zn = Rinv*(z.-xshift);
            cosp, sinp = [zn[1];zn[2]]./sqrt(zn[1]^2+zn[2]^2)
            ctmp = R1*[cosp;sinp;0];
            sint = R2/zn[3];
            cost = (abs(cosp)>1e-1)*(zn[1]/R2/cosp-R1/R2)+(abs(cosp)<=1e-1)*(zn[2]/R2/sinp-R1/R2)
            sint = (1-2*(sint<0))*sqrt(1-round(cost,digits=12)^2)
            return [cost;sint;cosp;sinp]
        end
    end
    if false # cilinder with hemispheres attached (pill or capsule-like shape) surface functions
        # #cilinder
        function Pgammafun(z)
            zn = Rinv*(z.-xshift);
            if abs(zn[3]) <= Zcil
                xy = (norm(zn[1:2])>1e-5)*(Rcil*zn[1:2]/norm(zn[1:2])) .+ (norm(zn[1:2])<=1e-5)*(Rcil*[0.;1])
                return Qrot*([xy[1];xy[2];zn[3]]) .+ xshift
            elseif zn[3]>Zcil
                xyz = (norm(zn .- ZCIL)>1e-5)*(Rcil*(zn .- ZCIL)/norm(zn .- ZCIL)) .+ (norm(zn .- ZCIL)<=1e-5)*(Rcil*[0.;0;1]) .+ ZCIL
                return Qrot*xyz .+ xshift
            else
                xyz = (norm(zn .+ ZCIL)>1e-5)*(Rcil*(zn .+ ZCIL)/norm(zn .+ ZCIL))  .+  (norm(zn .+ ZCIL)<=1e-5)*(Rcil*[0.;0;-1]) .- ZCIL
                return Qrot*xyz .+ xshift
            end
        end
        function normalz(z)
            z = Pgammafun(z);
            zn = Rinv*(z.-xshift);
            if abs(zn[3]) <= Zcil
                local xy = (zn[1:2]/norm(zn[1:2]))
                return Qrot*([xy[1];xy[2];0]/norm([xy[1];xy[2];0]))
            elseif zn[3]>Zcil
                local xyz = ((zn .- ZCIL)/norm(zn .- ZCIL))
                return Qrot*(xyz/norm(xyz))
            else
                local xyz = ((zn .+ ZCIL)/norm(zn .+ ZCIL))
                return Qrot*(xyz/norm(xyz))
            end
        end
        function insidepoint(z)
            # z = Pgammafun(z);
            zn = Rinv*(z.-xshift);
            if abs(zn[3]) <= Zcil
                return norm(zn[1:2])<Rcil
            elseif zn[3]>Zcil
                return norm(zn .- ZCIL)<Rcil
            else
                return norm(zn .+ ZCIL)<Rcil
            end
        end
        function far_insidepoint(z)
            # z = Pgammafun(z);
            zn = Rinv*(z.-xshift);
            if abs(zn[3]) <= Zcil
                return norm(zn[1:2])<Rcil-ε
            elseif zn[3]>Zcil
                return norm(zn .- ZCIL)<Rcil-ε
            else
                return norm(zn .+ ZCIL)<Rcil-ε
            end
        end
        function far_outsidepoint(z)
            # z = Pgammafun(z);
            zn = Rinv*(z.-xshift);
            if abs(zn[3]) <= Zcil
                return norm(zn[1:2])>Rcil+ε
            elseif zn[3]>Zcil
                return norm(zn .- ZCIL)>Rcil+ε
            else
                return norm(zn .+ ZCIL)>Rcil+ε
            end
        end
        function get_jac(z,η)
            # η = (1-2*insidepoint(z))*norm(z-Pgammafun(z));
            z = Pgammafun(z);
            zn = Rinv*(z.-xshift);
            # return Rcil*(Rcil+Zcil)/((Rcil+η)*(Rcil+η+Zcil))
            if abs(zn[3]) <= Zcil
                return 1-η/((Rcil+η))
            else
                return 1-2η/(Rcil+η)+η^2/(Rcil+η)^2
            end
        end
    end

    Pgammafun = Pgammafun_spheres
    insidepoint = insidepoint_spheres
    get_jac = get_jac_spheres
    far_insidepoint = far_insidepoint_spheres
    far_outsidepoint = far_outsidepoint_spheres

    Nepsl2h=Int(ceil(2*ε/h)) # how many discretization points fall in the tubular neighborhood with these values of ε and h
    h2 = h*h
    h3 = h2*h

    # println("Number of surface targets = $(length(surfTargetXYZ))")
    n_surf_trg = length(surfTargetXYZ)
    for i=1:n_surf_trg
        surfTargetXYZ[i] = Pgammafun(surfTargetXYZ[i])
    end
    
    rshift = xshift+shift;
    # discretizations in each direction, centered in rshift, with step h
    X = [(rshift[1]:(-h):Xlim[1])[end:-1:1]; (rshift[1]+h):h:Xlim[2]]; Nx = length(X)
    Y = [(rshift[2]:(-h):Ylim[1])[end:-1:1]; (rshift[2]+h):h:Ylim[2]]; Ny = length(Y)
    Z = [(rshift[3]:(-h):Zlim[1])[end:-1:1]; (rshift[3]+h):h:Zlim[2]]; Nz = length(Z)
    println("Computational box limits: X=$([X[1];X[end]]), Y=$([Y[1];Y[end]]), Z=$([Z[1];Z[end]])")

    Nvec = [Nx;Ny;Nz]

    println("Current run: h=$h, ε=$ε, Nepsl=$Nepsl2h")
    println("Number of discretization points in each direction:\n$Nx x $Ny x $Nz = $(Nx*Ny*Nz)")

    tau = h # cut-off parameter for function regularization in IBIMs
        
    epsl_ratio = epslE/epslI;
    @time M, Pg, source, indIJK_to_M, indM_to_IJK, dsignes, ab_single_t, iv1_t, w_K11_single_t, w_K22_single_t, w_K21_single_t, w_K12_single_t, target_normal = genCPM_corr_V2_PB(Pgammafun, insidepoint, far_insidepoint, far_outsidepoint, X, Y, Z, Nvec, ε, h; outflag=outflag, surfTargetXYZ=surfTargetXYZ, epsl_ratio=epsl_ratio, kappa_val=kappa_val)

    dPg = zeros(3,3)
    salvare_jac_2nd = Array{Float64}(undef,M);
    salvare_ker = Array{Float64}(undef,M);
    normal = Array{Array{Float64,1},1}(undef,M)
    β = Array{Float64}(undef,M)

    mysmoothf(Z) = 1.12678+sin(Z[1])*Z[2]+cos(Z[3])*0.2376*(Z[1]-0.3176)
    @time for m=1:M # for every node inside Tε, compute normal and jacobian
        i,j,k = indM_to_IJK[m,:]; # get ijk-indices from m-index
        # y = [ X[i]; Y[j]; Z[k] ] # node in 3D volume
        
        # computing the normal with potentially one-sided 2nd order scheme
        iind, i1 = get_indices_b(dsignes[i-2:i+2,j,k])
        jind, j1 = get_indices_b(dsignes[i,j-2:j+2,k])
        kind, k1 = get_indices_b(dsignes[i,j,k-2:k+2])
        normal[m] = [dot(dsignes[iind.+i,j,k],i1); dot(dsignes[i,jind.+j,k],j1);dot(dsignes[i,j,kind.+k],k1)]/h
        normal[m] /= norm(normal[m])

        # computing the Jacobian with potentially one-sided 2nd order scheme
        ind1, i11 = get_indices_b(Pg[i-2:i+2,j,k,1])
        ind2, i12 = get_indices_b(Pg[i-2:i+2,j,k,2])
        ind3, i13 = get_indices_b(Pg[i-2:i+2,j,k,3])
        dPg[1,1] = dot(i11,Pg[ind1.+i,j,k,1]);
        dPg[2,1] = dot(i12,Pg[ind2.+i,j,k,2]);
        dPg[3,1] = dot(i13,Pg[ind3.+i,j,k,3]);
        ind1, i11 = get_indices_b(Pg[i,j-2:j+2,k,1])
        ind2, i12 = get_indices_b(Pg[i,j-2:j+2,k,2])
        ind3, i13 = get_indices_b(Pg[i,j-2:j+2,k,3])
        dPg[1,2] = dot(i11,Pg[i,ind1.+j,k,1]);
        dPg[2,2] = dot(i12,Pg[i,ind2.+j,k,2]);
        dPg[3,2] = dot(i13,Pg[i,ind3.+j,k,3]);
        ind1, i11 = get_indices_b(Pg[i,j,k-2:k+2,1])
        ind2, i12 = get_indices_b(Pg[i,j,k-2:k+2,2])
        ind3, i13 = get_indices_b(Pg[i,j,k-2:k+2,3])
        dPg[1,3] = dot(i11,Pg[i,j,ind1.+k,1]);
        dPg[2,3] = dot(i12,Pg[i,j,ind2.+k,2]);
        dPg[3,3] = dot(i13,Pg[i,j,ind3.+k,3]);
        dPg /= h

        jac = svdvals(dPg) # SVD values of numerical Jacobian
        salvare_jac_2nd[m] = jac[1]*jac[2]

        salvare_ker[m] = (abs(dsignes[i,j,k])<ε)*w_ker(dsignes[i,j,k]/ε)/ε # zero outside the tubular neighborhood
        
        β[m] = mysmoothf(source[m]) # given density
    end
    
    salvare_wo = salvare_ker # missing h3
    salvare_w2 = salvare_ker.*salvare_jac_2nd # missing h3
    salvareu = salvare_w2
    
    surf_val = zeros(2)
    surf_val[1] = sum(salvare_wo)*h3
    surf_val[2] = sum(salvare_w2)*h3
    
    println("Val. surface area:\n$(surf_val)")

    # IBIM is with Jacobian=1
    println("Evaluating the potentials in IBIM (4) and CTR (4)")
    @time begin
        K11_IBIM = K11_PB_IBIM_target(β; source=source, normal=normal, salvare=salvare_wo, targets=surfTargetXYZ, tau=tau, kappa_val=kappa_val, theta_val=epsl_ratio)*h3
        K22_IBIM = K22_PB_IBIM_target(β; source=source, normal=normal, targetnormals=target_normal, salvare=salvare_wo, targets=surfTargetXYZ, tau=tau, kappa_val=kappa_val, theta_val=1/epsl_ratio)*h3
        K21_IBIM = K21_PB_IBIM_target(β; source=source, normal=normal, targetnormals=target_normal, salvare=salvare_wo, targets=surfTargetXYZ, tau=tau, kappa_val=kappa_val)*h3
        K12_IBIM = K12_PB_IBIM_target(β; source=source, salvare=salvare_wo, targets=surfTargetXYZ, tau=tau, kappa_val=kappa_val)*h3

        K11_corr = K11_PB_Q1corr_target(β; source=source, normal=normal, salvare=salvareu, targets=surfTargetXYZ, kappa_val=kappa_val, theta_val=epsl_ratio, h=h, iv1=iv1_t, wv1=w_K11_single_t)
        K22_corr = K22_PB_Q1corr_target(β; source=source, targetnormals=target_normal, salvare=salvareu, targets=surfTargetXYZ, kappa_val=kappa_val, theta_val=1/epsl_ratio, h=h, iv1=iv1_t, wv1=w_K22_single_t)
        K21_corr = K21_PB_Q1corr_target(β; source=source, normal=normal,  targetnormals=target_normal, salvare=salvareu, targets=surfTargetXYZ, kappa_val=kappa_val, h=h, iv1=iv1_t, wv1=w_K21_single_t)
        K12_corr = K12_PB_Q1corr_target(β; source=source, salvare=salvareu, targets=surfTargetXYZ, kappa_val=kappa_val, h=h, iv1=iv1_t, wv1=w_K12_single_t)
    end

    if ~false
        val_abs = zeros(8,n_surf_trg)
        val_abs[1,:] = K11_IBIM
        val_abs[2,:] = K22_IBIM
        val_abs[3,:] = K21_IBIM
        val_abs[4,:] = K12_IBIM
        val_abs[5,:] = K11_corr
        val_abs[6,:] = K22_corr
        val_abs[7,:] = K21_corr
        val_abs[8,:] = K12_corr

        val_err = zeros(4,n_surf_trg)
        # val_err[1,:] = abs.(λ1*ψ[1] .+ K11_IBIM_2*ψ[1] .- K12_IBIM_2*ψn[1] .- g1)./abs(g1)
        # val_err[2,:] = abs.(λ2*ψn[1] .+ K21_IBIM_2*ψ[1] .- K22_IBIM_2*ψn[1] .- g2)./abs(g2)
        # val_err[3,:] = abs.(λ1*ψ[1] .+ K11_corr*ψ[1] .- K12_corr*ψn[1] .- g1)./abs(g1)
        # val_err[4,:] = abs.(λ2*ψn[1] .+ K21_corr*ψ[1] .- K22_corr*ψn[1] .- g2)./abs(g2)

        mval_abs = mean(val_abs, dims=2)
        mval_err = mean(val_err, dims=2)
        println("Average of potentials over targets:\n",mval_abs[1:8])
        # println(mval_err)
    end

    if plotting_surface
        # println(surfTargetXYZ[1])
        # println(surfTargetXYZ)
        xp = [ sour[1] for sour in source ];
        yp = [ sour[2] for sour in source ];
        zp = [ sour[3] for sour in source ];
        xt = [ t[1] for t in surfTargetXYZ ];
        yt = [ t[2] for t in surfTargetXYZ ];
        zt = [ t[3] for t in surfTargetXYZ ];
        figure(50); clf()

        # c=scatter3D( xp, yp, zp, s=1, c=log10.(abs.(β)), marker=".")
        scatter3D( xp, yp, zp, s=1, c="k", marker=".")
        # colorbar(c)

        
        scatter3D( xt, yt, zt, s=3, c="r", marker=".")
        # colorbar(c)
        
        xlim([-1.0;1])
        ylim([-1.0;1])
        zlim([-1.0;1])
    end
 
    x=[ R, R1, R2]
    open(newdir_nrun*"/surfacePB_R_R1_R2_Rcil_data_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end
    
    x=[ a, b, c, xshift ]
    open(newdir_nrun*"/surfacePB_abc_xshift_data_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end
    x=[ kappa_val, epslE, epslI, R]
    open(newdir_nrun*"/surfacePB_kappa_epslEI_data_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end

    x=[ val_abs ]
    open(newdir_nrun*"/surfacePB_val_abs_data_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end
    x= size(val_abs)
    open(newdir_nrun*"/surfacePB_val_abs_size_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end
    x=[ val_err ]
    open(newdir_nrun*"/surfacePB_val_err_data_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end
    x= size(val_err)
    open(newdir_nrun*"/surfacePB_val_err_size_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end

    x=[ surf_val ]
    open(newdir_nrun*"/surfacePB_surf_val_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end

    x=[ h, ε, Nvec, 2*ε/h, tau ]
    open(newdir_nrun*"/surfacePB_h_epsl_data_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end

    println("Run over.")

    return h, val_abs, val_err, mval_abs, mval_err, surf_val
end