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
    @time M, Pg, source, indIJK_to_M, indM_to_IJK, dsignes, ab_single_t, iv1_t, w_K11_single_t, w_K22_single_t, w_K21_single_t, target_normal = genCPM_corr_V2_PB(Pgammafun, insidepoint, far_insidepoint, far_outsidepoint, X, Y, Z, Nvec, ε, h; outflag=outflag, surfTargetXYZ=surfTargetXYZ, epsl_ratio=epsl_ratio, kappa_val=kappa_val)

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
        K12_corr = K12_PB_Q1corr_target(β; source=source, salvare=salvareu, targets=surfTargetXYZ, kappa_val=kappa_val, h=h, iv1=iv1_t)
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

function PB_gen_shape_plane(N;
    w_ker::Function = Kcos, # averaging kernel
    surfTargetXYZ::Array{Array{Float64,1},1}=[[1;0.;0]], # surface targets: identified via 3D point, then projected
    outflag::Bool=true,  # interior or exterior problem: true for exterior
    shift::Array{Float64,1}=[0.;0.;0.], # shift of grid w.r.t. surface center (if any)
    epslI::Real=1.0, epslE::Real=1.0, kappa_val::Real=1.0,
    fε::Function=(t->2t), # function to get espilon from h
    plotting_surface::Bool=false, count::Int64=1
    )

    R=1.0
    # R1 = 0.7; R2 = 0.2; 
    R1 = 1.0; R2 = 0.7; 
    Rcil = 0.4; Zcil = 0.6;
    ZCIL = [0;0;Zcil]

    pl=4;
    h = 1.1*(2R+0.2)/(N)
    println(1.1*(2R+0.2)/(N))
    Xlim = [(-2.3)*1.1-pl*h; (1.5)*1.1+pl*h]
    Ylim = [(-1.5)*1.1-pl*h; (1.5)*1.1+pl*h]
    # Zlim = [(-1.5)*1.1-pl*h; (1.5)*1.1+pl*h]
    hei = 0.5
    Zlim = [(hei)-pl*h; (hei)+pl*h]
    
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
    Z = Array{Float64,1}(Zlim[1]:h:Zlim[2]); Nz = length(Z)
    println(Z)
    # Z = [1.05] 
    kz = 5
    println("Computational box limits: X=$([X[1];X[end]]), Y=$([Y[1];Y[end]]), Z=$([Z[1];Z[end]])")

    Nvec = [Nx;Ny;Nz]

    println("Current run: h=$h, ε=$ε, Nepsl=$Nepsl2h")
    println("Number of discretization points in each direction:\n$Nx x $Ny x $Nz = $(Nx*Ny*Nz)")

    tau = h # cut-off parameter for function regularization in IBIMs
        
    epsl_ratio = epslE/epslI;
    @time M, Pg, source, indIJK_to_M, indM_to_IJK, dsignes, ab_single_t, iv1_t, w_K11_single_t, w_K22_single_t, w_K21_single_t, target_normal, target_kappa = genCPM_corr_PB_plane(Pgammafun, insidepoint, far_insidepoint, far_outsidepoint, X, Y, Z, kz, Nvec, ε, h; outflag=outflag, surfTargetXYZ=surfTargetXYZ, epsl_ratio=epsl_ratio, kappa_val=kappa_val)

    dPg = zeros(3,3)
    salvare_jac_2nd = Array{Float64}(undef,M);
    salvare_ker = Array{Float64}(undef,M);
    normal = Array{Array{Float64,1},1}(undef,M)
    # errk2 = Array{Array{Float64,1},1}(undef,M)
    β = Array{Float64}(undef,M)
    errk2 = Array{Float64}(undef,M)
    errk3 = Array{Float64}(undef,M)
    k1n = Array{Float64}(undef,M)
    k2n = Array{Float64}(undef,M)
    M2 = [0]

    mysmoothf(Z) = 1.12678+sin(Z[1])*Z[2]+cos(Z[3])*0.2376*(Z[1]-0.3176)
    @time for m=1:M # for every node inside Tε, compute normal and jacobian
        i,j,k = indM_to_IJK[m,:]; # get ijk-indices from m-index
        y = [ X[i]; Y[j]; Z[k] ] # node in 3D volume
        
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

        # ka1, ka2 = get_kappa_spheres(y,dsignes[i,j,k])
        ka1, ka2 = get_kappa_spheres(y,0.)
        Jm = hessian3D_2nd(dsignes[i-1:i+1,j-1:j+1,k-1:k+1], h, h, h)
        eigvalF, eigvecF = eigen(Jm);
        eigvp = sortperm(abs.(eigvalF))[3:-1:1]
        # find first principal direction
        tau1 = eigvecF[:,eigvp[1]]; tau1 /= norm(tau1) # normalization
        # find second principal direction by using the normal
        tau2 = cross(normal[m], tau1); tau2 /= norm(tau2)
        k1_now = -eigvalF[eigvp[1]] # first eigenvalue (from the formula it comes with opposite sign)
        k2_now = -eigvalF[eigvp[2]] 
        # println("k1=$k1_now, k2=$k2_now, k=$(-eigvalF), eta=$(dsignes[i,j,k])")
        k1_now = k1_now/(1+dsignes[i,j,k]*k1_now)
        k2_now = k2_now/(1+dsignes[i,j,k]*k2_now)
        errk = abs.([ka1-k1_now; ka2-k2_now])
        # errk2[m] = norm([ka1-k1_now; ka2-k2_now],Inf)
        errk2[m] = abs(ka1-k1_now)
        errk3[m] = abs(ka2-k2_now)
        k1n[m] = ka1
        k2n[m] = ka2
        if norm(errk,Inf)>1e-1
            M2[1] += 1
            println("error = $errk")
            println("anl.: $ka1, $ka2")
            println("num.: $k1_now, $k2_now\n")
        end

        jac = svdvals(dPg) # SVD values of numerical Jacobian
        salvare_jac_2nd[m] = jac[1]*jac[2]

        salvare_ker[m] = (abs(dsignes[i,j,k])<ε)*w_ker(dsignes[i,j,k]/ε)/ε # zero outside the tubular neighborhood
        
        β[m] = mysmoothf(source[m]) # given density
    end

    println("-------------------------------")
    println(" inside     M = $M")
    println(" bad       M2 = $(M2[1])")
    println("         M/M2 = $(M/M2[1])")
    println(" percent M/M2 = $(round(100*M/M2[1],digits=4))%")
    println("-------------------------------")
    
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
        K22_IBIM = K22_PB_IBIM_target(β; source=source, targetnormals=target_normal, salvare=salvare_wo, targets=surfTargetXYZ, tau=tau, kappa_val=kappa_val, theta_val=1/epsl_ratio)*h3
        K21_IBIM = K21_PB_IBIM_target(β; source=source, normal=normal, targetnormals=target_normal, salvare=salvare_wo, targets=surfTargetXYZ, tau=tau, kappa_val=kappa_val)*h3
        K12_IBIM = K12_PB_IBIM_target(β; source=source, salvare=salvare_wo, targets=surfTargetXYZ, tau=tau, kappa_val=kappa_val)*h3

        K11_corr = K11_PB_Q1corr_target(β; source=source, normal=normal, salvare=salvareu, targets=surfTargetXYZ, kappa_val=kappa_val, theta_val=epsl_ratio, h=h, iv1=iv1_t, wv1=w_K11_single_t)
        K22_corr = K22_PB_Q1corr_target(β; source=source, targetnormals=target_normal, salvare=salvareu, targets=surfTargetXYZ, kappa_val=kappa_val, theta_val=1/epsl_ratio, h=h, iv1=iv1_t, wv1=w_K22_single_t)
        K21_corr = K21_PB_Q1corr_target(β; source=source, normal=normal,  targetnormals=target_normal, salvare=salvareu, targets=surfTargetXYZ, kappa_val=kappa_val, h=h, iv1=iv1_t, wv1=w_K21_single_t)
        K12_corr = K12_PB_Q1corr_target(β; source=source, salvare=salvareu, targets=surfTargetXYZ, kappa_val=kappa_val, h=h, iv1=iv1_t)
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


        figure(51); clf()

        # c=scatter3D( xp, yp, zp, s=1, c=log10.(abs.(β)), marker=".")
        sc=scatter3D( xp, yp, zp, s=2, c=log10.(errk2), marker=".", cmap=:jet, vmin=-5, vmax = 0.5)
        colorbar(sc)
        
        # scatter3D( xp, yp, zp, s=3, c="r", marker=".")
        # colorbar(c)
        
        xlim([-1.0;1])
        ylim([-1.0;1])
        zlim([-1.0;1])
        title("error k1")

        figure(52); clf()

        # c=scatter3D( xp, yp, zp, s=1, c=log10.(abs.(β)), marker=".")
        sc=scatter3D( xp, yp, zp, s=2, c=log10.(errk3), marker=".", cmap=:jet, vmin=-5, vmax = 0.5)
        colorbar(sc)
        
        # scatter3D( xp, yp, zp, s=3, c="r", marker=".")
        # colorbar(c)
        
        xlim([-1.0;1])
        ylim([-1.0;1])
        zlim([-1.0;1])
        title("error k2")


        figure(53); clf()

        # c=scatter3D( xp, yp, zp, s=1, c=log10.(abs.(β)), marker=".")
        sc=scatter3D( xp, yp, zp, s=2, c=k2n, marker=".", cmap=:jet)
        colorbar(sc)
        
        # scatter3D( xp, yp, zp, s=3, c="r", marker=".")
        # colorbar(c)
        
        xlim([-1.0;1])
        ylim([-1.0;1])
        zlim([-1.0;1])

        title("an. k1")
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

function PB_gen_shape_bad_kappas(N;
    w_ker::Function = Kcos, # averaging kernel
    surfTargetXYZ::Array{Array{Float64,1},1}=[[1;0.;0]], # surface targets: identified via 3D point, then projected
    outflag::Bool=true,  # interior or exterior problem: true for exterior
    shift::Array{Float64,1}=[0.;0.;0.], # shift of grid w.r.t. surface center (if any)
    epslI::Real=1.0, epslE::Real=1.0, kappa_val::Real=1.0,
    fε::Function=(t->2t), # function to get espilon from h
    plotting_surface::Bool=false, count::Int64=1
    )

    R=1.0
    # R1 = 0.7; R2 = 0.2; 
    R1 = 1.0; R2 = 0.7; 
    Rcil = 0.4; Zcil = 0.6;
    ZCIL = [0;0;Zcil]

    pl=4;
    h = 1.1*(2R+0.2)/(N)
    println(1.1*(2R+0.2)/(N))
    Xlim = [(-2.4)*1.1-pl*h; (1.6)*1.1+pl*h]
    Ylim = [(-1.5)*1.1-pl*h; (1.5)*1.1+pl*h]
    Zlim = [(-1.5)*1.1-pl*h; (1.5)*1.1+pl*h]
    # hei = 0.5
    # Zlim = [(hei)-pl*h; (hei)+pl*h]
    
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
    # Z = Array{Float64,1}(Zlim[1]:h:Zlim[2]); Nz = length(Z)
    # println(Z)
    # Z = [1.05] 
    # kz = 5
    println("Computational box limits: X=$([X[1];X[end]]), Y=$([Y[1];Y[end]]), Z=$([Z[1];Z[end]])")

    Nvec = [Nx;Ny;Nz]

    println("Current run: h=$h, ε=$ε, Nepsl=$Nepsl2h")
    println("Number of discretization points in each direction:\n$Nx x $Ny x $Nz = $(Nx*Ny*Nz)")

    tau = h # cut-off parameter for function regularization in IBIMs
        
    epsl_ratio = epslE/epslI;
    # @time M, Pg, source, indIJK_to_M, indM_to_IJK, dsignes, ab_single_t, iv1_t, w_K11_single_t, w_K22_single_t, w_K21_single_t, target_normal, target_kappa = genCPM_corr_PB_plane(Pgammafun, insidepoint, far_insidepoint, far_outsidepoint, X, Y, Z, kz, Nvec, ε, h; outflag=outflag, surfTargetXYZ=surfTargetXYZ, epsl_ratio=epsl_ratio, kappa_val=kappa_val)
    @time M, Pg, source, indIJK_to_M, indM_to_IJK, dsignes, ab_single_t, iv1_t, w_K11_single_t, w_K22_single_t, w_K21_single_t, target_normal = genCPM_corr_V2_PB(Pgammafun, insidepoint, far_insidepoint, far_outsidepoint, X, Y, Z, Nvec, ε, h; outflag=outflag, surfTargetXYZ=surfTargetXYZ, epsl_ratio=epsl_ratio, kappa_val=kappa_val)

    dPg = zeros(3,3)
    salvare_jac_2nd = Array{Float64}(undef,M);
    salvare_ker = Array{Float64}(undef,M);
    normal = Array{Array{Float64,1},1}(undef,M)
    # errk2 = Array{Array{Float64,1},1}(undef,M)
    β = Array{Float64}(undef,M)
    errk2 = Array{Float64}(undef,M)
    errk3 = Array{Float64}(undef,M)
    k1n = Array{Float64}(undef,M)
    k2n = Array{Float64}(undef,M)
    M2 = [0]

    mysmoothf(Z) = 1.12678+sin(Z[1])*Z[2]+cos(Z[3])*0.2376*(Z[1]-0.3176)
    @time for m=1:M # for every node inside Tε, compute normal and jacobian
        i,j,k = indM_to_IJK[m,:]; # get ijk-indices from m-index
        y = [ X[i]; Y[j]; Z[k] ] # node in 3D volume
        
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

        # ka1, ka2 = get_kappa_spheres(y,dsignes[i,j,k])
        ka1, ka2 = get_kappa_spheres(y,0.)
        Jm = hessian3D_2nd(dsignes[i-1:i+1,j-1:j+1,k-1:k+1], h, h, h)
        eigvalF, eigvecF = eigen(Jm);
        eigvp = sortperm(abs.(eigvalF))[3:-1:1]
        # find first principal direction
        tau1 = eigvecF[:,eigvp[1]]; tau1 /= norm(tau1) # normalization
        # find second principal direction by using the normal
        tau2 = cross(normal[m], tau1); tau2 /= norm(tau2)
        k1_now = -eigvalF[eigvp[1]] # first eigenvalue (from the formula it comes with opposite sign)
        k2_now = -eigvalF[eigvp[2]] 
        # println("k1=$k1_now, k2=$k2_now, k=$(-eigvalF), eta=$(dsignes[i,j,k])")
        k1_now = k1_now/(1+dsignes[i,j,k]*k1_now)
        k2_now = k2_now/(1+dsignes[i,j,k]*k2_now)
        errk = abs.([ka1-k1_now; ka2-k2_now])
        # errk2[m] = norm([ka1-k1_now; ka2-k2_now],Inf)
        errk2[m] = abs(ka1-k1_now)
        errk3[m] = abs(ka2-k2_now)
        k1n[m] = ka1
        k2n[m] = ka2
        if norm(errk,Inf)>1e-1
            M2[1] += 1
            # println("error = $errk")
            # println("anl.: $ka1, $ka2")
            # println("num.: $k1_now, $k2_now\n")
        end

        jac = svdvals(dPg) # SVD values of numerical Jacobian
        salvare_jac_2nd[m] = jac[1]*jac[2]

        salvare_ker[m] = (abs(dsignes[i,j,k])<ε)*w_ker(dsignes[i,j,k]/ε)/ε # zero outside the tubular neighborhood
        
        β[m] = mysmoothf(source[m]) # given density
    end

    percM = (100*M2[1]/M)
    println("-------------------------------")
    println(" inside     M = $M")
    println(" bad       M2 = $(M2[1])")
    println("         M2/M = $(M2[1]/M)")
    println(" percent M2/M = $(round(percM,digits=4))%")
    println("-------------------------------")
    
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
        K22_IBIM = K22_PB_IBIM_target(β; source=source, targetnormals=target_normal, salvare=salvare_wo, targets=surfTargetXYZ, tau=tau, kappa_val=kappa_val, theta_val=1/epsl_ratio)*h3
        K21_IBIM = K21_PB_IBIM_target(β; source=source, normal=normal, targetnormals=target_normal, salvare=salvare_wo, targets=surfTargetXYZ, tau=tau, kappa_val=kappa_val)*h3
        K12_IBIM = K12_PB_IBIM_target(β; source=source, salvare=salvare_wo, targets=surfTargetXYZ, tau=tau, kappa_val=kappa_val)*h3

        K11_corr = K11_PB_Q1corr_target(β; source=source, normal=normal, salvare=salvareu, targets=surfTargetXYZ, kappa_val=kappa_val, theta_val=epsl_ratio, h=h, iv1=iv1_t, wv1=w_K11_single_t)
        K22_corr = K22_PB_Q1corr_target(β; source=source, targetnormals=target_normal, salvare=salvareu, targets=surfTargetXYZ, kappa_val=kappa_val, theta_val=1/epsl_ratio, h=h, iv1=iv1_t, wv1=w_K22_single_t)
        K21_corr = K21_PB_Q1corr_target(β; source=source, normal=normal,  targetnormals=target_normal, salvare=salvareu, targets=surfTargetXYZ, kappa_val=kappa_val, h=h, iv1=iv1_t, wv1=w_K21_single_t)
        K12_corr = K12_PB_Q1corr_target(β; source=source, salvare=salvareu, targets=surfTargetXYZ, kappa_val=kappa_val, h=h, iv1=iv1_t)
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


        figure(51); clf()

        # c=scatter3D( xp, yp, zp, s=1, c=log10.(abs.(β)), marker=".")
        sc=scatter3D( xp, yp, zp, s=2, c=log10.(errk2), marker=".", cmap=:jet, vmin=-5, vmax = 0.5)
        colorbar(sc)
        
        # scatter3D( xp, yp, zp, s=3, c="r", marker=".")
        # colorbar(c)
        
        xlim([-1.0;1])
        ylim([-1.0;1])
        zlim([-1.0;1])
        title("error k1")

        figure(52); clf()

        # c=scatter3D( xp, yp, zp, s=1, c=log10.(abs.(β)), marker=".")
        sc=scatter3D( xp, yp, zp, s=2, c=log10.(errk3), marker=".", cmap=:jet, vmin=-5, vmax = 0.5)
        colorbar(sc)
        
        # scatter3D( xp, yp, zp, s=3, c="r", marker=".")
        # colorbar(c)
        
        xlim([-1.0;1])
        ylim([-1.0;1])
        zlim([-1.0;1])
        title("error k2")


        figure(53); clf()

        # c=scatter3D( xp, yp, zp, s=1, c=log10.(abs.(β)), marker=".")
        sc=scatter3D( xp, yp, zp, s=2, c=k2n, marker=".", cmap=:jet)
        colorbar(sc)
        
        # scatter3D( xp, yp, zp, s=3, c="r", marker=".")
        # colorbar(c)
        
        xlim([-1.0;1])
        ylim([-1.0;1])
        zlim([-1.0;1])

        title("an. k1")
    end
 
    x=[ R, R1, R2]
    open(newdir_nrun*"/surfacePB_R_R1_R2_Rcil_data_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end

    x=[ M, M2, percM]
    open(newdir_nrun*"/surfacePB_M_M2_data_$(nrun[1])_$(count).dat","w") do io
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

    return h, val_abs, val_err, mval_abs, mval_err, surf_val, [M; M2], percM
end

function PB_gen_shape_molecule_area(N;
    w_ker::Function = Kcos, # averaging kernel
    Xlim, Ylim, Zlim, Nx::Int64=65, h::Real=0.1, dsignes::Array{Float64,3}=zeros(4,4,4),
    fε::Function=(t->2t), # function to get espilon from h,
    cf::Real=1.0,
    plotting_surface::Bool=false, count::Int64=1
    )

    # Xlim, Ylim, Zlim, Nx, h, dsignes = read_protein("2022"; scaling=2)
    ε = fε(h)*cf

    Nepsl2h=Int(ceil(2*ε/h)) # how many discretization points fall in the tubular neighborhood with these values of ε and h
    h2 = h*h
    h3 = h2*h
    
    X = Array{Float64,1}(-(Nx-1)*h/2:h:(Nx-1)*h/2); 
    Y = Array{Float64,1}(-(Nx-1)*h/2:h:(Nx-1)*h/2); 
    Z = Array{Float64,1}(-(Nx-1)*h/2:h:(Nx-1)*h/2); 
    println("Computational box limits: X=$([X[1];X[end]]), Y=$([Y[1];Y[end]]), Z=$([Z[1];Z[end]])")

    Nvec = [Nx;Nx;Nx]

    println("Current run: h=$h, ε=$ε, Nepsl=$Nepsl2h")
    println("Number of discretization points in each direction:\n$Nx x $Nx x $Nx = $(Nx*Nx*Nx)")

    @time M, indIJK_to_M, indM_to_IJK, ds, normals, Jac_vec, source = genCPM_corr_molecule(dsignes, X, Y, Z, Nx, h, ε)
    # println(source)

    if plotting_surface
        # println(surfTargetXYZ[1])
        # println(surfTargetXYZ)
        xp = [ sour[1] for sour in source ];
        yp = [ sour[2] for sour in source ];
        zp = [ sour[3] for sour in source ];
        figure(50); clf()
        # c=scatter3D( xp, yp, zp, s=1, c=log10.(abs.(β)), marker=".")
        scatter3D( xp, yp, zp, s=1, c="k", marker=".")        
        xlim(Xlim)
        ylim(Xlim)
        zlim(Xlim)

        figure(51); clf()
        # c=scatter3D( xp, yp, zp, s=1, c=log10.(abs.(β)), marker=".")
        sc=scatter3D( xp, yp, zp, s=2, c=Jac_vec, marker=".", cmap=:jet, vmin=-5, vmax = 0.5)
        colorbar(sc)
        xlim(Xlim)
        ylim(Xlim)
        zlim(Xlim)
        title("jacobian")
    end

    salvare_jac_2nd = Array{Float64}(undef,M);
    salvare_ker = Array{Float64}(undef,M);

    @time @threads for m=1:M # for every node inside Tε, compute normal and jacobian
        i,j,k = indM_to_IJK[m,:]; # get ijk-indices from m-index
        # y = [ X[i]; Y[j]; Z[k] ] # node in 3D volume
        salvare_jac_2nd[m] = Jac_vec[m]
        salvare_ker[m] = (abs(dsignes[i,j,k])<ε)*w_ker(dsignes[i,j,k]/ε)/ε # zero outside the tubular neighborhood
    end
    
    salvare_wo = salvare_ker # missing h3
    salvare_w2 = salvare_ker.*salvare_jac_2nd # missing h3
    salvareu = salvare_w2
    
    surf_val = zeros(2)
    surf_val[1] = sum(salvare_wo)*h3
    surf_val[2] = sum(salvare_w2)*h3
    
    println("Val. surface area:\n$(surf_val)")
 
    x=[ surf_val ]
    open(newdir_nrun*"/moleculePB_surf_val_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end

    x=[ h, ε, Nvec, 2*ε/h ]
    open(newdir_nrun*"/moleculePB_h_epsl_data_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end

    println("Run over.")

    return surf_val
end

function PB_gen_shape_molecule_system(N;
    w_ker::Function = Kcos, # averaging kernel
    Xlim, Ylim, Zlim, Nx::Int64=65, h::Real=0.1, dsignes::Array{Float64,3}=zeros(4,4,4),
    fε::Function=(t->2t), # function to get espilon from h,
    cf::Real=1.0,
    epslI::Real=1.0, epslE::Real=1.0, kappa_val::Real=1.0,
    plotting_surface::Bool=false, count::Int64=1
    )

    tmp = readdlm("2022/1ycr.pqr")
    nSpheres = size(tmp,1)-2;
    x0Vec = [ tmp[i,6:8] for i=1:nSpheres]
    RadiiVec = tmp[1:end-2,10];
    qSpheres = tmp[1:end-2,9];

    # Xlim, Ylim, Zlim, Nx, h, dsignes = read_protein("2022"; scaling=2)
    ε = fε(h)*cf

    Nepsl2h=Int(ceil(2*ε/h)) # how many discretization points fall in the tubular neighborhood with these values of ε and h
    h2 = h*h
    h3 = h2*h
    
    X = Array{Float64,1}(-(Nx-1)*h/2:h:(Nx-1)*h/2); 
    Y = Array{Float64,1}(-(Nx-1)*h/2:h:(Nx-1)*h/2); 
    Z = Array{Float64,1}(-(Nx-1)*h/2:h:(Nx-1)*h/2); 
    println("Computational box limits: X=$([X[1];X[end]]), Y=$([Y[1];Y[end]]), Z=$([Z[1];Z[end]])")

    Nvec = [Nx;Nx;Nx]

    epsl_ratio = epslE/epslI;

    println("Current run: h=$h, ε=$ε, Nepsl=$Nepsl2h")
    println("Number of discretization points in each direction:\n$Nx x $Nx x $Nx = $(Nx*Nx*Nx)")

    # @time M, indIJK_to_M, indM_to_IJK, ds, normals, Jac_vec, source = genCPM_corr_molecule(dsignes, X, Y, Z, Nx, h, ε)
    @time M, indIJK_to_M, indM_to_IJK, iv1, w_K11_single, w_K22_single, w_K21_single, normal, ds, Jac_vec, source  = genCPM_corr_PB_system_molecule(dsignes, X, Y, Z, Nx, ε, h; epsl_ratio=epsl_ratio, kappa_val=kappa_val)
    # println(source)

    if plotting_surface
        # println(surfTargetXYZ[1])
        # println(surfTargetXYZ)
        xp = [ sour[1] for sour in source ];
        yp = [ sour[2] for sour in source ];
        zp = [ sour[3] for sour in source ];
        figure(50); clf()
        # c=scatter3D( xp, yp, zp, s=1, c=log10.(abs.(β)), marker=".")
        scatter3D( xp, yp, zp, s=1, c="k", marker=".")        
        xlim(Xlim)
        ylim(Xlim)
        zlim(Xlim)

        figure(51); clf()
        # c=scatter3D( xp, yp, zp, s=1, c=log10.(abs.(β)), marker=".")
        sc=scatter3D( xp, yp, zp, s=2, c=Jac_vec, marker=".", cmap=:jet, vmin=-5, vmax = 0.5)
        colorbar(sc)
        xlim(Xlim)
        ylim(Xlim)
        zlim(Xlim)
        title("jacobian")
    end

    salvare_jac_2nd = Array{Float64}(undef,M);
    salvare_ker = Array{Float64}(undef,M);

    @time @threads for m=1:M # for every node inside Tε, compute normal and jacobian
        i,j,k = indM_to_IJK[m,:]; # get ijk-indices from m-index
        # y = [ X[i]; Y[j]; Z[k] ] # node in 3D volume
        salvare_jac_2nd[m] = Jac_vec[m]
        salvare_ker[m] = (abs(dsignes[i,j,k])<ε)*w_ker(dsignes[i,j,k]/ε)/ε # zero outside the tubular neighborhood
    end
    
    salvare_wo = salvare_ker # missing h3
    salvare_w2 = salvare_ker.*salvare_jac_2nd # missing h3
    salvareu = salvare_w2
    
    surf_val = zeros(2)
    surf_val[1] = sum(salvare_wo)*h3
    surf_val[2] = sum(salvare_w2)*h3
    
    println("Val. surface area:\n$(surf_val)")
 
    x=[ surf_val ]
    open(newdir_nrun*"/moleculePB_surf_val_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end

    # solving the PB problem
    begin # defining the functions and linear maps for solving system
        K11_IBIM_f = (α -> K11_PB_IBIM_target(α; source=source, normal=normal, salvare=salvare_wo, targets=source, tau=tau, kappa_val=kappa_val, theta_val=epsl_ratio)*h3)
        K22_IBIM_f = (α -> K22_PB_IBIM_target(α; source=source, salvare=salvare_wo, targets=source, targetnormals=normal, tau=tau, kappa_val=kappa_val, theta_val=1/epsl_ratio)*h3)
        K21_IBIM_f = (α -> K21_PB_IBIM_target(α; source=source, normal=normal, salvare=salvare_wo, targets=source, targetnormals = normal, tau=tau, kappa_val=kappa_val)*h3)
        K12_IBIM_f = (α -> K12_PB_IBIM_target(α; source=source, salvare=salvare_wo, targets=source, tau=tau, kappa_val=kappa_val)*h3)

        K11_corr_f = (α -> K11_PB_Q1corr_target(α; source=source, normal=normal, salvare=salvareu, targets=source, kappa_val=kappa_val, theta_val=epsl_ratio, h=h, iv1=iv1, wv1=w_K11_single) )
        K22_corr_f = (α -> K22_PB_Q1corr_target(α; source=source, targetnormals=normal, salvare=salvareu, targets=source, kappa_val=kappa_val, theta_val=1/epsl_ratio, h=h, iv1=iv1, wv1=w_K22_single) )
        K21_corr_f = (α -> K21_PB_Q1corr_target(α; source=source, normal=normal,  targetnormals=normal, salvare=salvareu, targets=source, kappa_val=kappa_val, h=h, iv1=iv1, wv1=w_K21_single) )
        K12_corr_f = (α -> K12_PB_Q1corr_target(α; source=source, salvare=salvareu, targets=source, kappa_val=kappa_val, h=h, iv1=iv1) )

        function eq_IBIM(α)
            α1 = α[1:M]
            α2 = α[M+1:2M]
            return [λ1*α1 .+ K11_IBIM_f(α1) .- K12_IBIM_f(α2);λ2*α2 .+ K21_IBIM_f(α1) .- K22_IBIM_f(α2)]
        end

        function eq_corr(α)
            α1 = α[1:M]
            α2 = α[M+1:2M]
            return [λ1*α1 .+ K11_corr_f(α1) .- K12_corr_f(α2);λ2*α2 .+ K21_corr_f(α1) .- K22_corr_f(α2)]
        end
    end
    begin # creating the linear maps
        IBIM_L = LinearMap(eq_IBIM, 2M)
        corr_L = LinearMap(eq_corr, 2M)
    end
    begin # creating the RHS vectors
        g1 = zeros(M); g2 = zeros(M);
        for i=1:nSpheres
            g1 .+= qSpheres[i]*[ Gk_PB(source[m], x0Vec[i], 0) for m=1:M]/epslI
            g2 .+= qSpheres[i]*[ dGkdnx(source[m], x0Vec[i], normal[m], 0) for m=1:M]/epslI
        end
    end

    begin # solving the system with GMRES
        PsiAll_corr = gmres(corr_L, [g1;g2]; verbose=true, abstol=2e-5, reltol=2e-5)
        psi_corr = PsiAll_corr[1:M]
        psin_corr = PsiAll_corr[M+1:2M]
        x = PsiAll_corr
        open(newdir_nrun*"/plotting_corr_N_$(N)_$(nrun[1])_$(count).potent","w") do io
            writedlm(io,x)
        end
        
        PsiAll_IBIM = gmres(IBIM_L, [g1;g2]; verbose=true, abstol=2e-5, reltol=2e-5)
        psi_IBIM = PsiAll_IBIM[1:M]
        psin_IBIM = PsiAll_IBIM[M+1:2M]
        x = PsiAll_IBIM
        open(newdir_nrun*"/plotting_IBIM_N_$(N)_$(nrun[1])_$(count).potent","w") do io
            writedlm(io,x)
        end
    end
    begin # polarization energy computation 
        psi_rxn_IBIM(z) = h3*sum( salvareu .* ([ G0Gk(z,source[m],kappa_val)*psin_IBIM[m] for m=1:M] .- [dGny_diff(z,source[m],normal[m],kappa_val,epsl_ratio)*psi_IBIM[m] for m=1:M] ) )
        psi_rxn_corr(z) = h3*sum( salvareu .* ([ G0Gk(z,source[m],kappa_val)*psin_corr[m] for m=1:M] .- [dGny_diff(z,source[m],normal[m],kappa_val,epsl_ratio)*psi_corr[m] for m=1:M] ) )

        Gpol_IBIM = 0; Gpol_corr = 0;
        for i=1:nSpheres
            Gpol_IBIM += qSpheres[i]*psi_rxn_IBIM(x0Vec[i])
            Gpol_corr += qSpheres[i]*psi_rxn_corr(x0Vec[i])
        end
        Gpol_IBIM *= 0.5; Gpol_corr *= 0.5;

        x = [Gpol_IBIM; Gpol_corr]
        open(newdir_nrun*"/Gpol_values_N_$(N)_$(nrun[1])_$(count).txt","w") do io
            writedlm(io,x)
        end
    end

    x=[ h, ε, Nvec, 2*ε/h ]
    open(newdir_nrun*"/moleculePB_h_epsl_data_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end

    println("Run over.")

    return Gpol_IBIM, Gpol_corr, surf_val
end


function PB_gen_shape_system(h, N;
    w_ker::Function = Kcos, # averaging kernel
    surfTargetXYZ::Array{Array{Float64,1},1}=[[1;0.;0]], # surface targets: identified via 3D point, then projected
    outflag::Bool=true,  # interior or exterior problem: true for exterior
    shift::Array{Float64,1}=[0.;0.;0.], # shift of grid w.r.t. surface center (if any)
    epslI::Real=1.0, epslE::Real=1.0, kappa_val::Real=1.0,
    fε::Function=(t->2t), # function to get espilon from h
    plotting_surface::Bool=false, count::Int64=1,
    Zlim1::Real=-1,  Zlim2::Real=1, zcenter::Array{Float64,1}=[0.;0.;0.],
    hflag::Bool=false, Jflag::Bool=false
    )

    # # single sphere
    # x0Vec = [zeros(3)]
    # RadiiVec = [0.5]
    # R = RadiiVec[1]
    # qSpheres = [1.0]
    # nSpheres = length(RadiiVec)

    # # single torus
    x0Vec = [zeros(3)]
    RadiiVec = [0.5; 1.0]
    R = RadiiVec[1]
    qSpheres = [1.0]
    nSpheres = length(RadiiVec)-1

    # # two spheres conjoined
    # RadiiVec = [0.6; 0.3]
    # x0Vec = [[-RadiiVec[1]+0.2;0;0],[RadiiVec[2]-0.1;0;0]]
    # qSpheres = [1.0;1.5]
    # nSpheres = length(RadiiVec);
    # R3 = 0.2
    
    # R = 1.0
    # R = 0.5
    # R1 = 0.7; R2 = 0.2; 
    # Rcil = 0.4; Zcil = 0.6;
    # ZCIL = [0;0;Zcil]

    # pl=4;
    # h = 1.1*(2R+0.2)/(N)
    # h = (Zlim2-Zlim1)/(N-1)
    if hflag
        N = Int(ceil((Zlim2-Zlim1)/h))
    else
        h = (Zlim2-Zlim1)/(N)
    end
    println("h = $h, N = $N")
    ε = fε(h)
    
    # println(1.1*(2R+0.2)/(N))
    # Xlim = [(-2.3)*1.1-pl*h; (1.5)*1.1+pl*h]
    # Ylim = [(-1.5)*1.1-pl*h; (1.5)*1.1+pl*h]
    # Zlim = [(-1.5)*1.1-pl*h; (1.5)*1.1+pl*h]
    p_rad = 4;
    Xlim = [Zlim1-ε-p_rad*h; Zlim2+ε+p_rad*h].+zcenter[1]
    Ylim = [Zlim1-ε-p_rad*h; Zlim2+ε+p_rad*h].+zcenter[2]
    Zlim = [Zlim1-ε-p_rad*h; Zlim2+ε+p_rad*h].+zcenter[3]
    Xlim = round.(Xlim,digits=6)
    Ylim = round.(Ylim,digits=6)
    Zlim = round.(Zlim,digits=6)
    # Xlim = [(-0.7)*1.1-pl*h; (0.7)*1.1+pl*h]
    # Ylim = [(-0.7)*1.1-pl*h; (0.7)*1.1+pl*h]
    # Zlim = [(-0.7)*1.1-pl*h; (0.7)*1.1+pl*h]
    
    # println("epsl=",ε)
    
    Rx(a) = [1 0 0;0 cos(a) sin(a);0 -sin(a) cos(a)]
    Ry(a) = [cos(a) 0 sin(a);0 1 0;-sin(a) 0 cos(a)]
    Rz(a) = [cos(a) -sin(a) 0;sin(a) cos(a) 0;0 0 1]
    xshift = zeros(3)

    # a,b,c = 0.24402412255508432, 0.7454097947651017, 2.219760487439292
    a,b,c = 1.99487, 2.54097947651017, 4.219760487439292
    # a,b,c = 0., 0., 0.
    
    # a,b,c = 0., 0., 0.
    Qrot = Rz(c)*Ry(b)*Rx(a);
    Rinv = Rx(-a)*Ry(-b)*Rz(-c);
    x0 = xshift; #zeros(3); 

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

    # # single sphere
    # Pgammafun = (z -> Pgammafun_sphere(z,x0Vec[1],RadiiVec[1]) )
    # insidepoint = (z-> insidepoint_sphere(z,x0Vec[1],RadiiVec[1]) )
    # get_jac = (z-> get_jac_sphere(z,x0Vec[1],RadiiVec[1]) )
    # far_insidepoint = (z-> far_insidepoint_sphere(z,x0Vec[1],RadiiVec[1],ε) )
    # far_outsidepoint = (z-> far_outsidepoint_sphere(z,x0Vec[1],RadiiVec[1],ε) )

    # # single torus
    Pgammafun = (z -> Pgamma_torus(z, Qrot, Rinv, RadiiVec[1], RadiiVec[2], xshift) )
    insidepoint = (z-> insidepoint_torus(z, Rinv, RadiiVec[1], RadiiVec[2], xshift) ) 

    # # two spheres conjoined
    # delta = norm(x0Vec[1].-x0Vec[2])-RadiiVec[1]-RadiiVec[2];
    # L1 = RadiiVec[1]+R3
    # L2 = RadiiVec[2]+R3
    # L3 = RadiiVec[1]+RadiiVec[2]+delta
    # alpha = 0.5+(L1^2-L2^2)/(2L3^2)
    # x1a = x0Vec[1]
    # x2a = x0Vec[2]
    # x3a = x0Vec[1].+(alpha)*(x0Vec[2].-x0Vec[1])
    # a2 = sqrt(abs(L1^2-alpha^2*L3^2))
    # b2 = R3;
    # Xa = x0Vec[2].-x0Vec[1]
    # Ba = Xa/norm(Xa); 
    # Aa = [0;1;0.]
    # a3 = dot(Aa,Ba); a4 = norm(cross(Aa,Ba))
    # Ga = [a3 -a4 0;a4 a3 0;0 0 1];
    # u1a = Aa; u2a = Ba-a3*Aa; u2a /=norm(u2a); u3a = cross(Ba,Aa);
    # Fa = ([u1a u2a u3a])
    # Rottot = Fa*Ga*inv(Fa)
    # Rot3 = transpose(Rottot)

    # x4a = Rottot*[0;0;a2].+x3a
    # cone_compA = (1-2*(norm(x2a-x3a)>norm(x1a-x2a)))
    # sintC1 = a2/norm(x4a-x1a); costC1 = norm(x1a-x3a)/norm(x4a-x1a); tant = sintC1/costC1
    # cone_compB = (1-2*(norm(x1a-x3a)>norm(x1a-x2a)))
    # sintC2 = a2/norm(x4a-x2a); costC2 = norm(x2a-x3a)/norm(x4a-x2a); tant = sintC2/costC2
    # # # two touching spheres
    # Pgammafun = ( z -> Pgammafun_2spheres(z, x1a, x2a, Rot3, Rottot, sintC1, costC1, sintC2, costC2, x3a, cone_compA, cone_compB, a2, b2, RadiiVec, x0Vec) )
    # insidepoint = ( z -> insidepoint_2spheres(z, x1a, x2a, Rot3, Rottot, sintC1, costC1, sintC2, costC2, x3a, cone_compA, cone_compB, a2, b2, RadiiVec) )
    

    Nepsl2h=Int(ceil(2*ε/h))+1 # how many discretization points fall in the tubular neighborhood with these values of ε and h
    h2 = h*h
    h3 = h2*h

    # println("Number of surface targets = $(length(surfTargetXYZ))")
    n_surf_trg = length(surfTargetXYZ)
    for i=1:n_surf_trg
        surfTargetXYZ[i] = Pgammafun(surfTargetXYZ[i])
        # target_normal[i] = normalz_sphere(source[m],x0Vec[1],RadiiVec[1])
    end
    
    rshift = xshift+shift;
    # discretizations in each direction, centered in rshift, with step h
    # X = [(rshift[1]:(-h):Xlim[1])[end:-1:1]; (rshift[1]+h):h:Xlim[2]]; Nx = length(X)
    # Y = [(rshift[2]:(-h):Ylim[1])[end:-1:1]; (rshift[2]+h):h:Ylim[2]]; Ny = length(Y)
    # Z = [(rshift[3]:(-h):Zlim[1])[end:-1:1]; (rshift[3]+h):h:Zlim[2]]; Nz = length(Z)
    X = Array{Float64,1}(Xlim[1]:h:Xlim[2]); Nx = length(X)
    Y = Array{Float64,1}(Ylim[1]:h:Ylim[2]); Ny = length(Y)
    Z = Array{Float64,1}(Zlim[1]:h:Zlim[2]); Nz = length(Z)
    println("Computational box limits: X=$([X[1];X[end]]), Y=$([Y[1];Y[end]]), Z=$([Z[1];Z[end]])")

    Nvec = [Nx;Ny;Nz]

    println("Current run: h=$h, ε=$ε, Nepsl=$Nepsl2h")
    println("Number of discr. pts:\n$Nx x $Ny x $Nz = $(Nx*Ny*Nz)")
    
    if (Nx!=Ny || Nx!=Nz || Nz!=Ny)
        @warn("Not a perfect discretization cube")
    end
    
    tau = h # cut-off parameter for function regularization in IBIMs
        
    epsl_ratio = epslE/epslI;
    # @time M, Pg, source, indIJK_to_M, indM_to_IJK, dsignes, ab_single_t, iv1_t, w_K11_single_t, w_K22_single_t, w_K21_single_t, target_normal = genCPM_corr_V2_PB(Pgammafun, insidepoint, far_insidepoint, far_outsidepoint, X, Y, Z, Nvec, ε, h; outflag=outflag, surfTargetXYZ=surfTargetXYZ, epsl_ratio=epsl_ratio, kappa_val=kappa_val)

    @time M, Pg, source, indIJK_to_M, indM_to_IJK, dsignes, iv1, w_K11_single, w_K22_single, w_K21_single, normal, iv1_t, w_K11_single_t, w_K22_single_t, w_K21_single_t, target_normal = genCPM_corr_PB_system(Pgammafun, insidepoint, X, Y, Z, Nvec, ε, h; outflag=outflag, surfTargetXYZ=surfTargetXYZ, epsl_ratio=epsl_ratio, kappa_val=kappa_val)


    x = dsignes
    open(newdir_nrun*"/plotting_N_$(N)_$(nrun[1])_$(count).grid","w") do io
        writedlm(io,x)
    end
    xsource = zeros(M,3)
    for i=1:M
        xsource[i,:] = source[i]
    end
    x = xsource
    open(newdir_nrun*"/plotting_N_$(N)_$(nrun[1])_$(count).node","w") do io
        writedlm(io,x)
    end

    # for i=1:n_surf_trg
    #     # surfTargetXYZ[i] = Pgammafun(surfTargetXYZ[i])
    #     target_normal[i] = normalz_sphere(surfTargetXYZ[i],x0Vec[1],RadiiVec[1])
    # end

    # println("w_K11_single $w_K11_single")
    # println("w_K21_single $w_K21_single")
    # println("w_K22_single $w_K22_single")
    # println("iv1 $iv1")

    dPg = zeros(3,3)
    salvare_jac_2nd = Array{Float64}(undef,M);
    salvare_ker = Array{Float64}(undef,M);
    # normal = Array{Array{Float64,1},1}(undef,M)
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
        # normal[m] = normalz_sphere(source[m],x0Vec[1],RadiiVec[1])
        
        # β[m] = mysmoothf(source[m]) # given density
        β[m] = 1. # given density
    end
    
    salvare_wo = salvare_ker # missing h3
    salvare_w2 = salvare_ker.*salvare_jac_2nd # missing h3
    if Jflag
        salvareu = salvare_ker
    else
        salvareu = salvare_w2
    end
    # println(salvareu)
    
    surf_val = zeros(2)
    surf_val[1] = sum(salvare_wo)*h3
    surf_val[2] = sum(salvare_w2)*h3


    x = surf_val
    open(newdir_nrun*"/surface_values_N_$(N)_$(nrun[1])_$(count).txt","w") do io
        writedlm(io,x)
    end
    
    println("Val. surface area:\n$(surf_val)")

    surf_area = 4*pi*RadiiVec[1]^2

    surf_err_IBIM = abs(surf_area - surf_val[1])/abs(surf_area)
    surf_err_corr = abs(surf_area - surf_val[2])/abs(surf_area)
    println("Area: $surf_area")
    println("Vals: $surf_val")
    println("Surface errors: ")
    println(surf_err_IBIM," , ",surf_err_corr)

    # println(iv1)

    # evaluating the potential for known density β
    # IBIM is with Jacobian=1
    println("Evaluating the potentials in IBIM (4) and CTR (4)")
    @time begin
        K11_IBIM = K11_PB_IBIM_target(β; source=source, normal=normal, salvare=salvareu, targets=surfTargetXYZ, tau=tau, kappa_val=kappa_val, theta_val=epsl_ratio)*h3
        K22_IBIM = K22_PB_IBIM_target(β; source=source, targetnormals=target_normal, salvare=salvareu, targets=surfTargetXYZ, tau=tau, kappa_val=kappa_val, theta_val=1/epsl_ratio)*h3
        K21_IBIM = K21_PB_IBIM_target(β; source=source, normal=normal, targetnormals=target_normal, salvare=salvareu, targets=surfTargetXYZ, tau=tau, kappa_val=kappa_val)*h3
        K12_IBIM = K12_PB_IBIM_target(β; source=source, salvare=salvareu, targets=surfTargetXYZ, tau=tau, kappa_val=kappa_val)*h3

        K11_corr = K11_PB_Q1corr_target(β; source=source, normal=normal, salvare=salvareu, targets=surfTargetXYZ, kappa_val=kappa_val, theta_val=epsl_ratio, h=h, iv1=iv1_t, wv1=w_K11_single_t)
        K22_corr = K22_PB_Q1corr_target(β; source=source, targetnormals=target_normal, salvare=salvareu, targets=surfTargetXYZ, kappa_val=kappa_val, theta_val=1/epsl_ratio, h=h, iv1=iv1_t, wv1=w_K22_single_t)
        K21_corr = K21_PB_Q1corr_target(β; source=source, normal=normal,  targetnormals=target_normal, salvare=salvareu, targets=surfTargetXYZ, kappa_val=kappa_val, h=h, iv1=iv1_t, wv1=w_K21_single_t)
        K12_corr = K12_PB_Q1corr_target(β; source=source, salvare=salvareu, targets=surfTargetXYZ, kappa_val=kappa_val, h=h, iv1=iv1_t)

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

    begin # analytic values for single sphere
        λ1 = 0.5*(1+epsl_ratio); λ2 = 0.5*(1+1/epsl_ratio);

        psi_an = ones(M)*qSpheres[1]/(4*pi*epslE*(1+kappa_val*RadiiVec[1])*RadiiVec[1])
        psin_an = -ones(M)*qSpheres[1]/(4*pi*epslI*RadiiVec[1]^2)

        Gpol_an = qSpheres[1]^2*(1/(epslE*(1+kappa_val*RadiiVec[1]))-1/epslI)/(8*pi*RadiiVec[1])
        surf_area = 4*pi*RadiiVec[1]^2
    end

    if ~false
        g1 = qSpheres[1]/(4*pi*epslI*RadiiVec[1])
        g2 = -qSpheres[1]/(4*pi*epslI*RadiiVec[1]^2)

        val_abs = zeros(8,n_surf_trg)
        val_abs[1,:] = K11_IBIM
        val_abs[2,:] = K22_IBIM
        val_abs[3,:] = K21_IBIM
        val_abs[4,:] = K12_IBIM
        val_abs[5,:] = K11_corr
        val_abs[6,:] = K22_corr
        val_abs[7,:] = K21_corr
        val_abs[8,:] = K12_corr

        val_err2 = zeros(4,n_surf_trg)
        val_err2[1,:] = abs.(λ1*psi_an[1] .+ K11_IBIM*psi_an[1] .- K12_IBIM*psin_an[1] .- g1)./abs(g1)
        val_err2[2,:] = abs.(λ2*psin_an[1] .+ K21_IBIM*psi_an[1] .- K22_IBIM*psin_an[1] .- g2)./abs(g2)
        val_err2[3,:] = abs.(λ1*psi_an[1] .+ K11_corr*psi_an[1] .- K12_corr*psin_an[1] .- g1)./abs(g1)
        val_err2[4,:] = abs.(λ2*psin_an[1] .+ K21_corr*psi_an[1] .- K22_corr*psin_an[1] .- g2)./abs(g2)


        println("psi values: $(psi_an[1]), $(psin_an[1])")
        println("g1,g2 = $g1, $g2")
        avgK11_IBIM = mean(K11_IBIM)
        avgK11_corr = mean(K11_corr)
        avgK22_IBIM = mean(K22_IBIM)
        avgK22_corr = mean(K22_corr)
        avgK21_IBIM = mean(K21_IBIM)
        avgK21_corr = mean(K21_corr)
        avgK12_IBIM = mean(K12_IBIM)
        avgK12_corr = mean(K12_corr)
        println("Avg K11, IBIM/corr: $avgK11_IBIM, $avgK11_corr")
        println("Avg K22, IBIM/corr: $avgK22_IBIM, $avgK22_corr")
        println("Avg K21, IBIM/corr: $avgK21_IBIM, $avgK21_corr")
        println("Avg K12, IBIM/corr: $avgK12_IBIM, $avgK12_corr")

        mval_err2 = mean(val_err2, dims=2)
        println(mval_err2)

        # # println("Node potentials:")
        # @time K11_corrb = K11_PB_Q1corr_target(β; source=source, normal=normal, salvare=salvareu, targets=source, kappa_val=kappa_val, theta_val=epsl_ratio, h=h, iv1=iv1, wv1=w_K11_single)
        # # println(tmp1)
        # @time K22_corrb = K22_PB_Q1corr_target(β; source=source, targetnormals=normal, salvare=salvareu, targets=source, kappa_val=kappa_val, theta_val=1/epsl_ratio, h=h, iv1=iv1, wv1=w_K22_single)
        # # println(tmp1)
        # @time K21_corrb = K21_PB_Q1corr_target(β; source=source, normal=normal,  targetnormals=normal, salvare=salvareu, targets=source, kappa_val=kappa_val, h=h, iv1=iv1, wv1=w_K21_single)
        # # println(tmp1)
        # @time K12_corrb = K12_PB_Q1corr_target(β; source=source, salvare=salvareu, targets=source, kappa_val=kappa_val, h=h, iv1=iv1)
        # # println(tmp1)
        # # println(K22_corr)
        # # println(K21_corr)
        # # println(K12_corr)
        # avgK11_corrb = mean(K11_corrb)
        # avgK22_corrb = mean(K22_corrb)
        # avgK21_corrb = mean(K21_corrb)
        # avgK12_corrb = mean(K12_corrb)
        # println("Avg K11, corrb: $avgK11_corrb")
        # println("Avg K22, corrb: $avgK22_corrb")
        # println("Avg K21, corrb: $avgK21_corrb")
        # println("Avg K12, corrb: $avgK12_corrb")


        # val_err = zeros(4,M)
        # val_err[1,:] = abs.(λ1*psi_an[1] .+ K11_IBIM*psi_an[1] .- K12_IBIM*psin_an[1] .- g1)./abs(g1)
        # val_err[2,:] = abs.(λ2*psin_an[1] .+ K21_IBIM*psi_an[1] .- K22_IBIM*psin_an[1] .- g2)./abs(g2)
        # val_err[3,:] = abs.(λ1*psi_an[1] .+ K11_corrb*psi_an[1] .- K12_corrb*psin_an[1] .- g1)./abs(g1)
        # val_err[4,:] = abs.(λ2*psin_an[1] .+ K21_corrb*psi_an[1] .- K22_corrb*psin_an[1] .- g2)./abs(g2)

        mval_abs = mean(val_abs, dims=2)
        # mval_err = mean(val_err, dims=2)

        # println("Errors for sphere:")
        # println(mval_err)
        println("Average of potentials over targets:\n",mval_abs[1:8])
        # println(mval_err)
    end

    # solving the PB problem
    begin # defining the functions and linear maps for solving system
        K11_IBIM_f = (α -> K11_PB_IBIM_target(α; source=source, normal=normal, salvare=salvare_wo, targets=source, tau=tau, kappa_val=kappa_val, theta_val=epsl_ratio)*h3)
        K22_IBIM_f = (α -> K22_PB_IBIM_target(α; source=source, salvare=salvare_wo, targets=source, targetnormals=normal, tau=tau, kappa_val=kappa_val, theta_val=1/epsl_ratio)*h3)
        K21_IBIM_f = (α -> K21_PB_IBIM_target(α; source=source, normal=normal, salvare=salvare_wo, targets=source, targetnormals = normal, tau=tau, kappa_val=kappa_val)*h3)
        K12_IBIM_f = (α -> K12_PB_IBIM_target(α; source=source, salvare=salvare_wo, targets=source, tau=tau, kappa_val=kappa_val)*h3)

        K11_corr_f = (α -> K11_PB_Q1corr_target(α; source=source, normal=normal, salvare=salvareu, targets=source, kappa_val=kappa_val, theta_val=epsl_ratio, h=h, iv1=iv1, wv1=w_K11_single) )
        K22_corr_f = (α -> K22_PB_Q1corr_target(α; source=source, targetnormals=normal, salvare=salvareu, targets=source, kappa_val=kappa_val, theta_val=1/epsl_ratio, h=h, iv1=iv1, wv1=w_K22_single) )
        K21_corr_f = (α -> K21_PB_Q1corr_target(α; source=source, normal=normal,  targetnormals=normal, salvare=salvareu, targets=source, kappa_val=kappa_val, h=h, iv1=iv1, wv1=w_K21_single) )
        K12_corr_f = (α -> K12_PB_Q1corr_target(α; source=source, salvare=salvareu, targets=source, kappa_val=kappa_val, h=h, iv1=iv1) )

        function eq_IBIM(α)
            α1 = α[1:M]
            α2 = α[M+1:2M]
            return [λ1*α1 .+ K11_IBIM_f(α1) .- K12_IBIM_f(α2);λ2*α2 .+ K21_IBIM_f(α1) .- K22_IBIM_f(α2)]
        end

        function eq_corr(α)
            α1 = α[1:M]
            α2 = α[M+1:2M]
            return [λ1*α1 .+ K11_corr_f(α1) .- K12_corr_f(α2);λ2*α2 .+ K21_corr_f(α1) .- K22_corr_f(α2)]
        end
    end
    begin # creating the linear maps
        IBIM_L = LinearMap(eq_IBIM, 2M)
        corr_L = LinearMap(eq_corr, 2M)
    end
    begin # creating the RHS vectors
        g1 = zeros(M); g2 = zeros(M);
        for i=1:nSpheres
            g1 .+= qSpheres[i]*[ Gk_PB(source[m], x0Vec[i], 0) for m=1:M]/epslI
            g2 .+= qSpheres[i]*[ dGkdnx(source[m], x0Vec[i], normal[m], 0) for m=1:M]/epslI
        end
    end

    begin # solving the system with GMRES
        PsiAll_corr = gmres(corr_L, [g1;g2]; verbose=true, abstol=2e-5, reltol=2e-5)
        psi_corr = PsiAll_corr[1:M]
        psin_corr = PsiAll_corr[M+1:2M]
        x = PsiAll_corr
        open(newdir_nrun*"/plotting_corr_N_$(N)_$(nrun[1])_$(count).potent","w") do io
            writedlm(io,x)
        end
        
        PsiAll_IBIM = gmres(IBIM_L, [g1;g2]; verbose=true, abstol=2e-5, reltol=2e-5)
        psi_IBIM = PsiAll_IBIM[1:M]
        psin_IBIM = PsiAll_IBIM[M+1:2M]
        x = PsiAll_IBIM
        open(newdir_nrun*"/plotting_IBIM_N_$(N)_$(nrun[1])_$(count).potent","w") do io
            writedlm(io,x)
        end
    end
    begin # polarization energy computation 
        psi_rxn_IBIM(z) = h3*sum( salvareu .* ([ G0Gk(z,source[m],kappa_val)*psin_IBIM[m] for m=1:M] .- [dGny_diff(z,source[m],normal[m],kappa_val,epsl_ratio)*psi_IBIM[m] for m=1:M] ) )
        psi_rxn_corr(z) = h3*sum( salvareu .* ([ G0Gk(z,source[m],kappa_val)*psin_corr[m] for m=1:M] .- [dGny_diff(z,source[m],normal[m],kappa_val,epsl_ratio)*psi_corr[m] for m=1:M] ) )

        Gpol_IBIM = 0; Gpol_corr = 0;
        for i=1:nSpheres
            Gpol_IBIM += qSpheres[i]*psi_rxn_IBIM(x0Vec[i])
            Gpol_corr += qSpheres[i]*psi_rxn_corr(x0Vec[i])
        end
        Gpol_IBIM *= 0.5; Gpol_corr *= 0.5;

        x = [Gpol_IBIM; Gpol_corr]
        open(newdir_nrun*"/Gpol_values_N_$(N)_$(nrun[1])_$(count).txt","w") do io
            writedlm(io,x)
        end
    end

    begin # errors from analytic values
        sol_err_IBIM = sqrt( (h3*sum( salvareu .* (abs.(psi_an.-psi_IBIM).^2 .+ abs.(psin_an.-psin_IBIM).^2) )) / (h3*sum( salvareu .* (abs.(psi_an).^2 .+ abs.(psin_an).^2) )) )
        sol_err_corr = sqrt( (h3*sum( salvareu .* (abs.(psi_an.-psi_corr).^2 .+ abs.(psin_an.-psin_corr).^2) )) / (h3*sum( salvareu .* (abs.(psi_an).^2 .+ abs.(psin_an).^2) )) )
        
        surf_err_IBIM = abs(surf_area - surf_val[1])/abs(surf_area)
        surf_err_corr = abs(surf_area - surf_val[2])/abs(surf_area)

        en_err_IBIM = abs(Gpol_an-Gpol_IBIM)/abs(Gpol_an)
        en_err_corr = abs(Gpol_an-Gpol_corr)/abs(Gpol_an)

        errvec = [sol_err_IBIM surf_err_IBIM en_err_IBIM;sol_err_corr surf_err_corr en_err_corr]

        println("Error vector (IBIM/corr):")
        println("Solution error: $(errvec[:,1])")
        println("Surface  error: $(errvec[:,2])")
        println("Energy   error: $(errvec[:,3])")

        println("Gpol (IBIM/corr):")
        println(Gpol_IBIM)
        println(Gpol_corr)
        println("Surf-area (IBIM/corr):")
        println(surf_val[1])
        println(surf_val[2])
        # println()
    end
 
    # x=[ R, R1, R2]
    # open(newdir_nrun*"/surfacePB_R_R1_R2_Rcil_data_$(nrun[1])_$(count).dat","w") do io
    #     writedlm(io,x)
    # end

    x=[ RadiiVec, nSpheres]
    open(newdir_nrun*"/surfacePB_Radii_nSpheres_data_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end
    x=[ x0Vec, qSpheres]
    open(newdir_nrun*"/surfacePB_x0Vec_qSpheres_data_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end
    
    
    x=[ a, b, c, xshift ]
    open(newdir_nrun*"/surfacePB_abc_xshift_data_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end
    x=[ kappa_val, epslE, epslI]
    open(newdir_nrun*"/surfacePB_kappa_epslEI_data_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end

    # errors from solving the system
    x=[ errvec ]
    open(newdir_nrun*"/surfacePB_val_errvec_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end

    # potentials evaluated in targets
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

    return h, val_abs, mval_abs, surf_val, errvec
end

function PB_gen_shape_system_DEBUG(N;
    w_ker::Function = Kinf, # averaging kernel
    surfTargetXYZ::Array{Array{Float64,1},1}=[[1;0.;0]], # surface targets: identified via 3D point, then projected
    outflag::Bool=true,  # interior or exterior problem: true for exterior
    shift::Array{Float64,1}=[0.;0.;0.], # shift of grid w.r.t. surface center (if any)
    epslI::Real=1.0, epslE::Real=1.0, kappa_val::Real=1.0,
    fε::Function=(t->2t), # function to get espilon from h
    plotting_surface::Bool=false, count::Int64=1
    )

    RadiiVec = [0.5]
    x0Vec = [zeros(3)]
    # RadiiVec = [1.0; 0.7]
    # x0Vec = [[-RadiiVec[1];0;0],[RadiiVec[2];0;0]]
    qSpheres = [1.0;1.5]
    nSpheres = length(RadiiVec);
    
    # R = 1.0
    R = 0.6
    R1 = 0.7; R2 = 0.2; 
    Rcil = 0.4; Zcil = 0.6;
    ZCIL = [0;0;Zcil]

    pl=4;
    h = 1.1*(2R+0.2)/(N)
    println(1.1*(2R+0.2)/(N))
    # Xlim = [(-2.3)*1.1-pl*h; (1.5)*1.1+pl*h]
    # Ylim = [(-1.5)*1.1-pl*h; (1.5)*1.1+pl*h]
    # Zlim = [(-1.5)*1.1-pl*h; (1.5)*1.1+pl*h]
    Xlim = [(-0.6)*1.1-pl*h; (0.6)*1.1+pl*h]
    Ylim = [(-0.6)*1.1-pl*h; (0.6)*1.1+pl*h]
    Zlim = [(-0.6)*1.1-pl*h; (0.6)*1.1+pl*h]
    
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

    # single sphere
    Pgammafun = (z -> Pgammafun_sphere(z,x0Vec[1],RadiiVec[1]) )
    insidepoint = (z-> insidepoint_sphere(z,x0Vec[1],RadiiVec[1]) )
    get_jac = (z-> get_jac_sphere(z,x0Vec[1],RadiiVec[1]) )
    far_insidepoint = (z-> far_insidepoint_sphere(z,x0Vec[1],RadiiVec[1],ε) )
    far_outsidepoint = (z-> far_outsidepoint_sphere(z,x0Vec[1],RadiiVec[1],ε) )
    
    # # two touching spheres
    # Pgammafun = Pgammafun_spheres
    # insidepoint = insidepoint_spheres
    # get_jac = get_jac_spheres
    # far_insidepoint = far_insidepoint_spheres
    # far_outsidepoint = far_outsidepoint_spheres
    

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
    # @time M, Pg, source, indIJK_to_M, indM_to_IJK, dsignes, ab_single_t, iv1_t, w_K11_single_t, w_K22_single_t, w_K21_single_t, target_normal = genCPM_corr_V2_PB(Pgammafun, insidepoint, far_insidepoint, far_outsidepoint, X, Y, Z, Nvec, ε, h; outflag=outflag, surfTargetXYZ=surfTargetXYZ, epsl_ratio=epsl_ratio, kappa_val=kappa_val)

    @time M, Pg, source, indIJK_to_M, indM_to_IJK, dsignes, iv1, w_K11_single, w_K22_single, w_K21_single, normal, iv1_t, w_K11_single_t, w_K22_single_t, w_K21_single_t, target_normal = genCPM_corr_PB_system(Pgammafun, insidepoint, far_insidepoint, far_outsidepoint, X, Y, Z, Nvec, ε, h; outflag=outflag, surfTargetXYZ=surfTargetXYZ, epsl_ratio=epsl_ratio, kappa_val=kappa_val)

    # println("w_K11_single $w_K11_single")
    # println("w_K21_single $w_K21_single")
    # println("w_K22_single $w_K22_single")
    # println("iv1 $iv1")

    dPg = zeros(3,3)
    salvare_jac_2nd = Array{Float64}(undef,M);
    salvare_ker = Array{Float64}(undef,M);
    # normal = Array{Array{Float64,1},1}(undef,M)
    β = Array{Float64}(undef,M)

    mysmoothf(Z) = 1.12678+sin(Z[1])*Z[2]+cos(Z[3])*0.2376*(Z[1]-0.3176)
    @time for m=1:M # for every node inside Tε, compute normal and jacobian
        i,j,k = indM_to_IJK[m,:]; # get ijk-indices from m-index
        # y = [ X[i]; Y[j]; Z[k] ] # node in 3D volume
        
        # computing the normal with potentially one-sided 2nd order scheme
        # iind, i1 = get_indices_b(dsignes[i-2:i+2,j,k])
        # jind, j1 = get_indices_b(dsignes[i,j-2:j+2,k])
        # kind, k1 = get_indices_b(dsignes[i,j,k-2:k+2])
        # normal[m] = [dot(dsignes[iind.+i,j,k],i1); dot(dsignes[i,jind.+j,k],j1);dot(dsignes[i,j,kind.+k],k1)]/h
        # normal[m] /= norm(normal[m])

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

    println(salvareu)
    
    surf_val = zeros(2)
    surf_val[1] = sum(salvare_wo)*h3
    surf_val[2] = sum(salvare_w2)*h3
    
    println("Val. surface area:\n$(surf_val)")

    return h, iv1, w_K11_single, w_K22_single, w_K21_single, normal, salvareu, source, normal

    # evaluating the potential for known density β
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
        K12_corr = K12_PB_Q1corr_target(β; source=source, salvare=salvareu, targets=surfTargetXYZ, kappa_val=kappa_val, h=h, iv1=iv1_t)

        # println("Node potentials:")
        # @time tmp1 = K11_PB_Q1corr_target(β; source=source, normal=normal, salvare=salvareu, targets=source, kappa_val=kappa_val, theta_val=epsl_ratio, h=h, iv1=iv1, wv1=w_K11_single)
        # println(tmp1)
        # @time tmp1 = K22_PB_Q1corr_target(β; source=source, targetnormals=normal, salvare=salvareu, targets=source, kappa_val=kappa_val, theta_val=1/epsl_ratio, h=h, iv1=iv1, wv1=w_K22_single)
        # println(tmp1)
        # @time tmp1 = K21_PB_Q1corr_target(β; source=source, normal=normal,  targetnormals=normal, salvare=salvareu, targets=source, kappa_val=kappa_val, h=h, iv1=iv1, wv1=w_K21_single)
        # println(tmp1)
        # @time tmp1 = K12_PB_Q1corr_target(β; source=source, salvare=salvareu, targets=source, kappa_val=kappa_val, h=h, iv1=iv1)
        # println(tmp1)
        # println(K22_corr)
        # println(K21_corr)
        # println(K12_corr)
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

    # solving the PB problem
    begin # defining the functions and linear maps for solving system
        λ1 = 0.5*(1+epsl_ratio); λ2 = 0.5*(1+1/epsl_ratio);

        K11_IBIM_f = (α -> K11_PB_IBIM_target(α; source=source, normal=normal, salvare=salvare_wo, targets=source, tau=tau, kappa_val=kappa_val, theta_val=epsl_ratio)*h3)
        K22_IBIM_f = (α -> K22_PB_IBIM_target(α; source=source, normal=normal, salvare=salvare_wo, targets=source, targetnormals=normal, tau=tau, kappa_val=kappa_val, theta_val=1/epsl_ratio)*h3)
        K21_IBIM_f = (α -> K21_PB_IBIM_target(α; source=source, normal=normal, salvare=salvare_wo, targets=source, targetnormals = normal, tau=tau, kappa_val=kappa_val)*h3)
        K12_IBIM_f = (α -> K12_PB_IBIM_target(α; source=source, salvare=salvare_wo, targets=source, tau=tau, kappa_val=kappa_val)*h3)

        K11_corr_f = (α -> K11_PB_Q1corr_target(α; source=source, normal=normal, salvare=salvareu, targets=source, kappa_val=kappa_val, theta_val=epsl_ratio, h=h, iv1=iv1, wv1=w_K11_single) )
        K22_corr_f = (α -> K22_PB_Q1corr_target(α; source=source, targetnormals=normal, salvare=salvareu, targets=source, kappa_val=kappa_val, theta_val=1/epsl_ratio, h=h, iv1=iv1, wv1=w_K22_single) )
        K21_corr_f = (α -> K21_PB_Q1corr_target(α; source=source, normal=normal,  targetnormals=normal, salvare=salvareu, targets=source, kappa_val=kappa_val, h=h, iv1=iv1, wv1=w_K21_single) )
        K12_corr_f = (α -> K12_PB_Q1corr_target(α; source=source, salvare=salvareu, targets=source, kappa_val=kappa_val, h=h, iv1=iv1) )

        function eq_IBIM(α)
            α1 = α[1:M]
            α2 = α[M+1:2M]
            return [λ1*α1 .+ K11_IBIM_f(α1) .- K12_IBIM_f(α2);λ2*α2 .+ K21_IBIM_f(α1) .- K22_IBIM_f(α2)]
        end

        function eq_corr(α)
            α1 = α[1:M]
            α2 = α[M+1:2M]
            return [λ1*α1 .+ K11_corr_f(α1) .- K12_corr_f(α2);λ2*α2 .+ K21_corr_f(α1) .- K22_corr_f(α2)]
        end
    end
    begin # creating the linear maps
        IBIM_L = LinearMap(eq_IBIM, 2M)
        corr_L = LinearMap(eq_corr, 2M)
    end
    begin # creating the RHS vectors
        g1 = zeros(M); g2 = zeros(M);
        for i=1:nSpheres
            g1 .+= qSpheres[i]*[ Gk_PB(source[m], RadiiVec[i], 0) for m=1:M]/epslI
            g2 .+= qSpheres[i]*[ dGkdnx(source[m], RadiiVec[i], normal[m], 0) for m=1:M]/epslI
        end
    end

    begin # solving the system with GMRES
        PsiAll_corr = gmres(corr_L, [g1;g2]; verbose=true)
        PsiAll_IBIM = gmres(IBIM_L, [g1;g2]; verbose=true)
        psi_IBIM = PsiAll_IBIM[1:M]
        psin_IBIM = PsiAll_IBIM[M+1:2M]
        psi_corr = PsiAll_corr[1:M]
        psin_corr = PsiAll_corr[M+1:2M]
    end
    begin # polarization energy computation 
        psi_rxn_IBIM(z) = h3*sum( salvareu .* ([ G0Gk(z,source[m],kappa_val)*psin_IBIM[m] for m=1:M] .- [dGny_diff(z,source[m],normal[m],kappa_val,epsl_ratio)*psi_IBIM[m] for m=1:M] ) )
        psi_rxn_corr(z) = h3*sum( salvareu .* ([ G0Gk(z,source[m],kappa_val)*psin_corr[m] for m=1:M] .- [dGny_diff(z,source[m],normal[m],kappa_val,epsl_ratio)*psi_corr[m] for m=1:M] ) )

        Gpol_IBIM = 0; Gpol_corr = 0;
        for i=1:nSpheres
            Gpol_IBIM += qSpheres[i]*psi_rxn_IBIM(x0Vec[i])
            Gpol_corr += qSpheres[i]*psi_rxn_corr(x0Vec[i])
        end
        Gpol_IBIM *= 0.5; Gpol_corr *= 0.5;
    end

    begin # analytic values for single sphere
        psi_an = ones(M)*qSpheres[1]/(4*pi*epslE*(1+kappa_val*RadiiVec[1])*RadiiVec[1])
        psin_an = -ones(M)*qSpheres[1]/(4*pi*epslI*RadiiVec[1]^2)
        Gpol_an = qSpheres[1]^2*(1/(epslE*(1+kappa_val*RadiiVec[1]))-1/epslI)/(8*pi*RadiiVec[1])
        surf_area = 4*pi*Radii[1]^2
    end
    begin # errors from analytic values
        sol_err_IBIM = sqrt( (h3*sum( salvareu .* (abs.(psi_an.-psi_IBIM).^2 .+ abs.(psin_an.-psin_IBIM).^2) )) / (h3*sum( salvareu .* (abs.(psi_an).^2 .+ abs.(psin_an).^2) )) )
        sol_err_corr = sqrt( (h3*sum( salvareu .* (abs.(psi_an.-psi_corr).^2 .+ abs.(psin_an.-psin_corr).^2) )) / (h3*sum( salvareu .* (abs.(psi_an).^2 .+ abs.(psin_an).^2) )) )
        
        surf_err_IBIM = abs(surf_area - h3*sum(salvarewo))/abs(surf_area)
        surf_err_corr = abs(surf_area - h3*sum(salvareu))/abs(surf_area)

        en_err_IBIM = abs(Gpol_an-Gpol_IBIM)/abs(Gpol_an)
        en_err_corr = abs(Gpol_an-Gpol_corr)/abs(Gpol_an)

        errvec = [sol_err_IBIM surf_err_IBIM en_err_IBIM;sol_err_corr surf_err_corr en_err_corr]

        println("Error vector (IBIM/corr):")
        println("Solution error: $(errvec[:,1])")
        println("Surface  error: $(errvec[:,2])")
        println("Energy   error: $(errvec[:,3])")
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

    x=[ RadiiVec, nSpheres]
    open(newdir_nrun*"/surfacePB_Radii_nSpheres_data_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end
    x=[ x0Vec, qSpheres]
    open(newdir_nrun*"/surfacePB_x0Vec_qSpheres_data_$(nrun[1])_$(count).dat","w") do io
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

    # errors from solving the system
    x=[ errvec ]
    open(newdir_nrun*"/surfacePB_val_errvec_$(nrun[1])_$(count).dat","w") do io
        writedlm(io,x)
    end

    # potentials evaluated in targets
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

    return h, val_abs, val_err, mval_abs, mval_err, surf_val, errvec
end