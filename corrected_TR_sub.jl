# Subroutines for support of CTRvsIBIM for Poisson-Boltzmann

# standalone kernels (kappa=0 -> Laplace)
Gk_PB(x,y,kappa) = exp(-kappa*norm(x-y))/(4*pi*norm(x-y))
dGkdnx(x,y,nx,kappa) = -exp(-kappa*norm(x-y))*(1+kappa*norm(x-y))*dot(x-y,nx)/(4*pi*norm(x-y)^3)
dGkdny(x,y,ny,kappa) = exp(-kappa*norm(x-y))*(1+kappa*norm(x-y))*dot(x-y,ny)/(4*pi*norm(x-y)^3)

d2Gkdnxdny(x,y,nx,ny,kappa) = exp(-kappa*norm(x-y))*(1+kappa*norm(x-y))*(dot(nx,ny)-3*dot(x-y,nx)*dot(x-y,ny)/norm(x-y)^2)/(4*pi*norm(x-y)^3) -kappa^2*exp(-kappa*norm(x-y))*dot(x-y,nx)*dot(x-y,ny)/(4*pi*norm(x-y)^3)

# kernel differences as they appear in the formulations

# difference dG0dnx - theta dGkdnx
dGnx_diff(x,y,nx,kappa,theta) = (1-theta*exp(-kappa*norm(x-y))*(1+kappa*norm(x-y)))*(dot(y-x,nx)/norm(x-y)^3)/(4*pi) # notice (y-x) inside the dot, to get the minus sign
# usually theta = epsI/epsE

# difference dG0dny - theta dGkdny
dGny_diff(x,y,ny,kappa,theta) = (1-theta*exp(-kappa*norm(x-y))*(1+kappa*norm(x-y)))*(dot(x-y,ny)/norm(x-y)^3)/(4*pi) # notice (x-y) inside the dot, to retain the plus sign in constrast with the other kernel
# usually theta = epsE/epsI

# difference G0 - Gk
G0Gk(x,y,kappa) = (1-exp(-kappa*norm(x-y)))/(4*pi*norm(x-y))

# difference d2G0dnxdny - d2Gkdnydny, which is integrable
d2G_diff(x,y,nx,ny,kappa) = (dot(nx,ny)-3*dot(x-y,nx)*dot(x-y,ny)/norm(x-y)^2)*(1-(1+kappa*norm(x-y))*exp(-kappa*norm(x-y)))/(4*pi*norm(x-y)^3) - kappa^2*exp(-kappa*norm(x-y))*dot(x-y,nx)*dot(x-y,ny)/(4*pi*norm(x-y)^3)

# IBIM regularizations

K22_PB(x,y,nx,kappa,theta,tau) = (norm(x-y)>=tau)*dGnx_diff(x,y,nx,kappa,theta)
K21_PB(x,y,nx,ny,kappa,tau) = (norm(x-y)>=tau)*d2G_diff(x,y,nx,ny,kappa)

K11_PB(x,y,ny,kappa,theta,tau) = (norm(x-y)>=tau)*dGny_diff(x,y,ny,kappa,theta)
K12_PB(x,y,kappa,tau) = (norm(x-y)>=tau)*G0Gk(x,y,kappa) + (norm(x-y)<tau)*(exp(-kappa*tau)-1+kappa*tau)/(2*pi*kappa*tau^2)

# solution ψ* to the Single Ion Problem (sphere of radius r):
ψ_SIP(q,r,kappa,epsE) = q/(4*pi*epsE*(1+kappa*r)*r)
# q charge of the ion, epsE dielectric constant of the solvent
# its normal derivative dψ*/dn
ψn_SIP(q,r,epsI) = -q/(4*pi*epsI*r^2)
# epsI dielectric constant of the macromolecule

#################### 
# expansions

# two terms in the difference of the hypersingular kernels
F1_PB(x,y,nx,ny,k) = (dot(nx,ny)-3( dot(x.-y,nx)*dot(x.-y,ny)/(norm(x.-y)^2) ))*(1-exp(-k*norm(x.-y))-k*norm(x.-y)*exp(-k*norm(x.-y)))/(norm(x.-y)^3);
F2_PB(x,y,nx,ny,k) = exp(-k*norm(x.-y))*dot(x.-y,nx)*dot(x.-y,ny)/(norm(x.-y)^3);

Bfun_PB(x,v) = 0.5*(v[1]*x[1]^3/3+v[4]*x[2]^3/3+v[2]*x[1]^2*x[2]+v[3]*x[1]*x[2]^2)
Cfun_PB(x,v) = 0.5*[v[1]*x[1]^2+2*v[2]*x[1]*x[2]+v[3]*x[2]^2; v[4]*x[2]^2+2*v[3]*x[1]*x[2]+v[2]*x[1]^2]

#  x_p=h(Ay,dvec^Ty+η)=χ0(y)y + χ1(y)y^2 + ...
χ0_PB(x,D0,Amat) = D0*Amat*x # vector
χ1_PB(x,D0,Amat,dvec,Mmat,z,v) = D0*D0*Mmat*Amat*x*dot(dvec,x)+z*D0*Cfun_PB(D0*Amat*x,v) # vector

#  f(x_p)=ξ0(y)y^2 + ξ1(y)y^3 + ...
ψ0_PB(x,D0,Amat) = norm(D0*Amat*x); # scalar
ψ0_PB_a(x,D0,Amat) = ψ0_PB([cos(x);sin(x)],D0,Amat); # scalar
ψ1_PB(x,D0,Amat,dvec,Mmat,z,v) = dot(D0*Amat*x,χ1_PB(x,D0,Amat,dvec,Mmat,z,v))/ψ0_PB(x,D0,Amat) # scalar

# nx . ny = 1+a2(y) y^2 + a3(y)y^3 + ...
a2_PB(x,D0,Amat,Mmat) = -0.5*norm(Mmat*D0*Amat*x)^2
a3_PB(x,D0,Amat,dvec,Mmat,η,v) = -dot(Mmat*D0*Amat*x, dot(dvec,x)*Mmat*D0*D0*Mmat*Amat*x+([1 0;0 1]+η*Mmat*D0)*Cfun_PB(D0*Amat*x,v));

# ∇f(x_p) = c1(y)y + c2(y)y^2 + ...
c1_PB(x,D0,Amat,Mmat) = Mmat*D0*Amat*x; # vector
c2_PB(x,D0,Amat,dvec,Mmat,η,v) = dot(dvec,x)*Mmat*D0*D0*Mmat*Amat*x+([1 0;0 1]+η*Mmat*D0)*Cfun_PB(D0*Amat*x,v); # vector

ξ0_PB(x,D0,Amat,Mmat) = 0.5*dot(D0*Amat*x,Mmat*D0*Amat*x) # scalar
ξ1_PB(x,D0,Amat,dvec,Mmat,z,v) = dot(x,dvec)*dot(Amat*x, Mmat*D0*D0*Mmat*D0*Amat*x)+0.5*z*dot(D0*Cfun_PB(D0*Amat*x,v) , Mmat*D0*Amat*x) + 0.5*z*dot(D0*Amat*x, Mmat*D0*Cfun_PB(D0*Amat*x,v))+Bfun_PB(D0*Amat*x,v) # scalar

# cosθ_x = b1(y)y + b2(y)y^2 + ...
b1_PB(x,D0,Amat,Mmat) = -ξ0_PB(x,D0,Amat,Mmat)/ψ0_PB(x,D0,Amat) # scalar
b2_PB(x,D0,Amat,dvec,Mmat,z,v) = (ξ0_PB(x,D0,Amat,Mmat)*ψ1_PB(x,D0,Amat,dvec,Mmat,z,v)-ξ1_PB(x,D0,Amat,dvec,Mmat,z,v))/ψ0_PB(x,D0,Amat) # scalar

# cosθ_y = d1(y)y + d2(y)y^2 + ...
d1_PB(x,D0,Amat,Mmat) = (dot(χ0_PB(x,D0,Amat),c1_PB(x,D0,Amat,Mmat))-ξ0_PB(x,D0,Amat,Mmat))/ψ0_PB(x,D0,Amat) # scalar
d2_PB(x,D0,Amat,dvec,Mmat,z,v) = (dot(χ0_PB(x,D0,Amat),c2_PB(x,D0,Amat,dvec,Mmat,z,v))+dot(χ1_PB(x,D0,Amat,dvec,Mmat,z,v),c1_PB(x,D0,Amat,Mmat))-ξ1_PB(x,D0,Amat,dvec,Mmat,z,v)-(dot(χ0_PB(x,D0,Amat),c1_PB(x,D0,Amat,Mmat))-ξ0_PB(x,D0,Amat,Mmat))*ψ1_PB(x,D0,Amat,dvec,Mmat,z,v)/ψ0_PB(x,D0,Amat))/ψ0_PB(x,D0,Amat) # scalar

# exp(-k|x-y|) = 1 + e1(y)y + e2(y)y^2 + ... = 1 -k ψ0(y) y + (k^2/2 ψ0^2(y)-k ψ1(y))y^2 + ...
# |x-y|exp(-k|x-y|) = f1(y)y + f2(y)y^2 + ... = ψ0(y)y + (ψ0(y)e1(y)+ψ1(y))y^2 + ... =
# = ψ0(y)y + (-k ψ0^2(y)+ψ1(y))y^2 + ... 

# exp(-k|x-y|)cos\theta_x cos\theta_y/|x-y| = g1(y)y + g2(y)y^2 + ...
g1_PB(x,D0,Amat,Mmat) = d1_PB(x,D0,Amat,Mmat)*b1_PB(x,D0,Amat,Mmat)/ψ0_PB(x,D0,Amat) # scalar
g2_PB(x,D0,Amat,dvec,Mmat,z,v,k) = (d1_PB(x,D0,Amat,Mmat)*b2_PB(x,D0,Amat,dvec,Mmat,z,v)+d2_PB(x,D0,Amat,dvec,Mmat,z,v)*b1_PB(x,D0,Amat,Mmat)-k*ψ0_PB(x,D0,Amat)*d1_PB(x,D0,Amat,Mmat)*b1_PB(x,D0,Amat,Mmat))/ψ0_PB(x,D0,Amat) - d1_PB(x,D0,Amat,Mmat)*b1_PB(x,D0,Amat,Mmat)*ψ1_PB(x,D0,Amat,dvec,Mmat,z,v)/(ψ0_PB(x,D0,Amat)^2) # scalar

# cos\theta_y(1-ep^-1 exp(-k|x-y|)(1+k|x-y|))/|x-y|^2 = i0/|y| + i1 + ...
i0_PB(x,D0,Amat,Mmat,ep) = (1-ep)*d1_PB(x,D0,Amat,Mmat)/ψ0_PB(x,D0,Amat)^2
i0_PB_a(x,D0,Amat,Mmat,ep) = i0_PB([cos(x);sin(x)],D0,Amat,Mmat,ep)
i1_PB(x,D0,Amat,dvec,Mmat,z,v,ep) = ((1-ep)*d2_PB(x,D0,Amat,dvec,Mmat,z,v))/ψ0_PB(x,D0,Amat)^2 -2(1-ep)d1_PB(x,D0,Amat,Mmat)*ψ1_PB(x,D0,Amat,dvec,Mmat,z,v)/ψ0_PB(x,D0,Amat)^3
i1_PB_a(x,D0,Amat,dvec,Mmat,z,v,ep) = i1_PB([cos(x);sin(x)],D0,Amat,dvec,Mmat,z,v,ep)

# cos\theta_x(1-ep^-1 exp(-k|x-y|)(1+k|x-y|))/|x-y|^2 = j0/|y| + j1 + ...
j0_PB(x,D0,Amat,Mmat,ep) = (1-ep)*b1_PB(x,D0,Amat,Mmat)/ψ0_PB(x,D0,Amat)^2
j0_PB_a(x,D0,Amat,Mmat,ep) = j0_PB([cos(x);sin(x)],D0,Amat,Mmat,ep)
j1_PB(x,D0,Amat,dvec,Mmat,z,v,ep) = ((1-ep)*b2_PB(x,D0,Amat,dvec,Mmat,z,v))/ψ0_PB(x,D0,Amat)^2 -2(1-ep)b1_PB(x,D0,Amat,Mmat)*ψ1_PB(x,D0,Amat,dvec,Mmat,z,v)/ψ0_PB(x,D0,Amat)^3
j1_PB_a(x,D0,Amat,dvec,Mmat,z,v,ep) = j1_PB([cos(x);sin(x)],D0,Amat,dvec,Mmat,z,v,ep)

# simplified function for i0 and j0: note that the (1-ep) is not included in this expression
ij0_PB_a(x,D0,Amat,Mmat) = 0.5*(dot(D0*Amat*[cos(x);sin(x)], Mmat*D0*Amat*[cos(x);sin(x)]))/(norm(D0*Amat*[cos(x);sin(x)])^3)

# (1-exp(-k|x-y|))/|x-y| = ell0 + ell1 |y| + ...
ℓ0_PB(k) = k
ℓ1_PB(x,D0,Amat,k) = -0.5*k^2*ψ0_PB(x,D0,Amat)#(1+k)*ψ1_PB(x,D0,Amat,dvec,Mmat,z,v)/ψ0_PB(x,D0,Amat)-0.5*k^2*ψ0_PB(x,D0,Amat)
ℓ1_PB_a(x,D0,Amat,k) = ℓ1_PB([cos(x);sin(x)],D0,Amat,k)

#########################
### rotations along x,y,z axes
Rx(a) = [1 0 0;0 cos(a) sin(a);0 -sin(a) cos(a)]
Ry(a) = [cos(a) 0 sin(a);0 1 0;-sin(a) 0 cos(a)]
Rz(a) = [cos(a) -sin(a) 0;sin(a) cos(a) 0;0 0 1]

#########################################
"The exponential kernel, supported in (-1,1), with one vanishing moment"
function Kinf(r::Real)
    a=7.513931532835806
    # println("Kinf r is=",r)
    if abs(r)<1
      return a*exp(2.0/(r^2-1.0))
    else
      return 0.0
    end
end
#########################################


#########################################
# functions for quadruple correction
function cpm_derivatives(w, tau1, tau2, nvec, Pgamma, Mmat; z::Real=0.02, h::Real=0.01, k1::Real=-1, k2::Real=-1)
    # Pmat = zeros(3,3,3);
    # for i=1:3
    #     for j=1:3
    #         Pmat[i,j,:] = Pgamma(w.+((i-2)*h)*tau1.+((j-2)*h)*tau2.+z*nvec)
    #     end
    # end
    Pmat = zeros(5,5,3);
    for i=1:5
        for j=1:5
            # println(norm(w .+ ((i-3)*h)*tau1 .+ ((j-3)*h)*tau2 .+ z*nvec - (w .+ ((i-3)*h)*tau1 .+ ((j-2)*h)*tau2 .+ z*nvec)))
            # println(h)
            # println(v," ",w .+ ((i-3)*h)*tau1 .+ ((j-3)*h)*tau2 .+ z*nvec)
            v = Pgamma(w .+ ((i-3)*h)*tau1 .+ ((j-3)*h)*tau2 .+ z*nvec)
            Pmat[i,j,:] = Mmat*(v-w)
            # Pmat[i,j,:] = Pgamma(v)
        end
    end

    Dx = [1/12;-2/3;0;2/3;-1/12]
    Dxx = [-1/12;4/3;-5/2;4/3;-1/12]
    # Xx = (Pmat[3,2,1]-Pmat[1,2,1])/(2*h)
    # Xy = (Pmat[2,3,1]-Pmat[2,1,1])/(2*h)
    # Yx = (Pmat[3,2,2]-Pmat[1,2,2])/(2*h)
    # Yy = (Pmat[2,3,2]-Pmat[2,1,2])/(2*h)
    Xx = dot(Pmat[:,3,1],Dx)/h
    Xy = dot(Pmat[3,:,1],Dx)/h
    Yx = dot(Pmat[:,3,2],Dx)/h
    Yy = dot(Pmat[3,:,2],Dx)/h

    # Xxx = (Pmat[3,2,1]-2*Pmat[2,2,1]+Pmat[1,2,1])/(h^2)
    # Xyy = (Pmat[2,3,1]-2*Pmat[2,2,1]+Pmat[2,1,1])/(h^2)
    # Yxx = (Pmat[3,2,2]-2*Pmat[2,2,2]+Pmat[1,2,2])/(h^2)
    # Yyy = (Pmat[2,3,2]-2*Pmat[2,2,2]+Pmat[2,1,2])/(h^2)
    Xxx = dot(Dxx,Pmat[:,3,1])/(h^2)
    Xyy = dot(Dxx,Pmat[3,:,1])/(h^2)
    Yxx = dot(Dxx,Pmat[:,3,2])/(h^2)
    Yyy = dot(Dxx,Pmat[3,:,2])/(h^2)


    # Xxy = (Pmat[3,3,1]-Pmat[1,3,1]-Pmat[3,1,1]+Pmat[1,1,1])/(4*h^2)
    # Yxy = (Pmat[3,3,2]-Pmat[1,3,2]-Pmat[3,1,2]+Pmat[1,1,2])/(4*h^2)
    Xxy = dot(Dx, [dot(Dx,Pmat[:,1,1]);dot(Dx,Pmat[:,2,1]);dot(Dx,Pmat[:,3,1]);dot(Dx,Pmat[:,4,1]);dot(Dx,Pmat[:,5,1])])/(h^2)
    Yxy = dot(Dx, [dot(Dx,Pmat[:,1,2]);dot(Dx,Pmat[:,2,2]);dot(Dx,Pmat[:,3,2]);dot(Dx,Pmat[:,4,2]);dot(Dx,Pmat[:,5,2])])/(h^2)
    Xyx = Xxy
    Yyx = Yxy

    Mat1 = [Xx^2 2*Xx*Yx Yx^2 0; Xx*Xy (Xx*Yy+Xy*Yx) Yx*Yy 0; Xx*Xy (Xx*Yy+Xy*Yx) Yx*Yy 0; Xy^2 2*Xy*Yy Yy^2 0; 0 Xx^2 2*Xx*Yx Yx^2; 0 Xx*Xy (Xx*Yy+Xy*Yx) Yx*Yy; 0 Xx*Xy (Xx*Yy+Xy*Yx) Yx*Yy; 0 Xy^2 2*Xy*Yy Yy^2]

    MMat1 = Mat1'*Mat1

    Mat2 = [Xx^2 2*Xx*Yx Yx^2 0; Xy^2 2*Xy*Yy Yy^2 0; 0 Xx^2 2*Xx*Yx Yx^2; 0 Xy^2 2*Xy*Yy Yy^2]
    Mat3 = [Xx^2 2*Xx*Yx Yx^2 0; Xx*Xy (Xx*Yy+Xy*Yx) Yx*Yy 0; 0 Xx*Xy (Xx*Yy+Xy*Yx) Yx*Yy; 0 Xy^2 2*Xy*Yy Yy^2]

    bvec = [(1-z*k1)*Xxx; (1-z*k1)*Xxy; (1-z*k1)*Xyx; (1-z*k1)*Xyy; (1-z*k2)*Yxx; (1-z*k2)*Yxy; (1-z*k2)*Yyx; (1-z*k2)*Yyy]./z

    bvec1 = Mat1'*bvec

    bvec2 = bvec[[1;4;5;8]]
    bvec3 = bvec[[1;2;7;8]]

    v1 = MMat1\bvec1
    v2 = Mat2\bvec2
    v3 = Mat3\bvec3
    # println("v1-v2=",abs.(v1-v2))
    # println("v1-v3=",abs.(v1-v3))
    # println("v2-v3=",abs.(v2-v3))
    # return v1, v2, v3
    return v2
end
function cpm_derivatives_grid_all(w, tau1, tau2, nvec, Pgamma, Mmat; z::Real=0.02, h::Real=0.01, k1::Real=-1, k2::Real=-1)
    # Pmat = zeros(3,3,3);
    # for i=1:3
    #     for j=1:3
    #         Pmat[i,j,:] = Pgamma(w.+((i-2)*h)*tau1.+((j-2)*h)*tau2.+z*nvec)
    #     end
    # end
    Pmat = zeros(5,5,5,3);
    for i=1:5
        for j=1:5
            for k=1:5
                v = Pgamma(w .+ [i-3;j-3;k-3]*h .+ z*nvec)
                Pmat[i,j,k,:] = Mmat*(v-w)
                # Pmat[i,j,k,:] = v
            end
        end
    end

    Dx = [1/12;-2/3;0;2/3;-1/12]
    Dxx = [-1/12;4/3;-5/2;4/3;-1/12]
    Xx = dot(Dx,Pmat[:,3,3,1])/h
    Xy = dot(Dx,Pmat[3,:,3,1])/h
    Xz = dot(Dx,Pmat[3,3,:,1])/h

    Yx = dot(Dx,Pmat[:,3,3,2])/h
    Yy = dot(Dx,Pmat[3,:,3,2])/h
    Yz = dot(Dx,Pmat[3,3,:,2])/h
    
    # Zx = dot(Pmat[:,3,3,3],Dx)/h
    # Zy = dot(Pmat[3,:,3,3],Dx)/h
    # Zz = dot(Pmat[3,3,:,3],Dx)/h
    
    Xxx = dot(Dxx,Pmat[:,3,3,1])/(h^2)
    Xyy = dot(Dxx,Pmat[3,:,3,1])/(h^2)
    Xzz = dot(Dxx,Pmat[3,3,:,1])/(h^2)
    
    Yxx = dot(Dxx,Pmat[:,3,3,2])/(h^2)
    Yyy = dot(Dxx,Pmat[3,:,3,2])/(h^2)
    Yzz = dot(Dxx,Pmat[3,3,:,2])/(h^2)
    
    # Zxx = dot(Dxx,Pmat[:,3,3,3])/(h^2)
    # Zyy = dot(Dxx,Pmat[3,:,3,3])/(h^2)
    # Zzz = dot(Dxx,Pmat[3,3,:,3])/(h^2)
    
    Xxy = dot(Dx, [dot(Dx,Pmat[:,1,3,1]);dot(Dx,Pmat[:,2,3,1]);dot(Dx,Pmat[:,3,3,1]);dot(Dx,Pmat[:,4,3,1]);dot(Dx,Pmat[:,5,3,1])])/(h^2)
    Yxy = dot(Dx, [dot(Dx,Pmat[:,1,3,2]);dot(Dx,Pmat[:,2,3,2]);dot(Dx,Pmat[:,3,3,2]);dot(Dx,Pmat[:,4,3,2]);dot(Dx,Pmat[:,5,3,2])])/(h^2)
    Xyx = dot(Dx, [dot(Dx,Pmat[1,:,3,1]);dot(Dx,Pmat[2,:,3,1]);dot(Dx,Pmat[3,:,3,1]);dot(Dx,Pmat[4,:,3,1]);dot(Dx,Pmat[5,:,3,1])])/(h^2)
    Yyx = dot(Dx, [dot(Dx,Pmat[1,:,3,2]);dot(Dx,Pmat[2,:,3,2]);dot(Dx,Pmat[3,:,3,2]);dot(Dx,Pmat[4,:,3,2]);dot(Dx,Pmat[5,:,3,2])])/(h^2)
    # Zxy = dot(Dx, [dot(Dx,Pmat[:,1,3,3]);dot(Dx,Pmat[:,2,3,3]);dot(Dx,Pmat[:,3,3,3]);dot(Dx,Pmat[:,4,3,3]);dot(Dx,Pmat[:,5,3,3])])/(h^2)
    
    Xxz = dot(Dx, [dot(Dx,Pmat[:,3,1,1]);dot(Dx,Pmat[:,3,2,1]);dot(Dx,Pmat[:,3,3,1]);dot(Dx,Pmat[:,3,4,1]);dot(Dx,Pmat[:,3,5,1])])/(h^2)
    Yxz = dot(Dx, [dot(Dx,Pmat[:,3,1,2]);dot(Dx,Pmat[:,3,2,2]);dot(Dx,Pmat[:,3,3,2]);dot(Dx,Pmat[:,3,4,2]);dot(Dx,Pmat[:,3,5,2])])/(h^2)
    Xzx = dot(Dx, [dot(Dx,Pmat[1,3,:,1]);dot(Dx,Pmat[2,3,:,1]);dot(Dx,Pmat[3,3,:,1]);dot(Dx,Pmat[4,3,:,1]);dot(Dx,Pmat[5,3,:,1])])/(h^2)
    Yzx = dot(Dx, [dot(Dx,Pmat[1,3,:,2]);dot(Dx,Pmat[2,3,:,2]);dot(Dx,Pmat[3,3,:,2]);dot(Dx,Pmat[4,3,:,2]);dot(Dx,Pmat[5,3,:,2])])/(h^2)
    # Zxz = dot(Dx, [dot(Dx,Pmat[:,3,1,3]);dot(Dx,Pmat[:,3,2,3]);dot(Dx,Pmat[:,3,3,3]);dot(Dx,Pmat[:,3,4,3]);dot(Dx,Pmat[:,3,5,3])])/(h^2)
    
    Xyz = dot(Dx, [dot(Dx,Pmat[3,:,1,1]);dot(Dx,Pmat[3,:,2,1]);dot(Dx,Pmat[3,:,3,1]);dot(Dx,Pmat[3,:,4,1]);dot(Dx,Pmat[3,:,5,1])])/(h^2)
    Yyz = dot(Dx, [dot(Dx,Pmat[3,:,1,2]);dot(Dx,Pmat[3,:,2,2]);dot(Dx,Pmat[3,:,3,2]);dot(Dx,Pmat[3,:,4,2]);dot(Dx,Pmat[3,:,5,2])])/(h^2)
    Xzy = dot(Dx, [dot(Dx,Pmat[3,1,:,1]);dot(Dx,Pmat[3,2,:,1]);dot(Dx,Pmat[3,3,:,1]);dot(Dx,Pmat[3,4,:,1]);dot(Dx,Pmat[3,5,:,1])])/(h^2)
    Yzy = dot(Dx, [dot(Dx,Pmat[3,1,:,2]);dot(Dx,Pmat[3,2,:,2]);dot(Dx,Pmat[3,3,:,2]);dot(Dx,Pmat[3,4,:,2]);dot(Dx,Pmat[3,5,:,2])])/(h^2)
    # Zyz = dot(Dx, [dot(Dx,Pmat[3,:,1,3]);dot(Dx,Pmat[3,:,2,3]);dot(Dx,Pmat[3,:,3,3]);dot(Dx,Pmat[3,:,4,3]);dot(Dx,Pmat[3,:,5,3])])/(h^2)
    
    # Xzy = Xyz; Yzy = Yyz
    # Xzx = Xxz; Yzx = Yxz
    # Xyx = Xxy; Yyx = Yxy

    DXvec = [Xx; Xy; Xz]; DYvec = [Yx; Yy; Yz];# DZvec = [Zx; Zy; Zz]
    DXmat = [Xxx Xyx Xzx;Xxy Xyy Xzy;Xxz Xyz Xzz]
    DYmat = [Yxx Yyx Yzx;Yxy Yyy Yzy;Yxz Yyz Yzz]
    
    # DXmat = [Xxx Xxy Xxz;Xyx Xyy Xyz;Xzx Xzy Xzz]
    # DYmat = [Yxx Yxy Yxz;Yyx Yyy Yyz;Yzx Yzy Yzz]

    Xtx = dot(DXvec, tau1);
    Xty = dot(DXvec, tau2);
    Ytx = dot(DYvec, tau1);
    Yty = dot(DYvec, tau2);
    
    Xtxx = dot(tau1,DXmat*tau1)
    Xtyy = dot(tau2,DXmat*tau2)
    Ytxx = dot(tau1,DYmat*tau1)
    Ytyy = dot(tau2,DYmat*tau2)

    Xtxy = dot(tau2,DXmat*tau1)
    Ytxy = dot(tau2,DYmat*tau1)
    Xtyx = dot(tau1,DXmat*tau2)
    Ytyx = dot(tau1,DYmat*tau2)
    # println("X mixed second derivatives:")
    # println(abs(Xtxy-Xtyx))
    # println(abs(Ytxy-Ytyx),"\n")
    # Xyx = Xxy
    # Yyx = Yxy

    MmatX = [Xtx^2 Xtx*Ytx Xtx*Ytx Ytx^2;Xtx*Xty Xtx*Yty Xty*Ytx Ytx*Yty; Xtx*Xty Xty*Ytx Xtx*Yty Ytx*Yty; Xty^2 Xty*Yty Xty*Yty Yty^2]
    # println(cond(MmatX))

    bvec1 = [(1-z*k1)*Xtxx; (1-z*k1)*Xtxy; (1-z*k1)*Xtyx; (1-z*k1)*Xtyy]./z
    bvec2 = [(1-z*k2)*Ytxx; (1-z*k2)*Ytxy; (1-z*k2)*Ytyx; (1-z*k2)*Ytyy]./z

    v1 = MmatX\bvec1
    v2 = MmatX\bvec2
    # v3 = Mat3\bvec3
    # println("v1-v2=",abs.(v1-v2))
    # println("v1-v3=",abs.(v1-v3))
    # println("v2-v3=",abs.(v2-v3))
    # return v1, v2, v3
    return [v1;v2]
end
#########################################

#########################################
#########################
# projection mappings needed

function Pgamma_torus(z,Rv2,Rinv,R1,R2,xshift)
    zn = Rinv*(z.-xshift);
    phi = atan(zn[2],zn[1]);
    ctmp = R1*[cos(phi);sin(phi);0];
    return Rv2*((zn.-ctmp)/norm(zn.-ctmp)*R2.+ctmp).+xshift
end
function insidepoint_torus(z,Rinv,R1,R2,xshift)
    zn = Rinv*(z.-xshift);
    phi = atan(zn[2],zn[1]);
    ctmp = R1*[cos(phi);sin(phi);0];
    return norm(ctmp-zn)<R2
end

function Pgamma_fgen_4(X2::Array{T,1} where T<:Real;
    Rv2::Array{T,2} where T<:Real=[1 0 0;0 1 0;0 0 1],
    Rinv::Array{T,2} where T<:Real=[1 0 0;0 1 0;0 0 1],
    xshift::Array{T,1} where T<:Real=[0;0;0],
    k1::T where T<:Real=-1,
    k2::T where T<:Real=-1,
    fxxx::T where T<:Real=0,
    fyyy::T where T<:Real=0,
    fxxy::T where T<:Real=0,
    fxyy::T where T<:Real=0,
    fx4::T where T<:Real=0,
    fx3y::T where T<:Real=0,
    fx2y2::T where T<:Real=0,
    fxy3::T where T<:Real=0,
    fy4::T where T<:Real=0, TOL::Real=1e-9
    )

    local f(x,y) = 0.5*(k1*x^2+k2*y^2+fxxx*x^3/3+fyyy*y^3/3+fxxy*x^2*y+fxyy*x*y^2)+(fx4*x^4+fx3y*x^3*y+fx2y2*x^2*y^2+fxy3*x*y^3+fy4*y^4);
    local dfdx(x,y) = k1*x+0.5*(fxxx*x^2+2*fxxy*x*y+fxyy*y^2)+(4fx4*x^3+3fx3y*x^2*y+2fx2y2*x*y^2+fxy3*y^3)
    local dfdy(x,y) = k2*y+0.5*(fxxy*x^2+2*fxyy*x*y+fyyy*y^2)+(fx3y*x^3+2fx2y2*x^2*y+3fxy3*x*y^2+4fy4*y^3)
    local d2fdx2(x,y) = k1+fxxx*x+fxxy*y+(12fx4*x^2+6fx3y*x*y+2fx2y2*y^2)
    local d2fdy2(x,y) = k2+fyyy*y+fxyy*x+(2fx2y2*x^2+6fxy3*x*y+12fy4*y^2)
    local d2fdxy(x,y) = fxxy*x+fxyy*y+(3fx3y*x^2+4fx2y2*x*y+3fxy3*y^2)

    X = Rinv*(X2.-xshift);

    local F(Xb) = [Xb[1]-(X[3]-f(Xb[1],Xb[2]))*dfdx(Xb[1],Xb[2])-X[1]; Xb[2]-(X[3]-f(Xb[1],Xb[2]))*dfdy(Xb[1],Xb[2])-X[2]];
    local DF(Xb) = [1+(dfdx(Xb[1],Xb[2]))^2-(X[3]-f(Xb[1],Xb[2]))*d2fdx2(Xb[1],Xb[2]) dfdy(Xb[1],Xb[2])*dfdx(Xb[1],Xb[2])-(X[3]-f(Xb[1],Xb[2]))*d2fdxy(Xb[1],Xb[2]);dfdy(Xb[1],Xb[2])*dfdx(Xb[1],Xb[2])-(X[3]-f(Xb[1],Xb[2]))*d2fdxy(Xb[1],Xb[2]) 1+(dfdy(Xb[1],Xb[2]))^2-(X[3]-f(Xb[1],Xb[2]))*d2fdy2(Xb[1],Xb[2])]
    # local F(X) = [(X[1]-Z[1])-k1*X[1]*(-k1*0.5*X[1]^2-k2*0.5*X[2]^2+η-a*X[1]/c-b*X[2]/c);(X[2]-Z[2])-k2*X[2]*(-k1*0.5*X[1]^2-k2*0.5*X[2]^2+η-a*X[1]/c-b*X[2]/c)];
    # DF(X) = [1+k1*(k1*0.5*X[1]^2+k2*0.5*X[2]^2-Z[3])+k1^2*X[1]^2 k1*k2*X[1]*X[2];k1*k2*X[1]*X[2] 1+k2*(k1*0.5*X[1]^2+k2*0.5*X[2]^2-Z[3])+k2^2*X[2]^2];
    local d = 1;
    # TOL = 1e-11;
    # local X = Z[1:2];
    # local X = [Z[1];2*sqrt(Z[3]-k1*0.5*Z[1]^2)/k2];
    local Xb = [X[1]/(1-X[3]*k1);X[2]/(1-X[3]*k2)];
    # println("X=$X, res=$(norm(F(X),Inf))")
    # X = [0.;0.]
    local niter=0;
    while (norm(d,Inf)>TOL && niter<25)
        d = lu(DF(Xb))\F(Xb);
        Xb -= d;
        niter += 1;
        # println("i=$niter, X=$X, d=$(norm(d,Inf)), res=$(norm(F(X),Inf))")
    end
    if niter==25
        @warn println("No convergence for Newton's method, residual=$d")
    end
    # println("x,y=$(Xb[1]), $(Xb[2])")
    # return [X[1];X[2];0.5*(k1*X[1]^2+k2*X[2]^2)]
    return xshift.+Rv2*[Xb;f(Xb[1],Xb[2])]
end

function insidepoint_fgen_4(X2::Array{T,1} where T<:Real;
    # Rv2::Array{T,2} where T<:Real=[1 0 0;0 1 0;0 0 1],
    Rinv::Array{T,2} where T<:Real=[1 0 0;0 1 0;0 0 1],
    xshift::Array{T,1} where T<:Real=[0;0;0],
    k1::T where T<:Real=-1,
    k2::T where T<:Real=-1,
    fxxx::T where T<:Real=0,
    fyyy::T where T<:Real=0,
    fxxy::T where T<:Real=0,
    fxyy::T where T<:Real=0,
    fx4::T where T<:Real=0,
    fx3y::T where T<:Real=0,
    fx2y2::T where T<:Real=0,
    fxy3::T where T<:Real=0,
    fy4::T where T<:Real=0
    )
    local f(x,y) = 0.5*(k1*x^2+k2*y^2+fxxx*x^3/3+fyyy*y^3/3+fxxy*x^2*y+fxyy*x*y^2)+(fx4*x^4+fx3y*x^3*y+fx2y2*x^2*y^2+fxy3*x*y^3+fy4*y^4);
    
    X = Rinv*(X2.-xshift);
    return (X[3]<f(X[1],X[2]))
end

function Pgamma_circle(x,R,x0)
  return (norm(x.-x0)>1e-3)*((x.-x0)*R/norm(x.-x0).+x0) + (norm(x.-x0)<=1e-3)*(x0 .+ R*x0/norm(x0))
end

# functions for CPM onto two touching spheres

function isInTriangle(p1,p2,p3,x)
  v = [p1[1]-p3[1] p2[1]-p3[1];p1[2]-p3[2] p2[2]-p3[2]]\[x[1]-p3[1];x[2]-p3[2]];
  return ((v[1]>=0) && (v[2]>=0) && ((v[1]+v[2])<=1))
end

function Pgammafun_spheres(z)
  # println(z)
  Rvec = [1.0;0.7;0.3];
  v= [-Rvec[1];0;0]; 
  w = [Rvec[2];0;0];

  if z[1]>=-Rvec[1] && z[1]<=Rvec[2]
      L1 = Rvec[3]+Rvec[2];
      L2 = Rvec[2]+Rvec[1];
      L3 = Rvec[1]+Rvec[3];
      angles = [ acos((L2^2+L3^2-L1^2)/(2*L2*L3)); acos((L1^2+L3^2-L2^2)/(2*L1*L3)); acos((L1^2+L2^2-L3^2)/(2*L1*L2))]
  
      # finding plane passing through v,w and z
      Amat = [w[1]-z[1] w[3]-z[3];-z[1]+v[1] v[3]-z[3]]
      avec = Amat\[-w[2]+z[2];-v[2]+z[2]]
  
      # finding direction normal to w-v, normalized
      s3 = z[3];
      Bmat = [avec[1] 1;w[1]-v[1] w[2]-v[2]]
      bvec = [z[1]*avec[1]+avec[2]*(z[3]-s3)+z[2];(w[1]-v[1])*v[1]+(w[2]-v[2])*v[2]+(v[3]-s3)*(w[3]-v[3])]# +cos(angles[1])
      svec = Bmat\bvec; s = [svec[1];svec[2];s3].-v; s /= norm(s)
  
      # finding point where triangle height touches the base vw
      par = (L3^2-L1^2)/(2*L2^2)+0.5
      vwx = v + par*(w-v)
      
      heig = sqrt(abs(L1^2-L2^2*(1-par)^2)) # 2eight
      u = vwx .+ (s)*heig # center of sphere rolling above the generated surface
      
      # 2D basis of the plane, origin in v.
      e1 = (w.-v)/norm(w.-v)
      e2 = s
      cvec = [e1[1] e2[1]; e1[2] e2[2]]\[u[1]-v[1];u[2]-v[2]]
      z2D = [e1[1] e2[1]; e1[2] e2[2]]\[z[1]-v[1];z[2]-v[2]]
      p1 = [0;0]
      p2 = [norm(w-v);0]
      p3 = cvec
  
      if isInTriangle(p1,p2,p3,z2D)
          return Pgamma_circle(z,Rvec[3],u)
      elseif isInTriangle(p1,p2,-p3,z2D)
          return Pgamma_circle(z,Rvec[3],vwx .- (s)*heig)
      else
          pz1 = Pgamma_circle(z,Rvec[1],v)
          pz2 = Pgamma_circle(z,Rvec[2],w)
          return (norm(z-pz1)<norm(z-pz2))*pz1 .+ (norm(z-pz1)>=norm(z-pz2))*pz2
      end
  else
      pz1 = Pgamma_circle(z,Rvec[1],v)
      pz2 = Pgamma_circle(z,Rvec[2],w)
      return (norm(z-pz1)<norm(z-pz2))*pz1 .+ (norm(z-pz1)>=norm(z-pz2))*pz2
  end
end

function insidepoint_spheres(z)
  Rvec = [1.0;0.7;0.3];
  v= [-Rvec[1];0;0]; 
  w = [Rvec[2];0;0];

  if z[1]>=-Rvec[1] && z[1]<=Rvec[2]
      L1 = Rvec[3]+Rvec[2];
      L2 = Rvec[2]+Rvec[1];
      L3 = Rvec[1]+Rvec[3];
      angles = [ acos((L2^2+L3^2-L1^2)/(2*L2*L3)); acos((L1^2+L3^2-L2^2)/(2*L1*L3)); acos((L1^2+L2^2-L3^2)/(2*L1*L2))]
  
      # finding plane passing through v,w and z
      Amat = [w[1]-z[1] w[3]-z[3];-z[1]+v[1] v[3]-z[3]]
      avec = Amat\[-w[2]+z[2];-v[2]+z[2]]
  
      # finding direction normal to w-v, normalized
      s3 = z[3];
      Bmat = [avec[1] 1;w[1]-v[1] w[2]-v[2]]
      bvec = [z[1]*avec[1]+avec[2]*(z[3]-s3)+z[2];(w[1]-v[1])*v[1]+(w[2]-v[2])*v[2]+(v[3]-s3)*(w[3]-v[3])]# +cos(angles[1])
      svec = Bmat\bvec; s = [svec[1];svec[2];s3].-v; s /= norm(s)
  
      # finding point where triangle height touches the base vw
      par = (L3^2-L1^2)/(2*L2^2)+0.5
      vwx = v + par*(w-v)
      
      heig = sqrt(abs(L1^2-L2^2*(1-par)^2)) # height
      u = vwx .+ (s)*heig # center of sphere rolling above the generated surface
      
      # 2D basis of the plane, origin in v.
      e1 = (w.-v)/norm(w.-v)
      e2 = s
      cvec = [e1[1] e2[1]; e1[2] e2[2]]\[u[1]-v[1];u[2]-v[2]]
      z2D = [e1[1] e2[1]; e1[2] e2[2]]\[z[1]-v[1];z[2]-v[2]]
      p1 = [0;0]
      p2 = [norm(w-v);0]
      p3 = cvec
      
      if isInTriangle(p1,p2,p3,z2D)
          return norm(z-u)>Rvec[3] #Pgamma_circle(z,Rvec[3],u)
      elseif isInTriangle(p1,p2,-p3,z2D)
          return norm(z .- vwx .+ (s)*heig)>Rvec[3] #Pgamma_circle(z,Rvec[3],vwx .- (s)*h)
      else
          pz1 = Pgamma_circle(z,Rvec[1],v)
          pz2 = Pgamma_circle(z,Rvec[2],w)
          return (norm(z-pz1)<norm(z-pz2))*(norm(z-v)<Rvec[1]) .+ (norm(z-pz1)>=norm(z-pz2))*(norm(z-w)<Rvec[2])
      end
  else
      pz1 = Pgamma_circle(z,Rvec[1],v)
      pz2 = Pgamma_circle(z,Rvec[2],w)
      return (norm(z-pz1)<norm(z-pz2))*(norm(z-v)<Rvec[1]) .+ (norm(z-pz1)>=norm(z-pz2))*(norm(z-w)<Rvec[2])
  end
end

function get_jac_spheres(z,η)
  Rvec = [1.0;0.7;0.3];
  v= [-Rvec[1];0;0]; 
  w = [Rvec[2];0;0];

  if z[1]>=-Rvec[1] && z[1]<=Rvec[2]
      L1 = Rvec[3]+Rvec[2];
      L2 = Rvec[2]+Rvec[1];
      L3 = Rvec[1]+Rvec[3];
      angles = [ acos((L2^2+L3^2-L1^2)/(2*L2*L3)); acos((L1^2+L3^2-L2^2)/(2*L1*L3)); acos((L1^2+L2^2-L3^2)/(2*L1*L2))]
  
      # finding plane passing through v,w and z
      Amat = [w[1]-z[1] w[3]-z[3];-z[1]+v[1] v[3]-z[3]]
      avec = Amat\[-w[2]+z[2];-v[2]+z[2]]
  
      # finding direction normal to w-v, normalized
      s3 = z[3];
      Bmat = [avec[1] 1;w[1]-v[1] w[2]-v[2]]
      bvec = [z[1]*avec[1]+avec[2]*(z[3]-s3)+z[2];(w[1]-v[1])*v[1]+(w[2]-v[2])*v[2]+(v[3]-s3)*(w[3]-v[3])]# +cos(angles[1])
      svec = Bmat\bvec; s = [svec[1];svec[2];s3].-v; s /= norm(s)
  
      # finding point where triangle height touches the base vw
      par = (L3^2-L1^2)/(2*L2^2)+0.5
      vwx = v + par*(w-v)
      
      heig = sqrt(abs(L1^2-L2^2*(1-par)^2)) # height
      u = vwx .+ (s)*heig # center of sphere rolling above the generated surface
      
      # 2D basis of the plane, origin in v.
      e1 = (w.-v)/norm(w.-v)
      e2 = s
      cvec = [e1[1] e2[1]; e1[2] e2[2]]\[u[1]-v[1];u[2]-v[2]]
      z2D = [e1[1] e2[1]; e1[2] e2[2]]\[z[1]-v[1];z[2]-v[2]]
      p1 = [0;0]
      p2 = [norm(w-v);0]
      p3 = cvec
      
      if isInTriangle(p1,p2,p3,z2D)
          η = (1-2*(norm(z-u)>Rvec[3]))*(norm(z-Pgamma_circle(z,Rvec[3],u))) 
          return 1+2η*(1/(2Rvec[3]))+(η^2)*(1/Rvec[3])^2
          #Pgamma_circle(z,Rvec[3],u)
      elseif isInTriangle(p1,p2,-p3,z2D)
          return (1-2(norm(z-(vwx .- (s)*h)))>Rvec[3])*(norm(z-Pgamma_circle(z,Rvec[3],(vwx .- (s)*h)))) #Pgamma_circle(z,Rvec[3],vwx .- (s)*h)
          return 1+2η*(1/(2Rvec[3]))+η^2(1/Rvec[3]^2)
      else
          pz1 = Pgamma_circle(z,Rvec[1],v)
          pz2 = Pgamma_circle(z,Rvec[2],w)
          if pz1<pz2
              η = (1-2(norm(z-v)<Rvec[1]))*(norm(z-Pgamma_circle(z,Rvec[1],v)))
              return 1+2η*(-1/(2Rvec[1]))+η^2(1/Rvec[1]^2)
          else
              η = (1-2(norm(z-w)<Rvec[2]))*(norm(z-Pgamma_circle(z,Rvec[2],w)))
              return 1+2η*(-1/(2Rvec[2]))+η^2(1/Rvec[2]^2)
          end
          # return (norm(z-pz1)<norm(z-pz2))*(norm(z-v)<Rvec[1]) .+ (norm(z-pz1)>=norm(z-pz2))*(norm(z-w)<Rvec[2])
      end
  else
      pz1 = Pgamma_circle(z,Rvec[1],v)
      pz2 = Pgamma_circle(z,Rvec[2],w)
      if pz1<pz2
          η = (1-2(norm(z-v)<Rvec[1]))*(norm(z-Pgamma_circle(z,Rvec[1],v)))
          return 1+2η*(-1/(2Rvec[1]))+η^2(1/Rvec[1]^2)
      else
          η = (1-2(norm(z-w)<Rvec[2]))*(norm(z-Pgamma_circle(z,Rvec[2],w)))
          return 1+2η*(-1/(2Rvec[2]))+η^2(1/Rvec[2]^2)
      end
  end
end

function far_insidepoint_spheres(z)
  Rvec = [1.0;0.7;0.3];
  v= [-Rvec[1];0;0]; 
  w = [Rvec[2];0;0];

  if z[1]>=-Rvec[1] && z[1]<=Rvec[2]
      L1 = Rvec[3]+Rvec[2];
      L2 = Rvec[2]+Rvec[1];
      L3 = Rvec[1]+Rvec[3];
      angles = [ acos((L2^2+L3^2-L1^2)/(2*L2*L3)); acos((L1^2+L3^2-L2^2)/(2*L1*L3)); acos((L1^2+L2^2-L3^2)/(2*L1*L2))]
  
      # finding plane passing through v,w and z
      Amat = [w[1]-z[1] w[3]-z[3];-z[1]+v[1] v[3]-z[3]]
      avec = Amat\[-w[2]+z[2];-v[2]+z[2]]
  
      # finding direction normal to w-v, normalized
      s3 = z[3];
      Bmat = [avec[1] 1;w[1]-v[1] w[2]-v[2]]
      bvec = [z[1]*avec[1]+avec[2]*(z[3]-s3)+z[2];(w[1]-v[1])*v[1]+(w[2]-v[2])*v[2]+(v[3]-s3)*(w[3]-v[3])]# +cos(angles[1])
      svec = Bmat\bvec; s = [svec[1];svec[2];s3].-v; s /= norm(s)
  
      # finding point where triangle height touches the base vw
      par = (L3^2-L1^2)/(2*L2^2)+0.5
      vwx = v + par*(w-v)
      
      heig = sqrt(abs(L1^2-L2^2*(1-par)^2)) # height
      u = vwx .+ (s)*heig # center of sphere rolling above the generated surface
      
      # 2D basis of the plane, origin in v.
      e1 = (w.-v)/norm(w.-v)
      e2 = s
      cvec = [e1[1] e2[1]; e1[2] e2[2]]\[u[1]-v[1];u[2]-v[2]]
      z2D = [e1[1] e2[1]; e1[2] e2[2]]\[z[1]-v[1];z[2]-v[2]]
      p1 = [0;0]
      p2 = [norm(w-v);0]
      p3 = cvec
      
      if isInTriangle(p1,p2,p3,z2D)
          return norm(z-u)>Rvec[3]+ε #Pgamma_circle(z,Rvec[3],u)
      elseif isInTriangle(p1,p2,-p3,z2D)
          return norm(z .- vwx .+ (s)*heig)>Rvec[3]+ε #Pgamma_circle(z,Rvec[3],vwx .- (s)*h)
      else
          pz1 = Pgamma_circle(z,Rvec[1],v)
          pz2 = Pgamma_circle(z,Rvec[2],w)
          return (norm(z-pz1)<norm(z-pz2))*(norm(z-v)<Rvec[1]-ε) .+ (norm(z-pz1)>=norm(z-pz2))*(norm(z-w)<Rvec[2]-ε)
      end
  else
      pz1 = Pgamma_circle(z,Rvec[1],v)
      pz2 = Pgamma_circle(z,Rvec[2],w)
      return (norm(z-pz1)<norm(z-pz2))*(norm(z-v)<Rvec[1]-ε) .+ (norm(z-pz1)>=norm(z-pz2))*(norm(z-w)<Rvec[2]-ε)
  end
end

function far_outsidepoint_spheres(z)
  Rvec = [1.0;0.7;0.3];
  v= [-Rvec[1];0;0]; 
  w = [Rvec[2];0;0];

  if z[1]>=-Rvec[1] && z[1]<=Rvec[2]
      L1 = Rvec[3]+Rvec[2];
      L2 = Rvec[2]+Rvec[1];
      L3 = Rvec[1]+Rvec[3];
      angles = [ acos((L2^2+L3^2-L1^2)/(2*L2*L3)); acos((L1^2+L3^2-L2^2)/(2*L1*L3)); acos((L1^2+L2^2-L3^2)/(2*L1*L2))]
  
      # finding plane passing through v,w and z
      Amat = [w[1]-z[1] w[3]-z[3];-z[1]+v[1] v[3]-z[3]]
      avec = Amat\[-w[2]+z[2];-v[2]+z[2]]
  
      # finding direction normal to w-v, normalized
      s3 = z[3];
      Bmat = [avec[1] 1;w[1]-v[1] w[2]-v[2]]
      bvec = [z[1]*avec[1]+avec[2]*(z[3]-s3)+z[2];(w[1]-v[1])*v[1]+(w[2]-v[2])*v[2]+(v[3]-s3)*(w[3]-v[3])]# +cos(angles[1])
      svec = Bmat\bvec; s = [svec[1];svec[2];s3].-v; s /= norm(s)
  
      # finding point where triangle height touches the base vw
      par = (L3^2-L1^2)/(2*L2^2)+0.5
      vwx = v + par*(w-v)
      
      heig = sqrt(abs(L1^2-L2^2*(1-par)^2)) # height
      u = vwx .+ (s)*heig # center of sphere rolling above the generated surface
      
      # 2D basis of the plane, origin in v.
      e1 = (w.-v)/norm(w.-v)
      e2 = s
      cvec = [e1[1] e2[1]; e1[2] e2[2]]\[u[1]-v[1];u[2]-v[2]]
      z2D = [e1[1] e2[1]; e1[2] e2[2]]\[z[1]-v[1];z[2]-v[2]]
      p1 = [0;0]
      p2 = [norm(w-v);0]
      p3 = cvec
      
      if isInTriangle(p1,p2,p3,z2D)
          return norm(z-u)<Rvec[3]-ε #Pgamma_circle(z,Rvec[3],u)
      elseif isInTriangle(p1,p2,-p3,z2D)
          return norm(z .- vwx .+ (s)*heig)<Rvec[3]-ε #Pgamma_circle(z,Rvec[3],vwx .- (s)*h)
      else
          pz1 = Pgamma_circle(z,Rvec[1],v)
          pz2 = Pgamma_circle(z,Rvec[2],w)
          return (norm(z-pz1)<norm(z-pz2))*(norm(z-v)>Rvec[1]+ε) .+ (norm(z-pz1)>=norm(z-pz2))*(norm(z-w)>Rvec[2]+ε)
      end
  else
      pz1 = Pgamma_circle(z,Rvec[1],v)
      pz2 = Pgamma_circle(z,Rvec[2],w)
      return (norm(z-pz1)<norm(z-pz2))*(norm(z-v)>Rvec[1]+ε) .+ (norm(z-pz1)>=norm(z-pz2))*(norm(z-w)>Rvec[2]+ε)
  end
end

function Pgammafun_2spheres(z, x1a, x2a, Rot3, Rottot, sintC1, costC1, sintC2, costC2, x3a, cone_compA, cone_compB, a2, b2, RadiiVec, x0Vec)
  p2a = Rx(-pi*0.5)*Rot3*(z.-x1a)
  p3a = Rx(pi*0.5)*Rot3*(z.-x2a)
  p2va = isInCone(sintC1, costC1, norm(x1a-x3a),p2a, cone_compA)
  p3va = isInCone(sintC2, costC2, norm(x2a-x3a),p3a, cone_compB)
  if ( (cone_compA*cone_compB>0) && (p2va || p3va) ) || ( (cone_compA<0) && (p3va && ~p2va) ) || ( (cone_compB<0) && (p2va && ~p3va) )
      RotX = Rottot*Rx(-pi*0.5)
      RotXinv = transpose(RotX)
      return Pgamma_torus(z,RotX,RotXinv,a2,b2,x3a)
  else
      Pz1 = Pgamma_circle(z, RadiiVec[1], x0Vec[1])
      Pz2 = Pgamma_circle(z, RadiiVec[2], x0Vec[2])
      if norm(z.-x1a)<RadiiVec[1]
          return Pz1
      elseif norm(z.-x2a)<RadiiVec[2]
          return Pz2
      else
          return (norm(z-Pz1)<norm(z-Pz2))*Pz1 .+ (norm(z-Pz1)>=norm(z-Pz2))*Pz2
      end
  end
end
function insidepoint_2spheres(z, x1a, x2a, Rot3, Rottot, sintC1, costC1, sintC2, costC2, x3a, cone_compA, cone_compB, a2, b2, RadiiVec)
  p2a = Rx(-pi*0.5)*Rot3*(z.-x1a)
  p3a = Rx(pi*0.5)*Rot3*(z.-x2a)
  p2va = isInCone(sintC1, costC1, norm(x1a-x3a),p2a, cone_compA)
  p3va = isInCone(sintC2, costC2, norm(x2a-x3a),p3a, cone_compB)
  if ( (cone_compA*cone_compB>0) && (p2va || p3va) ) || ( (cone_compA<0) && (p3va && ~p2va) ) || ( (cone_compB<0) && (p2va && ~p3va) )
      RotX = Rottot*Rx(-pi*0.5)
      RotXinv = transpose(RotX)
      return !insidepoint_torus(z,RotXinv,a2,b2,x3a)
      # return Pgamma_torus(z,RotX,RotXinv,a2,b2,x3a)
  elseif norm(z.-x1a)<RadiiVec[1]
      return true
  elseif norm(z.-x2a)<RadiiVec[2]
      return true
  else
      return false
  end
end

# functions for CPM onto single sphere

function Pgammafun_sphere(z,x0,R)
  return (norm(z-x0)>1e-8)*((z-x0)/norm(z-x0)*R+x0) + (norm(z-x0)<=1e-8)*(x0+R*[0.;0;1.]);
end
function normalz_sphere(z,x0,R)
  z = Pgammafun_sphere(z,x0,R);
  return (norm(z-x0)>1e-8)*((z-x0)/norm(z-x0)) + (norm(z-x0)<=1e-8)*([0;0;1.]);
end
function insidepoint_sphere(z,x0,R)
  # z = Pgammafun_sphere(z,x0,R)
  return norm(z-x0)<R
end
function far_insidepoint_sphere(z,x0,R,ε)
  # z = Pgammafun_sphere(z,x0,R)
  return norm(z-x0)<R-ε
end
function far_outsidepoint_sphere(z,x0,R,ε)
  # z = Pgammafun_sphere(z,x0,R)
  return norm(z-x0)>R+ε
end


#########################
# weight computation functions

function biquintic(X,Y,wXY,xy)
    c = [ones(1,6); Y[:]'; Y[:]'.^2; Y[:]'.^3; Y[:]'.^4; Y[:]'.^5] \ [1; xy[2]; xy[2]^2; xy[2]^3; xy[2]^4; xy[2]^5]
    b = [ones(1,6); X[:]'; X[:]'.^2; X[:]'.^3; X[:]'.^4; X[:]'.^5] \ [1; xy[1]; xy[1]^2; xy[1]^3; xy[1]^4; xy[1]^5]
    return b'*wXY*c
end

###################################
# Fourier coefficients for ω weights w1
###################################

# path = path#*"2020/Weights/S2"

allAlpha_s2 = Array{Float64,1}(undef,101)
filepath = path*"/allalpha_s2.dat"
allAlpha_s2 .= (readdlm(filepath,Float64))[:,1]
Na_s2 = length(allAlpha_s2)
h_Alpha_s2 = allAlpha_s2[2]-allAlpha_s2[1]

allBeta_s2 = Array{Float64,1}(undef,101)
filepath = path*"/allbeta_s2.dat"
allBeta_s2 .= (readdlm(filepath,Float64))[:,1]
Nb_s2 = length(allBeta_s2)
h_Beta_s2 = allBeta_s2[2]-allBeta_s2[1]

nf_s2_plus = 22;
Nf_s2_plus = 2*nf_s2_plus+1;

filepath = path*"/allweights_s2_ab_plus.dat"
w0_all_s2_ab_plus = readdlm(filepath)
w0_all_s2_ab_plus = reshape(w0_all_s2_ab_plus, Nf_s2_plus, Na_s2, Nb_s2)

# allAlpha_s2_w4 = Array{Float64,1}(undef,101)
# filepath = path*"/allalpha_s2_w4.dat"
# allAlpha_s2_w4 .= (readdlm(filepath,Float64))[:,1]
# Na_s2_w4 = length(allAlpha_s2_w4)
# h_Alpha_s2_w4 = allAlpha_s2_w4[2]-allAlpha_s2_w4[1]

# allBeta_s2_w4 = Array{Float64,1}(undef,101)
# filepath = path*"/allbeta_s2_w4.dat"
# allBeta_s2_w4 .= (readdlm(filepath,Float64))[:,1]
# Nb_s2_w4 = length(allBeta_s2_w4)
# h_Beta_s2_w4 = allBeta_s2_w4[2]-allBeta_s2_w4[1]

# filepath = path*"/allweights_s2_ab_w4_plus.dat"
# w0_all_s2_ab_w4_plus = readdlm(filepath)
# w0_all_s2_ab_w4_plus = reshape(w0_all_s2_ab_w4_plus, Nf_s2_plus, Na_s2, Nb_s2, 4)

# nf_s2_w4_plus = 22;
# Nf_s2_w4_plus = 2*nf_s2_w4_plus+1;

global AmatFI2 = ones(45,45)
global tmatFI2 = (LinRange(0,π,46))[1:45]
global hmatFI2 = tmatFI2[2]-tmatFI2[1];
tvec2 = (LinRange(0,pi,46))[1:45];
for i=1:22
    AmatFI2[:,2*i:2*i+1] = [cos.(2*(i)*tvec2) sin.(2*(i)*tvec2)]
end
global FmatFI2 = lu(AmatFI2)


function weightsinterp_s2biquintic_ab_plus(abArray::Array{Array{Float64,1},1})
  interpfun = biquintic
  approxweights = Array{Float64,2}(undef,Nf_s2_plus,length(abArray))
  for (iab,ab) in enumerate(abArray)
    α = ab[1]
    β = ab[2]

    i = Int(ceil((α-allAlpha_s2[1])/h_Alpha_s2))+1
    j = Int(ceil((β-allBeta_s2[1])/h_Beta_s2))+1
    i = min( Na_s2-3 , max( i , 3 ) )
    j = min( Nb_s2-3 , max( j , 3 ) )
    inti = i-2:i+3
    intj = j-2:j+3
    w0 = w0_all_s2_ab_plus[:,inti,intj]

    for k=1:Nf_s2_plus
      approxweights[k,iab] = interpfun(allAlpha_s2[inti],allBeta_s2[intj],w0[k,:,:],[α;β]);
    end
  end

  return approxweights
end

function weightsinterp_s2biquintic_ab_w4_plus(abArray::Array{Array{Float64,1},1})
    interpfun = biquintic

    approxweights = Array{Float64,3}(undef,Nf_s2_w4_plus,length(abArray), 4)
    for (iab,ab) in enumerate(abArray)
        α = ab[1]
        β = ab[2]

        i = Int(ceil((α-allAlpha_s2_w4[1])/h_Alpha_s2_w4))+1
        j = Int(ceil((β-allBeta_s2_w4[1])/h_Beta_s2_w4))+1
        i = min( Na_s2_w4-3 , max( i , 3 ) )
        j = min( Nb_s2_w4-3 , max( j , 3 ) )
        inti = i-2:i+3
        intj = j-2:j+3
        w0 = w0_all_s2_ab_w4_plus[:,inti,intj,:]

        for i4=1:4
        for k=1:Nf_s2_w4_plus
            approxweights[k,iab,i4] = interpfun(allAlpha_s2_w4[inti],allBeta_s2_w4[intj],w0[k,:,:,i4],[α;β]);
        end
        end
    end

    return approxweights
end

####################
#### interpolation, k=0, single and quadruple
####################

#### k=0, single correction, interpolation
function w_k0_ptilde1(ab; lfun::Function=(x->ones(size(x))) )
    omall_w1 = weightsinterp_s2biquintic_ab_plus([ab])
    allvec_w1 = FmatFI2\(lfun.(tmatFI2));
    return dot(omall_w1,allvec_w1)
end
#### k=0, quadruple correction, interpolation
function w_k0_ptilde4(ab; lfun::Function=(x->ones(size(x))) )
    omall = weightsinterp_s2biquintic_ab_w4_plus([ab])
    allvec = FmatFI2\(lfun.(tmatFI2));
    return  [dot(omall[:,1,1],allvec); dot(omall[:,1,2],allvec); dot(omall[:,1,3],allvec); dot(omall[:,1,4],allvec)]
end

####################
#### "exact" limit computation, needed for debugging and checking accuracy
####################

#### k=0, single correction
function weights2D_k0_single_TOL(α,β; lfun::Function=(x->ones(size(x))), checking::Bool=false, TOL::Real=1e-5)
    # function weights2D_nonsmooth(α,β; θ=0,ϕ=0, checking::Bool=false, Integr::Real=0.0)
    α=Float64(α)
    β=Float64(β)
    Integr = 0.0
    a=Float64(1.9);
    order=8;

    s(x::Array{Float64,1},h::Float64) = (norm(x,Inf)>=abs(h*0.5-1e-12))/norm(x);
    # s(x::Array{Float64,1},h::Float64)=sqrt(x[1]^2+x[2]^2)^(-1)*(norm(x,Inf)>=abs(h*0.5-1e-12));
    g(x::Array{Float64,1})=exp(-norm(x)^order);
    ffun(x::Array{Float64,1},h)=s(x,h)*lfun(atan(x[2],x[1]))*g(x);
    Hold = 0.5;

    # Int1,_=quadgk(x->exp.(-x.^8),BigFloat(0.0),a; rtol=1e-12,order=order);
    # Int1 = big(9.417426998497014880874021440506391502415192067235218146676389256935818128440567e-01);
    Int1 = (9.4174269984970148e-01)
    # Int2,_=quadgk(x->lfun.(x),BigFloat(0.0),2*pi; rtol=1e-12);

    d = 1;
    Int2 = 0.0;
    m = 10;

    t=0:0.001:2*pi;
    lt = lfun.(t)
    m1 = maximum(lt)
    m2 = minimum(lt)
    while d>TOL
        # println(m)
        Int2old = Int2
        xtmp = (LinRange(0,2*pi,m+1))[1:m]
        htmp = xtmp[2]-xtmp[1]
        Int2 = htmp*sum(lfun.(xtmp))
        d = abs(Int2-Int2old)
        if checking
            println(d)
            figure(200);
            subplot(1,2,1)
            plot(zeros(2),[-Hold;Hold],"-k")
            plot([-Hold;Hold],zeros(2),"-k")
            subplot(1,2,2)
            plot(t,lt)
            plot(pi*0.25*ones(2),[m1;m2],"--k")
            plot(pi*0.5*ones(2),[m1;m2],"--k")
            plot(pi*ones(2),[m1;m2],"--k")
        end
        m *= 2
    end

    Integr = [Int1*Int2]

    Ivals1=[(0.0)];
    d1 = 1.0
    # TOL = 1e-12;
    l=1;
    # h=H[1]
    w0 = zeros(15);
    h = Hold;
    while d1>TOL
        M = Int(ceil(max(a/h-α,a/h-β)))
        nodes = [ [(i - α)*h ; (j - β)*h] for i in -M:M for j in [-M:-1;1:M] ]
        nodes = [ nodes; [ [(i - α)*h ; (j - β)*h] for i in [-M:-1;1:M] for j in [0] ] ]
        nodes = nodes[norm.(nodes,2).<a]
        Ivals1[1]=(h^2)*sum(ffun.(nodes,0.0));
        w0[l+1]=(Integr[1]-Ivals1[1])/(h*g(-[α*h;β*h]))# + s(h)*(α<1e-15);
        d1 = abs(w0[l]-w0[l+1])
        if checking
            println(d1," ",h*g(-[α*h;β*h]))
            # println()
            subplot(1,2,1)
            plot([α*h],[β*h],"*r")
        end
        l+=1;
        h*=0.5;
    end
    return 2*w0[l]-w0[l-1]
end

#### k=0, quadruple correction
function weights2D_k0_quadruple_TOL(α,β; lfun::Function=(x->ones(size(x))), checking::Bool=false, TOL::Real=1e-5)
    # function weights2D_nonsmooth(α,β; θ=0,ϕ=0, checking::Bool=false, Integr::Real=0.0)
    # α=Float64(α)
    # β=Float64(β)
    Integr = 0.0
    a=Float64(1.9);
    order=8;

    s(x::Array{Float64,1},h::Float64)=(norm(x,Inf)>=abs(h*0.5-1e-12))/norm(x);
    # s(x::Array{Float64,1},h::Float64)=sqrt(x[1]^2+x[2]^2)^(-1)*(norm(x,Inf)>=abs(h*0.5-1e-12));
    g(x::Array{Float64,1})=exp(-norm(x)^order);
    ffun(x::Array{Float64,1},h)=s(x,h)*lfun(atan(x[2],x[1]))*g(x);
    Hold = 0.5;

    Int11 = 9.417426998497014880874021440506391502415192067235218146676389256935818128440567e-01;
    Int12 = (4.53201238527738538991336076967878085039119309541285104185138983388879772057945e-01)
    Int13 = (2.963045230520751135425275146800503449978101982744337237045235381856304187181172e-01)
    # Int2,_=quadgk(x->lfun.(x),BigFloat(0.0),2*pi; rtol=1e-12);
    lfun21(x) = lfun(x)
    lfun22(x) = lfun(x)*cos(x)
    lfun23(x) = lfun(x)*sin(x)
    lfun24(x) = lfun(x)*sin(x)*cos(x)

    d = 1;
    Int2 = zeros(4);
    Int2old = zeros(4);
    m = 10;
    while d>TOL
        # println(m)
        Int2old .= Int2;
        xtmp = (LinRange(0,2*pi,m+1))[1:m]
        htmp = xtmp[2]-xtmp[1]
        Int2[1] = htmp*sum(lfun21.(xtmp))
        Int2[2] = htmp*sum(lfun22.(xtmp))
        Int2[3] = htmp*sum(lfun23.(xtmp))
        Int2[4] = htmp*sum(lfun24.(xtmp))
        d = norm(Int2-Int2old,Inf)
        # println(d)
        m *= 2
    end
    # println(Int2)

    # Integr = [2*Int1*Int2]
    Integr = [ Int11*Int2[1]; Int12*Int2[2]; Int12*Int2[3]; Int13*Int2[4] ]
    h = Hold
    # TOL = 1e-10;
    i = 1
    d = 1.0
    om = zeros(10,4);
    while d > TOL
        M = Int(round(a/h));
        X = [ [(j1 - α)*h ; (j2 - β)*h] for j1 in -M:M for j2 in [-M:-1; 2:M] ]
        X = [ X; [ [(j1 - α)*h ; (j2 - β)*h] for j1 in [-M:-1; 2:M] for j2 in [0;1] ] ]
        # X = X[norm.(X,2).<a]
        FFvec = ffun.(X,0.0)
        X1 = [ x[1] for x in X]
        X2 = [ x[2] for x in X]

        h2 = h^2
        #for i=1:1
        local T, Tx, Ty, Txy = h2*sum(FFvec), h2*sum(FFvec.*X1), h2*sum(FFvec.*X2), h2*sum(FFvec.*X1.*X2);

        # A = [g((α-1)*h) g((α)*h); g((α-1)*h)*(α-1)*h g((α)*h)*α*h];
        xij = [-α*h;-β*h]
        Xij = [xij,xij+[0;h],xij+[h;h],xij+[h;0]]

        A = zeros(4,4);
        for j = 1:4
            A[1,j] = g(Xij[j])
            A[2,j] = g(Xij[j])*Xij[j][1]
            A[3,j] = g(Xij[j])*Xij[j][2]
            A[4,j] = g(Xij[j])*Xij[j][1]*Xij[j][2]
            # A[5,j] = g(Xij[j])*Xij[j][1]^2
            # A[6,j] = g(Xij[j])*Xij[j][2]^2
        end
        b = [ Integr[1]-T; Integr[2]-Tx; Integr[3]-Ty; Integr[4]-Txy ]/h;
        om[i+1,1:4] = lu(A)\b;
        d = norm( [ om[i,1]-om[i+1,1]; om[i,2]-om[i+1,2]; om[i,3]-om[i+1,3]; om[i,4]-om[i+1,4] ] , Inf)
        if checking
            println(d)
        end
        i += 1
        h /= 2
        # println(abs(om[i,1]-om[i+1,1])," ",abs(om[i,2]-om[i+1,2]))
    end

    return [2*om[i,1]-om[i-1,1]; 2*om[i,2]-om[i-1,2]; 2*om[i,3]-om[i-1,3]; 2*om[i,4]-om[i-1,4]]
end

#### k=1, single correction # the weights are always 1
function weights2D_k1_single_TOL(α,β; lfun::Function=(x->1.0), checking::Bool=false, TOL::Real=1e-5)
    # function weights2D_nonsmooth(α,β; θ=0,ϕ=0, checking::Bool=false, Integr::Real=0.0)
    α=Float64(α)
    β=Float64(β)
    Integr = 0.0
    a=Float64(1.9);
    order=8;

    # s(x::Array{Float64,1},h::Float64)=sqrt(cos(θ)^2*(x[1]^2+x[2]^2)+sin(θ)^2*(sin(ϕ)*x[1]-cos(ϕ)*x[2])^2)^(-1)*(norm(x,Inf)>=abs(h*0.5-1e-12));
    # s(x::Array{Float64,1},h::Float64)=sqrt(x[1]^2+x[2]^2)^(-1)*(norm(x,Inf)>=abs(h*0.5-1e-12));
    g(x::Array{Float64,1})=exp(-norm(x)^order);
    ffun(x::Array{Float64,1},h)=lfun(atan(x[2],x[1]))*g(x)*(norm(x,Inf)>=abs(h*0.5-1e-12));
    Hold = 0.5;

    # Int1,_=quadgk(x->exp.(-x.^8),BigFloat(0.0),a; rtol=1e-12,order=order);
    # Int1 = big(9.417426998497014880874021440506391502415192067235218146676389256935818128440567e-01);
    # Int1 = (9.4174269984970148e-01)
    Int1 = (4.53201238527738538991336076967878085039119309541285104185138983388879772057945e-01)

    d = 1;
    Int2 = 0.0;
    m = 20;

    # t=0:0.001:2*pi;
    # lt = lfun.(t)
    # m1 = maximum(lt)
    # m2 = minimum(lt)
    while d>TOL
        # println(m)
        Int2old = Int2
        xtmp = (LinRange(0,2*pi,m+1))[1:m]
        htmp = xtmp[2]-xtmp[1]
        Int2 = htmp*sum(lfun.(xtmp))
        d = abs(Int2-Int2old)
        if checking
            println(d)
            # figure(200);
            # subplot(1,2,1)
            # plot(zeros(2),[-Hold;Hold],"-k")
            # plot([-Hold;Hold],zeros(2),"-k")
            # subplot(1,2,2)
            # plot(t,lt)
            # plot(pi*0.25*ones(2),[m1;m2],"--k")
            # plot(pi*0.5*ones(2),[m1;m2],"--k")
            # plot(pi*ones(2),[m1;m2],"--k")
        end
        m *= 2
    end

    Integr = [Int1*Int2]

    if checking
        println("Ref integr: $(Integr[1])")
        println("integral calculated, now onto weight")
    end

    Ivals1=[(0.0)];
    d1 = 1.0
    # TOL = 1e-12;
    l=1;
    # h=H[1]
    w0 = zeros(15);
    h = Hold;
    while d1>TOL
        M = Int(ceil(max(a/h-α,a/h-β)))
        nodes = [ [(i - α)*h ; (j - β)*h] for i in -M:M for j in [-M:-1;1:M] ]
        nodes = [ nodes; [ [(i - α)*h ; (j - β)*h] for i in [-M:-1;1:M] for j in [0] ] ]
        nodes = nodes[norm.(nodes,2).<a]
        Ivals1[1]=(h^2)*sum(ffun.(nodes,0.0));
        w0[l+1]=(Integr[1]-Ivals1[1])/(h^2*g(-[α*h;β*h]))# + s(h)*(α<1e-15);
        d1 = abs(w0[l]-w0[l+1])
        if checking
            println(d1," ",h^2*g(-[α*h;β*h]))
            # println()
            # subplot(1,2,1)
            # plot([α*h],[β*h],"*r")
        end
        l+=1;
        h*=0.5;
    end
    return 2*w0[l]-w0[l-1]
end

#### k=1, quadruple correction
function weights2D_k1_quadruple_TOL(α,β; lfun::Function=(x->ones(size(x))), checking::Bool=false, TOL::Real=1e-5)
    # function weights2D_nonsmooth(α,β; θ=0,ϕ=0, checking::Bool=false, Integr::Real=0.0)
    # α=Float64(α)
    # β=Float64(β)
    Integr = 0.0
    a = Float64(1.9);
    order=8;

    g(x::Array{Float64,1}) = exp(-norm(x)^order);
    ffun(x::Array{Float64,1},h) = lfun(atan(x[2],x[1]))*g(x)*(norm(x,Inf)>=abs(h*0.5-1e-12));
    Hold = 0.5;

    # Int11 = 9.417426998497014880874021440506391502415192067235218146676389256935818128440567e-01;
    Int11 = (4.53201238527738538991336076967878085039119309541285104185138983388879772057945e-01)
    Int12 = (2.963045230520751135425275146800503449978101982744337237045235381856304187181172e-01)
    Int13 = (2.2155673136318950341227093541764314784969368201529939626618235063916257459140e-01)
    # Int2,_=quadgk(x->lfun.(x),BigFloat(0.0),2*pi; rtol=1e-12);
    lfun21(x) = lfun(x)
    lfun22(x) = lfun(x)*cos(x)
    lfun23(x) = lfun(x)*sin(x)
    lfun24(x) = lfun(x)*sin(x)*cos(x)

    d = 1;
    Int2 = zeros(4);
    Int2old = zeros(4);
    m = 20;
    while d>TOL
        # println(m)
        Int2old .= Int2;
        xtmp = (LinRange(0,2*pi,m+1))[1:m]
        htmp = xtmp[2]-xtmp[1]
        Int2[1] = htmp*sum(lfun21.(xtmp))
        Int2[2] = htmp*sum(lfun22.(xtmp))
        Int2[3] = htmp*sum(lfun23.(xtmp))
        Int2[4] = htmp*sum(lfun24.(xtmp))
        d = norm(Int2-Int2old,Inf)
        # println(d)
        m *= 2
    end
    # println(Int2)

    # Integr = [2*Int1*Int2]
    Integr = [ Int11*Int2[1]; Int12*Int2[2]; Int12*Int2[3]; Int13*Int2[4] ]
    h = Hold
    # TOL = 1e-10;
    i = 1
    d = 1.0
    om = zeros(10,4);
    while d > TOL
        M = Int(round(a/h));
        X = [ [(j1 - α)*h ; (j2 - β)*h] for j1 in -M:M for j2 in [-M:-1; 2:M] ]
        X = [ X; [ [(j1 - α)*h ; (j2 - β)*h] for j1 in [-M:-1; 2:M] for j2 in [0;1] ] ]
        # X = X[norm.(X,2).<a]
        FFvec = ffun.(X,0.0)
        X1 = [ x[1] for x in X]
        X2 = [ x[2] for x in X]

        h2 = h^2
        #for i=1:1
        local T, Tx, Ty, Txy = h2*sum(FFvec), h2*sum(FFvec.*X1), h2*sum(FFvec.*X2), h2*sum(FFvec.*X1.*X2);

        # A = [g((α-1)*h) g((α)*h); g((α-1)*h)*(α-1)*h g((α)*h)*α*h];
        xij = [-α*h;-β*h]
        Xij = [xij,xij+[0;h],xij+[h;h],xij+[h;0]]

        A = zeros(4,4);
        for j = 1:4
            A[1,j] = g(Xij[j])
            A[2,j] = g(Xij[j])*Xij[j][1]
            A[3,j] = g(Xij[j])*Xij[j][2]
            A[4,j] = g(Xij[j])*Xij[j][1]*Xij[j][2]
            # A[5,j] = g(Xij[j])*Xij[j][1]^2
            # A[6,j] = g(Xij[j])*Xij[j][2]^2
        end
        b = [ Integr[1]-T; Integr[2]-Tx; Integr[3]-Ty; Integr[4]-Txy ]/h2;
        om[i+1,1:4] = lu(A)\b;
        d = norm( [ om[i,1]-om[i+1,1]; om[i,2]-om[i+1,2]; om[i,3]-om[i+1,3]; om[i,4]-om[i+1,4] ] , Inf)
        if checking
            println(d)
        end
        i += 1
        h /= 2
        # println(abs(om[i,1]-om[i+1,1])," ",abs(om[i,2]-om[i+1,2]))
    end

    return [2*om[i,1]-om[i-1,1]; 2*om[i,2]-om[i-1,2]; 2*om[i,3]-om[i-1,3]; 2*om[i,4]-om[i-1,4]]
end

#### k=2, single correction
function weights2D_k2_single_TOL(α,β; lfun::Function=(x->ones(size(x))), checking::Bool=false, TOL::Real=1e-5)
    # function weights2D_nonsmooth(α,β; θ=0,ϕ=0, checking::Bool=false, Integr::Real=0.0)
    α=Float64(α)
    β=Float64(β)
    Integr = 0.0
    # θ = θ%π
    # θ = abs(θ)*(abs(θ)<π*0.5) + (π-abs(θ))*(abs(θ)>π*0.5);
    # ϕ = ϕ%(2*pi)
    # ϕ = ϕ*(ϕ>=0) + (ϕ+2*π)*(ϕ<0);
    # ϕ = ϕ%pi
    # println("θ=$θ")
    a=Float64(1.9);
    order=8;

    # s(x::Array{Float64,1},h::Float64)=sqrt(cos(θ)^2*(x[1]^2+x[2]^2)+sin(θ)^2*(sin(ϕ)*x[1]-cos(ϕ)*x[2])^2)^(-1)*(norm(x,Inf)>=abs(h*0.5-1e-12));
    # s(x::Array{Float64,1},h::Float64)=sqrt(x[1]^2+x[2]^2)^(-1)*(norm(x,Inf)>=abs(h*0.5-1e-12));
    g(x::Array{Float64,1})=exp(-norm(x)^order);
    ffun(x::Array{Float64,1},h)=norm(x)*lfun(atan(x[2],x[1]))*g(x)*(norm(x,Inf)>=abs(h*0.5-1e-12));
    Hold = 0.5;

    # Int1,_=quadgk(x->exp.(-x.^8),BigFloat(0.0),a; rtol=1e-12,order=order);
    # Int1 = big(9.417426998497014880874021440506391502415192067235218146676389256935818128440567e-01);
    # Int1 = (9.4174269984970148e-01)
    # Int1 = (0.4532012385277385389)
    Int1 = (2.963045230520751135425275146800503449978101982744337237045235381856304187181172e-01)
    # Int2,_=quadgk(x->lfun.(x),BigFloat(0.0),2*pi; rtol=1e-12);
    # lfun2(x) = lfun(x-θ_0)/(sqrt(cos(θ)^2 + (sin(θ)^2)*(sin(ϕ-x)^2)))
    lfun2(x) = lfun(x)

    d = 1;
    Int2 = 0.0;
    m = 20;

    # t=0:0.001:2*pi;
    # lt = lfun.(t)
    # m1 = maximum(lt)
    # m2 = minimum(lt)
    while d>TOL
        # println(m)
        Int2old = Int2
        xtmp = (LinRange(0,2*pi,m+1))[1:m]
        htmp = xtmp[2]-xtmp[1]
        Int2 = htmp*sum(lfun2.(xtmp))
        d = abs(Int2-Int2old)
        if checking
            println(d)
            # figure(200);
            # subplot(1,2,1)
            # plot(zeros(2),[-Hold;Hold],"-k")
            # plot([-Hold;Hold],zeros(2),"-k")
            # subplot(1,2,2)
            # plot(t,lt)
            # plot(pi*0.25*ones(2),[m1;m2],"--k")
            # plot(pi*0.5*ones(2),[m1;m2],"--k")
            # plot(pi*ones(2),[m1;m2],"--k")
        end
        m *= 2
    end

    Integr = [Int1*Int2]

    if checking
        println("Ref integr: $(Integr[1])")
        println("integral calculated, now onto weight")
    end

    Ivals1=[(0.0)];
    d1 = 1.0
    # TOL = 1e-12;
    l=1;
    # h=H[1]
    w0 = zeros(15);
    h = Hold;
    h3 = h^3
    h2 = h^2
    while d1>TOL
        M = Int(ceil(max(a/h-α,a/h-β)))
        nodes = [ [(i - α)*h ; (j - β)*h] for i in -M:M for j in [-M:-1;1:M] ]
        nodes = [ nodes; [ [(i - α)*h ; (j - β)*h] for i in [-M:-1;1:M] for j in [0] ] ]
        nodes = nodes[norm.(nodes,2).<a]
        Ivals1[1]=(h2)*sum(ffun.(nodes,0.0));
        w0[l+1]=(Integr[1]-Ivals1[1])/(h3*g(-[α*h;β*h]))# + s(h)*(α<1e-15);
        d1 = abs(w0[l]-w0[l+1])
        if checking
            println(d1," ",h3*g(-[α*h;β*h]))
            # println()
            # subplot(1,2,1)
            # plot([α*h],[β*h],"*r")
        end
        l+=1;
        h*=0.5;
        h3 = h^3
        h2 = h^2
    end
    return 2*w0[l]-w0[l-1]
end


#################################################################
#################################################################
# CPM code for general surface. only info needed is projection mapping and inside/outside information
#################################################################
#################################################################

function genCPM_corr_V2_PB( Pgammafun::Function, insidepoint::Function, far_insidepoint::Function, far_outsidepoint::Function, X::Array{Float64,1}, Y::Array{Float64,1}, Z::Array{Float64,1}, Nvec::Array{Int,1}, epsl::Real, h::Real; outflag::Bool=true, surfTargetXYZ=[[1;0;0]], epsl_ratio::Real=1.0, kappa_val::Real=1.0 )

  secondt_const = 0.25/pi

  @time begin
    # inside or outside in the sense: exterior or interior domains are considered; nodes far enough away are considered, so no singularity or near singularity for solution evaluation

    insidefun(z) = true;
    if outflag
      insidefun = far_outsidepoint
    else
      insidefun = far_insidepoint
    end

    Nx, Ny, Nz = Nvec
    M_t = length(surfTargetXYZ)

    indIJK_to_M = -ones(Int64,Nx,Ny,Nz); # to go from ijk to m-index; initialized as -1
    Pg = Array{Float64,4}(undef,Nx,Ny,Nz,3) # 3D projection array for every index i,j,k,1:3 for each component
    dsignes = Array{Float64,3}(undef,Nx,Ny,Nz) # for every index ijk, signed distance (for curvature)
    Mtemp = Int(ceil(4.4*Nx^2.9/5)); # initialization size for how many interior points to the domain
    
    PΓ = Array{Array{Float64,1},1}(undef,Mtemp) # projection mapping for m-indices
    indM_to_IJK = zeros(Int64,Mtemp,3); # indices from m to ijk
    # d = Array{Float64,1}(undef,Mtemp) # unsigned distance, m-indices
    # ds = Array{Float64,1}(undef,Mtemp) # signed distance, m-indices
    
    m=0;
    for i=1:Nx
      for j=1:Ny
        for k=1:Nz # a->R2, b->R1
          zpt = [X[i];Y[j];Z[k]] # given node
          Pzpt = Pgammafun(zpt)
          Pg[i,j,k,:] = Pzpt # projection
          dsignes[i,j,k] = (1-2*insidepoint(zpt))*norm(Pzpt-zpt) # signed distance
          if ( norm(zpt-Pzpt) < epsl )
            m += 1; # inner node counter
            indIJK_to_M[i,j,k] = m; # node inside the tubular neighborhood
            indM_to_IJK[m,:] = [i;j;k];
            PΓ[m] = Pzpt
            # d[m] = norm(zpt-Pzpt)
            # ds[m] = dsignes[i,j,k]
          end
        end
      end
    end

    M = m; # finalized the correct number of interior nodes
    # cast to the correct size
    indM_to_IJK = indM_to_IJK[1:M,:]
    PΓ = PΓ[1:M]
    # d = d[1:M]
    # ds = ds[1:M]
  end

  # here starts the correction precomputation
  println("CPM mapping completed; interior M=$M.")

  # println( "Size of arrays to use (bytes):" )
  # println( "                           Pg: ",sizeof( Pg ) )
  # println( "                           PΓ: ",sizeof( PΓ ) )
  # println( "                  indIJK_to_M: ",sizeof( indIJK_to_M ) )
  # println( "                  indM_to_IJK: ",sizeof( indM_to_IJK ) )
  # println( "                      dsignes: ",sizeof( dsignes ) )

  begin
    M_int = Int(ceil(4*epsl/h))

    # _t is because it's for specifit targets, not for solving the general problem
    # for every target, an array of couples (α,β) is considered (initialized size M_int), corresponding to all the nodes to correct
    ab_single_t = [ zeros(Float64,2) for i=1:M_t, j=1:M_int ] # couples -1/2≦(α,β)≦1/2 for single correction
    iv1_t = [ -ones(Int64,M_int) for i=1:M_t ] # m-indices of the closest node needed for correction 
  
    w_K11_single_t = [ -ones(Float64,M_int) for i=1:M_t ] # for every target, weights for every correction needed; 1 weight because single correction
    w_K22_single_t = [ -ones(Float64,M_int) for i=1:M_t ] # for every target, weights for every correction needed; 1 weight because single correction
    w_K21_single_t = [ -ones(Float64,M_int) for i=1:M_t ] # for every target, weights for every correction needed; 1 weight because single correction
    # w_K12_single_t = [ -ones(Float64,M_int) for i=1:M_t ] # for every target, weights for every correction needed; 1 weight because single correction
      
    eta_vec = [ zeros(Float64,M_int) for i=1:M_t ] # only need it in this function

    target_normal = [ zeros(3) for i=1:M_t]
  end
  println("Initialization of correction arrays completed. Number of target points M_t=$M_t.")

  # println( "Size of arrays to use (bytes):" )
  # println( "                  ab_single_t: ",sizeof( ab_single_t ) )
  # println( "                        iv1_t: ",sizeof( iv1_t ) )
  # println( "               w_K11_single_t: ",sizeof( w_K11_single_t ) )
  # println( "               w_K22_single_t: ",sizeof( w_K22_single_t ) )
  # println( "               w_K21_single_t: ",sizeof( w_K21_single_t ) )
  # println( "               w_K12_single_t: ",sizeof( w_K12_single_t ) )

  @time begin
    for m = 1:M_t # one for every target
      x_t, y_t, z_t = surfTargetXYZ[m]
      zpt = [x_t;y_t;z_t]

      Pzpt = Pgammafun(zpt)
      sign_now = (1-2*insidepoint(zpt))
      dszpt = sign_now*norm(Pzpt-zpt) # signed distance

      # initialization of permuation of zpt, Pzpt
      Permzpt = zeros(3); Permzpt .= Pzpt

      #### calculating the curvatures and principal directions on the surface using the signed distance function (specifically its  Hessian, approximated with 2nd order in h)

      # build Hessian from scratch for target points:
      htmp = h
      Np = 4
      Mtmp2 = zeros(2Np+1,2Np+1,2Np+1)
      dsignes_loc = zeros(2Np+1,2Np+1,2Np+1)
      for itmp=1:2Np+1
        for jtmp=1:2Np+1
          for ktmp=1:2Np+1
            zpttmp = zpt .+ htmp*[itmp-Np-1;jtmp-Np-1;ktmp-Np-1] # discretization grid centered in zpt
            Pzpttmp = Pgammafun(zpttmp)
            dsignes_loc[itmp,jtmp,ktmp] = (1-2*insidepoint(zpttmp))*norm(Pzpttmp-zpttmp)
            Mtmp2[itmp,jtmp,ktmp] = (1-2*insidepoint(zpttmp))*norm(Pzpttmp-zpttmp) # signed distance
          end
        end
      end
      iind, i1 = get_indices_b(dsignes_loc[(Np-1):(Np+3),Np+1,Np+1]) # potentially one-sided 2nd order finite differences
      jind, j1 = get_indices_b(dsignes_loc[Np+1,(Np-1):(Np+3),Np+1])
      kind, k1 = get_indices_b(dsignes_loc[Np+1,Np+1,(Np-1):(Np+3)])
      norm_now = [dot(dsignes_loc[iind.+(Np+1),Np+1,Np+1],i1); dot(dsignes_loc[Np+1,jind.+(Np+1),Np+1],j1);dot(dsignes_loc[Np+1,Np+1,kind.+(Np+1)],k1)]/h
      norm_now /= norm(norm_now)
      target_normal[m] = norm_now # normal found via gradient of signed distance function

      Jm = hessian3D_2nd(Mtmp2[(Np):(Np+2),(Np):(Np+2),(Np):(Np+2)], htmp, htmp, htmp) # from 
      
      #### Attention: only works if target is one of the grid nodes.
      # i_now, j_now, k_now = indM_to_IJK[m,:] # only works if target is one of the grid nodes.
      # Jm = hessian3D( dsignes[i_now-2:i_now+2,j_now-2:j_now+2,k_now-2:k_now+2], dx, dy, dz);
      #### Attention end

      eigvalF, eigvecF = eigen(Jm);
      eigvp = sortperm(abs.(eigvalF))[3:-1:1]
      # find first principal direction
      tau1 = eigvecF[:,eigvp[1]]; tau1 /= norm(tau1) # normalization
      # find second principal direction by using the normal
      tau2 = cross(norm_now, tau1); tau2 /= norm(tau2)
      k1_now = -eigvalF[eigvp[1]] # first eigenvalue (from the formula it comes with opposite sign)
      k2_now = -eigvalF[eigvp[2]] # second eigenvalue, second largest. If it's smaller in absolute value than 0 eigenvalue corresponding to the normal direction, still fine, still zero.

      # k1_now and k2_now are the principal curvatures on the parallel surface, and tau1, tau2 are the principal directions so that the coordinate system on the surface is (tau1, tau2, norm_now)

      # transform the curvatures from the parallel surface to the original surface (zero level set)
      k1_now /= (1+dszpt*k1_now)
      k2_now /= (1+dszpt*k2_now)

      Mmat = [k1_now 0;0 k2_now] # M matrix from (3.25) in https://arxiv.org/abs/2203.04854
      D0(z) = [1/(1-k1_now*z) 0;0 1/(1-k2_now*z)] # matrix D_0 at the beginning of page 23

      Nmax = maximum(Nvec);
      W1 = zeros(Nmax); W2 = zeros(Nmax); W3 = zeros(Nmax)
      W1[1:Nx] .= X; W2[1:Ny] .= Y; W3[1:Nz] .= Z; # initialization: coord. system (x,y,z), no permutation
      
      theta = acos(round(norm_now[3]/norm(norm_now),digits=13)) # theta,phi directions of the normal, singularity line
      phi = atan( norm_now[2], norm_now[1] )
      sqrt2 = sqrt(2); pb = [1;2;3]; pf = [1;2;3]
      e1 = [1;0;0]; e2 = [0;1;0]; e3 = [0;0;1] # initialized coord system versors
      e1_now = [0;0;0]; e2_now = [0;0;0]; e3_now = [0;0;0];
      
      if ( (abs(tan(theta)) >= sqrt2 ) && (abs(tan(phi)))>1-1e-8 ) # new coordinate systems if normal is too tilted; angles are updated
        pb .= [2;3;1]; pf .= [3;1;2]
        theta = acos( norm_now[pf[3]]/norm(norm_now) );
        phi = ( abs(theta%(π))>1e-6 && abs(theta%(π))<(π-1e-6) )*atan( norm_now[pf[2]], norm_now[pf[1]] )
        W1[1:Nz] .= Z; W2[1:Nx] .= X; W3[1:Ny] .= Y;
      elseif ( (abs(tan(theta)) >= sqrt2 ) && ( abs(tan(phi))<1 ) )
        pb .= [3;1;2]; pf .= [2;3;1]
        theta = acos( norm_now[pf[3]]/norm(norm_now) );
        phi = (abs(theta%(π))>1e-6 && abs(theta%(π))<(π-1e-6))*atan( norm_now[pf[2]], norm_now[pf[1]] )
        W1[1:Ny] .= Y; W2[1:Nz] .= Z; W3[1:Nx] .= X;
      end
      Permzpt .= Permzpt[pf];
      Nvecp = Nvec[pf]
      
      e1_now .= e1[pb]; e2_now .= e2[pb]; e3_now .= e3[pb]
      e1_tmp = e1_now[pf]; e2_tmp = e2_now[pf]; nvec_tmp = norm_now[pf];
      
      Amat = [e1_tmp[1] e1_tmp[2];e2_tmp[1] e2_tmp[2]] # matrix A  and vector d from (3.21)
    #   dvec = [nvec_tmp[1]; nvec_tmp[2]] # not needed for single correction
      
      # a,b,c calculated with θn,ϕn angles of normal direction in 3D, because I need to calculate distance from a point to that line
      ac_angle = tan(theta)*cos(phi)
      bc_angle = tan(theta)*sin(phi) # values needed for permutated setting

      # range of indices to check, so that all the nodes in the tubular neighborhood are considered
      Mx3p = min( Nvecp[3], Int(  ceil((Permzpt[3] + epsl - W3[1])/h)))
      Mx3m = max( 1, Int( floor((Permzpt[3] - epsl - W3[1])/h)))
      
      indIm = Mx3m:Mx3p # range of indices in z-direction (after permutation)
      indIc = 1:length(indIm) # how many indices in total, from 1 to ...

      # intersection of line with planes in permutated basis
      Int_zpt = [ [ac_angle*( W3[im]-Permzpt[3])+Permzpt[1]; bc_angle*(W3[im]-Permzpt[3])+Permzpt[2]; W3[im]] for im in indIm ]

      # indices corresponding to (bottom left) closest grid node to the intersection
      Indint_p = [ [Int(floor((Int_zpt[im][1]-W1[1])/h))+1; Int(floor((Int_zpt[im][2]-W2[1])/h))+1; indIm[im] ] for im in indIc ]

      # α,β for bottom left closest grid node
      tmp1 = [ [( Int_zpt[im][1] - ((Indint_p[im][1]-1)*h + W1[1]) )/h; ( Int_zpt[im][2] - ((Indint_p[im][2]-1)*h + W2[1]) )/h] for im in indIc ]
      
      # indices (only first two) corresponding to (actual) closest point to  to intersection
      Indint_p_w1 = [ [Indint_p[im][1]+(tmp1[im][1]>0.5); Indint_p[im][2]+(tmp1[im][2]>0.5)] for im in indIc ]
    
      # α,β for closest grid node
      ab_single_t[m, indIc] = [ [( Int_zpt[im][1] - ((Indint_p_w1[im][1]-1)*h + W1[1]) )/h; ( Int_zpt[im][2] - ((Indint_p_w1[im][2]-1)*h + W2[1]) )/h] for im in indIc ]

      # 3D indices of single correction nodes, correctly permutated to standard xyz system (pb, backward)
      Ind_tmp = [ ([ Indint_p_w1[im][1]; Indint_p_w1[im][2]; Indint_p[im][3] ])[pb] for im in indIc ]
      # third is from Indint_p, because the third is not saved into Indint_p_w1

      # 1D array of subindices of indIc for which the indices of the corr. nodes are not out of bounds, important for not passing bad values to indIJK_to_M
      ind_check = indIc[ [ prod( (ind_tmp.<=Nvec) .& (ind_tmp.>=1) ) for ind_tmp in Ind_tmp ].>0 ]

      # 1D array of m-indices of valid correction points; indices given to indIJK_to_M to get m-index
      (iv1_t[m])[ind_check] = [ indIJK_to_M[ind_tmp[1],ind_tmp[2],ind_tmp[3]] for ind_tmp in Ind_tmp[ind_check] ]

      for itmp in ind_check # indices sets which are not out of bounds
          η = (1-2*insidepoint((Int_zpt[itmp])[pb]))*norm(Permzpt-Int_zpt[itmp]) 
          (eta_vec[m])[itmp] = η
          α, β = ab_single_t[m, itmp]
          D0_now = D0(η)

          # weights for the different kernels
          w_K11_single_t[ m ][ itmp ] = w_k0_ptilde1([α;β]; lfun=(t->i0_PB_a(t,D0_now,Amat, Mmat, epsl_ratio)))*secondt_const
          w_K22_single_t[ m ][ itmp ] = w_k0_ptilde1([α;β]; lfun=(t->j0_PB_a(t,D0_now,Amat, Mmat, 1/epsl_ratio)))*(-secondt_const)
          w_K21_single_t[ m ][ itmp ] = w_k0_ptilde1([α;β]; lfun=(t->1.0./ψ0_PB_a(t,D0_now,Amat)))*(secondt_const)*0.5*kappa_val^2
          # w_K12_single_t[ m ][ itmp ] = secondt_const*kappa_val # weight is always one
          # technically useless array
      end
      
      # only keep the values corresponding to nodes inside the neighborhod
      p1tmp = (iv1_t[m]) .> 0
      iv1_t[m] = (iv1_t[m])[p1tmp]
      
      eta_vec[ m ] = (eta_vec[ m ])[p1tmp]     
      
      w_K11_single_t[ m ] = (w_K11_single_t[ m ])[p1tmp]
      w_K22_single_t[ m ] = (w_K22_single_t[ m ])[p1tmp]
      # w_K12_single_t[ m ] = (w_K12_single_t[ m ])[p1tmp]
      w_K21_single_t[ m ] = (w_K21_single_t[ m ])[p1tmp]

    end
  end
  println("Correction arrays filled - targets")
  println("Only targets considered")

  println("--- Precomputations over ---")

  return M, Pg, PΓ, indIJK_to_M, indM_to_IJK, dsignes, ab_single_t,
  iv1_t, w_K11_single_t, w_K22_single_t, w_K21_single_t, target_normal
end

function genCPM_corr_PB_system( Pgammafun::Function, insidepoint::Function, X::Array{Float64,1}, Y::Array{Float64,1}, Z::Array{Float64,1}, Nvec::Array{Int,1}, epsl::Real, h::Real; outflag::Bool=true, surfTargetXYZ=[[1;0;0]], epsl_ratio::Real=1.0, kappa_val::Real=1.0 )

  secondt_const = 0.25/pi

  @time begin
    # inside or outside in the sense: exterior or interior domains are considered; nodes far enough away are considered, so no singularity or near singularity for solution evaluation

    # insidefun(z) = true;
    # if outflag
    #   insidefun = far_outsidepoint
    # else
    #   insidefun = far_insidepoint
    # end

    Nx, Ny, Nz = Nvec
    M_t = length(surfTargetXYZ)

    indIJK_to_M = -ones(Int64,Nx,Ny,Nz); # to go from ijk to m-index; initialized as -1
    Pg = Array{Float64,4}(undef,Nx,Ny,Nz,3) # 3D projection array for every index i,j,k,1:3 for each component
    dsignes = Array{Float64,3}(undef,Nx,Ny,Nz) # for every index ijk, signed distance (for curvature)
    Mtemp = Int(ceil(4.4*Nx^2.9/3)); # initialization size for how many interior points to the domain
    
    PΓ = Array{Array{Float64,1},1}(undef,Mtemp) # projection mapping for m-indices
    indM_to_IJK = zeros(Int64,Mtemp,3); # indices from m to ijk
    # d = Array{Float64,1}(undef,Mtemp) # unsigned distance, m-indices
    # ds = Array{Float64,1}(undef,Mtemp) # signed distance, m-indices
    
    m=0;
    for i=1:Nx
      for j=1:Ny
        for k=1:Nz # a->R2, b->R1
          zpt = [X[i];Y[j];Z[k]] # given node
          Pzpt = Pgammafun(zpt)
          Pg[i,j,k,:] = Pzpt # projection
          dsignes[i,j,k] = (1-2*insidepoint(zpt))*norm(Pzpt-zpt) # signed distance
          if ( norm(zpt-Pzpt) < epsl )
            m += 1; # inner node counter
            indIJK_to_M[i,j,k] = m; # node inside the tubular neighborhood
            indM_to_IJK[m,:] = [i;j;k];
            PΓ[m] = Pzpt
            # d[m] = norm(zpt-Pzpt)
            # ds[m] = dsignes[i,j,k]
          end
        end
      end
    end

    M = m; # finalized the correct number of interior nodes
    # cast to the correct size
    indM_to_IJK = indM_to_IJK[1:M,:]
    PΓ = PΓ[1:M]
    # d = d[1:M]
    # ds = ds[1:M]
  end

  # here starts the correction precomputation
  println("CPM mapping completed; interior M=$M.")

  # println( "Size of arrays to use (bytes):" )
  # println( "                           Pg: ",sizeof( Pg ) )
  # println( "                           PΓ: ",sizeof( PΓ ) )
  # println( "                  indIJK_to_M: ",sizeof( indIJK_to_M ) )
  # println( "                  indM_to_IJK: ",sizeof( indM_to_IJK ) )
  # println( "                      dsignes: ",sizeof( dsignes ) )

  begin
    M_int = Int(ceil(4*epsl/h))

    # _t is because it's for specifit targets, not for solving the general problem
    # for every target, an array of couples (α,β) is considered (initialized size M_int), corresponding to all the nodes to correct
    ab_single = [ zeros(Float64,2) for i=1:M, j=1:M_int ] # couples -1/2≦(α,β)≦1/2 for single correction
    iv1 = [ -ones(Int64,M_int) for i=1:M ] # m-indices of the closest node needed for correction 
  
    w_K11_single = [ -ones(Float64,M_int) for i=1:M ] # for every target, weights for every correction needed; 1 weight because single correction
    w_K22_single = [ -ones(Float64,M_int) for i=1:M ] # for every target, weights for every correction needed; 1 weight because single correction
    w_K21_single = [ -ones(Float64,M_int) for i=1:M ] # for every target, weights for every correction needed; 1 weight because single correction
    # w_K12_single_t = [ -ones(Float64,M_int) for i=1:M ] # for every target, weights for every correction needed; 1 weight because single correction
      
    eta_vec = [ zeros(Float64,M_int) for i=1:M ] # only need it in this function

    normal = [ zeros(3) for i=1:M]
  end
  println("Initialization of correction arrays completed. Number of nodes M=$M.")

  # println( "Size of arrays to use (bytes):" )
  # println( "                  ab_single_t: ",sizeof( ab_single_t ) )
  # println( "                        iv1: ",sizeof( iv1 ) )
  # println( "               w_K11_single: ",sizeof( w_K11_single ) )
  # println( "               w_K22_single: ",sizeof( w_K22_single ) )
  # println( "               w_K21_single: ",sizeof( w_K21_single ) )
  # println( "               w_K12_single_t: ",sizeof( w_K12_single_t ) )

  @time begin
    for m = 1:M # one for every node
      x_t, y_t, z_t = PΓ[m]
      zpt = [x_t;y_t;z_t]

      i0,j0,k0 = indM_to_IJK[m,:]
      zptP = [X[i0];Y[j0];Z[k0]]

      Pzpt = zpt
      sign_now = (1-2*insidepoint(zpt))
      dszpt = sign_now*norm(Pzpt-zpt) # signed distance

      # initialization of permuation of zpt, Pzpt
      Permzpt = zeros(3); Permzpt .= Pzpt

      #### calculating the curvatures and principal directions on the surface using the signed distance function (specifically its  Hessian, approximated with 2nd order in h)

      # build Hessian from scratch for target points:
      # htmp = h
      # Np = 4
      # Mtmp2 = zeros(2Np+1,2Np+1,2Np+1)
      # dsignes_loc = zeros(2Np+1,2Np+1,2Np+1)
      # for itmp=1:2Np+1
      #   for jtmp=1:2Np+1
      #     for ktmp=1:2Np+1
      #       zpttmp = zpt .+ htmp*[itmp-Np-1;jtmp-Np-1;ktmp-Np-1] # discretization grid centered in zpt
      #       Pzpttmp = Pgammafun(zpttmp)
      #       dsignes_loc[itmp,jtmp,ktmp] = (1-2*insidepoint(zpttmp))*norm(Pzpttmp-zpttmp)
      #       Mtmp2[itmp,jtmp,ktmp] = (1-2*insidepoint(zpttmp))*norm(Pzpttmp-zpttmp) # signed distance
      #     end
      #   end
      # end
      
      # iind, i1 = get_indices_b(dsignes_loc[(Np-1):(Np+3),Np+1,Np+1]) # potentially one-sided 2nd order finite differences
      # jind, j1 = get_indices_b(dsignes_loc[Np+1,(Np-1):(Np+3),Np+1])
      # kind, k1 = get_indices_b(dsignes_loc[Np+1,Np+1,(Np-1):(Np+3)])
      # norm_now = [dot(dsignes_loc[iind.+(Np+1),Np+1,Np+1],i1); dot(dsignes_loc[Np+1,jind.+(Np+1),Np+1],j1);dot(dsignes_loc[Np+1,Np+1,kind.+(Np+1)],k1)]/h
      # norm_now /= norm(norm_now)
      # normal[m] = norm_now # normal found via gradient of signed distance function

      # Jm = hessian3D_2nd(Mtmp2[(Np):(Np+2),(Np):(Np+2),(Np):(Np+2)], htmp, htmp, htmp) # from 
      
      #### Attention: only works if target is one of the grid nodes.
      i_now, j_now, k_now = indM_to_IJK[m,:] # only works if target is one of the grid nodes.
      Jm = hessian3D_2nd( dsignes[i_now-1:i_now+1, j_now-1:j_now+1, k_now-1:k_now+1], h, h, h);

      iind, i1 = get_indices_b(dsignes[i_now-2:i_now+2,j_now,k_now]) # potentially one-sided 2nd order finite differences
      jind, j1 = get_indices_b(dsignes[i_now,j_now-2:j_now+2,k_now])
      kind, k1 = get_indices_b(dsignes[i_now,j_now,k_now-2:k_now+2])
      norm_now = [dot(dsignes[i_now.+iind,j_now,k_now],i1); dot(dsignes[i_now,j_now.+jind,k_now],j1);dot(dsignes[i_now,j_now,k_now.+kind],k1)]/h
      norm_now /= norm(norm_now)
      normal[m] = norm_now # normal found via gradient of signed distance function
      #### Attention end

      eigvalF, eigvecF = eigen(Jm);
      eigvp = sortperm(abs.(eigvalF))[3:-1:1]
      # find first principal direction
      tau1 = eigvecF[:,eigvp[1]]; tau1 /= norm(tau1) # normalization
      # find second principal direction by using the normal
      tau2 = cross(norm_now, tau1); tau2 /= norm(tau2)
      k1_now = -eigvalF[eigvp[1]] # first eigenvalue (from the formula it comes with opposite sign)
      k2_now = -eigvalF[eigvp[2]] # second eigenvalue, second largest. If it's smaller in absolute value than 0 eigenvalue corresponding to the normal direction, still fine, still zero.

      # k1_now and k2_now are the principal curvatures on the parallel surface, and tau1, tau2 are the principal directions so that the coordinate system on the surface is (tau1, tau2, norm_now)

      # transform the curvatures from the parallel surface to the original surface (zero level set)
      k1_now /= (1+dszpt*k1_now)
      k2_now /= (1+dszpt*k2_now)

      Mmat = [k1_now 0;0 k2_now] # M matrix from (3.25) in https://arxiv.org/abs/2203.04854
      D0(z) = [1/(1-k1_now*z) 0;0 1/(1-k2_now*z)] # matrix D_0 at the beginning of page 23

      Nmax = maximum(Nvec);
      W1 = zeros(Nmax); W2 = zeros(Nmax); W3 = zeros(Nmax)
      W1[1:Nx] .= X; W2[1:Ny] .= Y; W3[1:Nz] .= Z; # initialization: coord. system (x,y,z), no permutation
      
      theta = acos(round(norm_now[3]/norm(norm_now),digits=13)) # theta,phi directions of the normal, singularity line
      phi = atan( norm_now[2], norm_now[1] )
      sqrt2 = sqrt(2); pb = [1;2;3]; pf = [1;2;3]
      e1 = [1;0;0]; e2 = [0;1;0]; e3 = [0;0;1] # initialized coord system versors
      e1_now = [0;0;0]; e2_now = [0;0;0]; e3_now = [0;0;0];
      
      if ( (abs(tan(theta)) >= sqrt2 ) && (abs(tan(phi)))>1-1e-8 ) # new coordinate systems if normal is too tilted; angles are updated
        pb .= [2;3;1]; pf .= [3;1;2]
        theta = acos( norm_now[pf[3]]/norm(norm_now) );
        phi = ( abs(theta%(π))>1e-6 && abs(theta%(π))<(π-1e-6) )*atan( norm_now[pf[2]], norm_now[pf[1]] )
        W1[1:Nz] .= Z; W2[1:Nx] .= X; W3[1:Ny] .= Y;
      elseif ( (abs(tan(theta)) >= sqrt2 ) && ( abs(tan(phi))<1 ) )
        pb .= [3;1;2]; pf .= [2;3;1]
        theta = acos( norm_now[pf[3]]/norm(norm_now) );
        phi = (abs(theta%(π))>1e-6 && abs(theta%(π))<(π-1e-6))*atan( norm_now[pf[2]], norm_now[pf[1]] )
        W1[1:Ny] .= Y; W2[1:Nz] .= Z; W3[1:Nx] .= X;
      end
      Permzpt .= Permzpt[pf];
      zptPerm = zptP[pf];
      Nvecp = Nvec[pf]
      
      e1_now .= e1[pb]; e2_now .= e2[pb]; e3_now .= e3[pb]
      t1_tmp = tau1[pf]; t2_tmp = tau2[pf]; nvec_tmp = norm_now[pf];
      
      Amat = [t1_tmp[1] t1_tmp[2];t2_tmp[1] t2_tmp[2]] # matrix A  and vector d from (3.21)
      #   dvec = [nvec_tmp[1]; nvec_tmp[2]] # not needed for single correction
      
      # a,b,c calculated with θn,ϕn angles of normal direction in 3D, because I need to calculate distance from a point to that line
      ac_angle = tan(theta)*cos(phi)
      bc_angle = tan(theta)*sin(phi) # values needed for permutated setting

      # range of indices to check, so that all the nodes in the tubular neighborhood are considered
      Mx3p = min( Nvecp[3], Int(  ceil((Permzpt[3] + epsl - W3[1])/h)))
      Mx3m = max( 1, Int( floor((Permzpt[3] - epsl - W3[1])/h)))
      # Mx0 = Int(  round((zptPerm[3] - W3[1])/h))
      Mx0 = ([i0;j0;k0])[pf[3]]
      # if ((Mx0>Mx3p) && (Mx0<Mx3m))
        # println("Mx0=$Mx0, between $(Mx3m)-$(Mx3p)")
      # end

      indIm = setdiff(Mx3m:Mx3p, Mx0) # range of indices in z-direction (after permutation)
      indIc = 1:length(indIm) # how many indices in total, from 1 to ...

      # intersection of line with planes in permutated basis
      Int_zpt = [ [ac_angle*( W3[im]-Permzpt[3])+Permzpt[1]; bc_angle*(W3[im]-Permzpt[3])+Permzpt[2]; W3[im]] for im in indIm ]
      # Int_zpt = [Int_zpt; [zptPerm] ]

      # indices corresponding to (bottom left) closest grid node to the intersection
      Indint_p = [ [Int(floor((Int_zpt[im][1]-W1[1])/h))+1; Int(floor((Int_zpt[im][2]-W2[1])/h))+1; indIm[im] ] for im in indIc ]

      # α,β for bottom left closest grid node
      tmp1 = [ [( Int_zpt[im][1] - ((Indint_p[im][1]-1)*h + W1[1]) )/h; ( Int_zpt[im][2] - ((Indint_p[im][2]-1)*h + W2[1]) )/h] for im in indIc ]
      
      # indices (only first two) corresponding to (actual) closest point to  to intersection
      Indint_p_w1 = [ [Indint_p[im][1]+(tmp1[im][1]>0.5); Indint_p[im][2]+(tmp1[im][2]>0.5)] for im in indIc ]
    
      # α,β for closest grid node
      ab_single[m, indIc] = [ [( Int_zpt[im][1] - ((Indint_p_w1[im][1]-1)*h + W1[1]) )/h; ( Int_zpt[im][2] - ((Indint_p_w1[im][2]-1)*h + W2[1]) )/h] for im in indIc ]
      # ab_single[m, length(indIc)+1] = [0;0.]

      # 3D indices of single correction nodes, correctly permutated to standard xyz system (pb, backward)
      Ind_tmp = [ ([ Indint_p_w1[im][1]; Indint_p_w1[im][2]; Indint_p[im][3] ])[pb] for im in indIc ]
      # third is from Indint_p, because the third is not saved into Indint_p_w1

      # 1D array of subindices of indIc for which the indices of the corr. nodes are not out of bounds, important for not passing bad values to indIJK_to_M
      ind_check = indIc[ [ prod( (ind_tmp.<=Nvec) .& (ind_tmp.>=1) ) for ind_tmp in Ind_tmp ].>0 ]

      # 1D array of m-indices of valid correction points; indices given to indIJK_to_M to get m-index
      (iv1[m])[ind_check] = [ indIJK_to_M[ind_tmp[1],ind_tmp[2],ind_tmp[3]] for ind_tmp in Ind_tmp[ind_check] ]
      (iv1[m])[length(indIc)+1] = m

      for itmp in ind_check # indices sets which are not out of bounds
          η = (1-2*insidepoint((Int_zpt[itmp])[pb]))*norm(Permzpt-Int_zpt[itmp]) 
          (eta_vec[m])[itmp] = η
          α, β = ab_single[m, itmp]
          D0_now = D0(η)

          # weights for the different kernels
          # w_K11_single[ m ][ itmp ] = w_k0_ptilde1([α;β]; lfun=(t->i0_PB_a(t,D0_now,Amat, Mmat, epsl_ratio)))*secondt_const
          # w_K22_single[ m ][ itmp ] = w_k0_ptilde1([α;β]; lfun=(t->j0_PB_a(t,D0_now,Amat, Mmat, 1/epsl_ratio)))*(-secondt_const)
          # w_K21_single[ m ][ itmp ] = w_k0_ptilde1([α;β]; lfun=(t->1.0./ψ0_PB_a(t,D0_now,Amat)))*(secondt_const)*0.5*kappa_val^2

          wjj_test = w_k0_ptilde1([α;β]; lfun=(t->ij0_PB_a(t,D0_now,Amat, Mmat)));
          # w11_test = wjj_test*secondt_const*(1-epsl_ratio)
          # w22_test = wjj_test*secondt_const*(1- 1/epsl_ratio)
          # w21_test = w_k0_ptilde1([α;β]; lfun=(t->1/norm(D0_now*Amat*[cos(t);sin(t)])))*secondt_const*0.5*kappa_val^2

          w_K11_single[ m ][ itmp ] = wjj_test*secondt_const*(1-epsl_ratio)
          w_K22_single[ m ][ itmp ] = wjj_test*secondt_const*(1- 1/epsl_ratio)
          w_K21_single[ m ][ itmp ] = w_k0_ptilde1([α;β]; lfun=(t->1/norm(D0_now*Amat*[cos(t);sin(t)])))*secondt_const*0.5*kappa_val^2
          # errval1 = [abs(w11_test-w_K11_single[ m ][ itmp ]);abs(w22_test-w_K22_single[ m ][ itmp ]);abs(w21_test-w_K21_single[ m ][ itmp ])]
          # if ~prod(errval1 .< 1e-9)
          #   println("weight errors")
          #   println("0=",abs(w11_test-w_K11_single[ m ][ itmp ]))
          #   println("0=",abs(w22_test-w_K22_single[ m ][ itmp ]))
          #   println("0=",abs(w21_test-w_K21_single[ m ][ itmp ]))  
          # end
          # w_K12_single_t[ m ][ itmp ] = secondt_const*kappa_val # weight is always one
          # technically useless array
          # ij0_PB_a(x,D0,Amat,Mmat)
      end
      itmp = length(indIc)+1
      # (eta_vec[m])[itmp] = dsignes[i0,j0,k0]
      η = dsignes[i0,j0,k0]
      α = 0.; β = 0.;
      D0_now = D0(η)
      # weights for the different kernels
      w_K11_single[ m ][ itmp ] = w_k0_ptilde1([α;β]; lfun=(t->i0_PB_a(t,D0_now,Amat, Mmat, epsl_ratio)))*secondt_const
      w_K22_single[ m ][ itmp ] = w_k0_ptilde1([α;β]; lfun=(t->j0_PB_a(t,D0_now,Amat, Mmat, 1/epsl_ratio)))*(-secondt_const)
      w_K21_single[ m ][ itmp ] = w_k0_ptilde1([α;β]; lfun=(t->1.0./ψ0_PB_a(t,D0_now,Amat)))*(secondt_const)*0.5*kappa_val^2

      
      # only keep the values corresponding to nodes inside the neighborhod
      p1tmp = (iv1[m]) .> 0
      iv1[m] = (iv1[m])[p1tmp]

      # iv2 = (iv1[m])[1:end-1]
      # for i=1:length(iv2)
      #   if (iv2[i]-m==0)
      #     println("Mx0=$Mx0, between $(Mx3m)-$(Mx3p)")
      #     println(indIm)
      #     println("m=$m, indices $(iv1[m])")
      #   end
      # end
      
      eta_vec[ m ] = (eta_vec[ m ])[p1tmp]     
      
      w_K11_single[ m ] = (w_K11_single[ m ])[p1tmp]
      w_K22_single[ m ] = (w_K22_single[ m ])[p1tmp]
      # w_K12_single_t[ m ] = (w_K12_single_t[ m ])[p1tmp]
      w_K21_single[ m ] = (w_K21_single[ m ])[p1tmp]

    end
  end
  println("Correction arrays filled - nodes")

  # setup for evaluation in specific target points on the surface
  begin
    M_int = Int(ceil(4*epsl/h))

    # _t is because it's for specifit targets, not for solving the general problem
    # for every target, an array of couples (α,β) is considered (initialized size M_int), corresponding to all the nodes to correct
    ab_single_t = [ zeros(Float64,2) for i=1:M_t, j=1:M_int ] # couples -1/2≦(α,β)≦1/2 for single correction
    iv1_t = [ -ones(Int64,M_int) for i=1:M_t ] # m-indices of the closest node needed for correction 
  
    w_K11_single_t = [ -ones(Float64,M_int) for i=1:M_t ] # for every target, weights for every correction needed; 1 weight because single correction
    w_K22_single_t = [ -ones(Float64,M_int) for i=1:M_t ] # for every target, weights for every correction needed; 1 weight because single correction
    w_K21_single_t = [ -ones(Float64,M_int) for i=1:M_t ] # for every target, weights for every correction needed; 1 weight because single correction
    # w_K12_single_t = [ -ones(Float64,M_int) for i=1:M_t ] # for every target, weights for every correction needed; 1 weight because single correction
      
    eta_vec = [ zeros(Float64,M_int) for i=1:M_t ] # only need it in this function

    target_normal = [ zeros(3) for i=1:M_t]
  end
  println("Initialization of correction arrays completed. Number of target points M_t=$M_t.")

  # println( "Size of arrays to use (bytes):" )
  # println( "                  ab_single_t: ",sizeof( ab_single_t ) )
  # println( "                        iv1_t: ",sizeof( iv1_t ) )
  # println( "               w_K11_single_t: ",sizeof( w_K11_single_t ) )
  # println( "               w_K22_single_t: ",sizeof( w_K22_single_t ) )
  # println( "               w_K21_single_t: ",sizeof( w_K21_single_t ) )
  # println( "               w_K12_single_t: ",sizeof( w_K12_single_t ) )

  @time begin
    for m = 1:M_t # one for every target
      x_t, y_t, z_t = surfTargetXYZ[m]
      zpt = [x_t;y_t;z_t]

      Pzpt = Pgammafun(zpt)
      sign_now = (1-2*insidepoint(zpt))
      dszpt = sign_now*norm(Pzpt-zpt) # signed distance

      # initialization of permuation of zpt, Pzpt
      Permzpt = zeros(3); Permzpt .= Pzpt

      #### calculating the curvatures and principal directions on the surface using the signed distance function (specifically its  Hessian, approximated with 2nd order in h)

      # build Hessian from scratch for target points:
      htmp = h
      Np = 4
      Mtmp2 = zeros(2Np+1,2Np+1,2Np+1)
      dsignes_loc = zeros(2Np+1,2Np+1,2Np+1)
      for itmp=1:2Np+1
        for jtmp=1:2Np+1
          for ktmp=1:2Np+1
            zpttmp = zpt .+ htmp*[itmp-Np-1;jtmp-Np-1;ktmp-Np-1] # discretization grid centered in zpt
            Pzpttmp = Pgammafun(zpttmp)
            dsignes_loc[itmp,jtmp,ktmp] = (1-2*insidepoint(zpttmp))*norm(Pzpttmp-zpttmp)
            Mtmp2[itmp,jtmp,ktmp] = (1-2*insidepoint(zpttmp))*norm(Pzpttmp-zpttmp) # signed distance
          end
        end
      end
      iind, i1 = get_indices_b(dsignes_loc[(Np-1):(Np+3),Np+1,Np+1]) # potentially one-sided 2nd order finite differences
      jind, j1 = get_indices_b(dsignes_loc[Np+1,(Np-1):(Np+3),Np+1])
      kind, k1 = get_indices_b(dsignes_loc[Np+1,Np+1,(Np-1):(Np+3)])
      norm_now = [dot(dsignes_loc[iind.+(Np+1),Np+1,Np+1],i1); dot(dsignes_loc[Np+1,jind.+(Np+1),Np+1],j1);dot(dsignes_loc[Np+1,Np+1,kind.+(Np+1)],k1)]/h
      norm_now /= norm(norm_now)
      target_normal[m] = norm_now # normal found via gradient of signed distance function

      Jm = hessian3D_2nd(Mtmp2[(Np):(Np+2),(Np):(Np+2),(Np):(Np+2)], htmp, htmp, htmp) # from 
      
      #### Attention: only works if target is one of the grid nodes.
      # i_now, j_now, k_now = indM_to_IJK[m,:] # only works if target is one of the grid nodes.
      # Jm = hessian3D( dsignes[i_now-2:i_now+2,j_now-2:j_now+2,k_now-2:k_now+2], dx, dy, dz);
      #### Attention end

      eigvalF, eigvecF = eigen(Jm);
      eigvp = sortperm(abs.(eigvalF))[3:-1:1]
      # find first principal direction
      tau1 = eigvecF[:,eigvp[1]]; tau1 /= norm(tau1) # normalization
      # find second principal direction by using the normal
      tau2 = cross(norm_now, tau1); tau2 /= norm(tau2)
      k1_now = -eigvalF[eigvp[1]] # first eigenvalue (from the formula it comes with opposite sign)
      k2_now = -eigvalF[eigvp[2]] # second eigenvalue, second largest. If it's smaller in absolute value than 0 eigenvalue corresponding to the normal direction, still fine, still zero.

      # k1_now and k2_now are the principal curvatures on the parallel surface, and tau1, tau2 are the principal directions so that the coordinate system on the surface is (tau1, tau2, norm_now)

      # transform the curvatures from the parallel surface to the original surface (zero level set)
      k1_now /= (1+dszpt*k1_now)
      k2_now /= (1+dszpt*k2_now)

      Mmat = [k1_now 0;0 k2_now] # M matrix from (3.25) in https://arxiv.org/abs/2203.04854
      D0(z) = [1/(1-k1_now*z) 0;0 1/(1-k2_now*z)] # matrix D_0 at the beginning of page 23

      Nmax = maximum(Nvec);
      W1 = zeros(Nmax); W2 = zeros(Nmax); W3 = zeros(Nmax)
      W1[1:Nx] .= X; W2[1:Ny] .= Y; W3[1:Nz] .= Z; # initialization: coord. system (x,y,z), no permutation
      
      theta = acos(round(norm_now[3]/norm(norm_now),digits=13)) # theta,phi directions of the normal, singularity line
      phi = atan( norm_now[2], norm_now[1] )
      sqrt2 = sqrt(2); pb = [1;2;3]; pf = [1;2;3]
      e1 = [1;0;0]; e2 = [0;1;0]; e3 = [0;0;1] # initialized coord system versors
      e1_now = [0;0;0]; e2_now = [0;0;0]; e3_now = [0;0;0];
      
      if ( (abs(tan(theta)) >= sqrt2 ) && (abs(tan(phi)))>1-1e-8 ) # new coordinate systems if normal is too tilted; angles are updated
        pb .= [2;3;1]; pf .= [3;1;2]
        theta = acos( norm_now[pf[3]]/norm(norm_now) );
        phi = ( abs(theta%(π))>1e-6 && abs(theta%(π))<(π-1e-6) )*atan( norm_now[pf[2]], norm_now[pf[1]] )
        W1[1:Nz] .= Z; W2[1:Nx] .= X; W3[1:Ny] .= Y;
      elseif ( (abs(tan(theta)) >= sqrt2 ) && ( abs(tan(phi))<1 ) )
        pb .= [3;1;2]; pf .= [2;3;1]
        theta = acos( norm_now[pf[3]]/norm(norm_now) );
        phi = (abs(theta%(π))>1e-6 && abs(theta%(π))<(π-1e-6))*atan( norm_now[pf[2]], norm_now[pf[1]] )
        W1[1:Ny] .= Y; W2[1:Nz] .= Z; W3[1:Nx] .= X;
      end
      Permzpt .= Permzpt[pf];
      Nvecp = Nvec[pf]
      
      e1_now .= e1[pb]; e2_now .= e2[pb]; e3_now .= e3[pb]
      t1_tmp = tau1[pf]; t2_tmp = tau2[pf]; nvec_tmp = norm_now[pf];
      # t1_tmp = e1_now[pf]; t2_tmp = e2_now[pf]; nvec_tmp = norm_now[pf];
      
      Amat = [t1_tmp[1] t1_tmp[2];t2_tmp[1] t2_tmp[2]] # matrix A  and vector d from (3.21)
    #   dvec = [nvec_tmp[1]; nvec_tmp[2]] # not needed for single correction
      
      # a,b,c calculated with θn,ϕn angles of normal direction in 3D, because I need to calculate distance from a point to that line
      ac_angle = tan(theta)*cos(phi)
      bc_angle = tan(theta)*sin(phi) # values needed for permutated setting

      # range of indices to check, so that all the nodes in the tubular neighborhood are considered
      Mx3p = min( Nvecp[3], Int(  ceil((Permzpt[3] + epsl - W3[1])/h)))
      Mx3m = max( 1, Int( floor((Permzpt[3] - epsl - W3[1])/h)))
      
      indIm = Mx3m:Mx3p # range of indices in z-direction (after permutation)
      indIc = 1:length(indIm) # how many indices in total, from 1 to ...

      # intersection of line with planes in permutated basis
      Int_zpt = [ [ac_angle*( W3[im]-Permzpt[3])+Permzpt[1]; bc_angle*(W3[im]-Permzpt[3])+Permzpt[2]; W3[im]] for im in indIm ]

      # indices corresponding to (bottom left) closest grid node to the intersection
      Indint_p = [ [Int(floor((Int_zpt[im][1]-W1[1])/h))+1; Int(floor((Int_zpt[im][2]-W2[1])/h))+1; indIm[im] ] for im in indIc ]

      # α,β for bottom left closest grid node
      tmp1 = [ [( Int_zpt[im][1] - ((Indint_p[im][1]-1)*h + W1[1]) )/h; ( Int_zpt[im][2] - ((Indint_p[im][2]-1)*h + W2[1]) )/h] for im in indIc ]
      
      # indices (only first two) corresponding to (actual) closest point to  to intersection
      Indint_p_w1 = [ [Indint_p[im][1]+(tmp1[im][1]>0.5); Indint_p[im][2]+(tmp1[im][2]>0.5)] for im in indIc ]
    
      # α,β for closest grid node
      ab_single_t[m, indIc] = [ [( Int_zpt[im][1] - ((Indint_p_w1[im][1]-1)*h + W1[1]) )/h; ( Int_zpt[im][2] - ((Indint_p_w1[im][2]-1)*h + W2[1]) )/h] for im in indIc ]

      # 3D indices of single correction nodes, correctly permutated to standard xyz system (pb, backward)
      Ind_tmp = [ ([ Indint_p_w1[im][1]; Indint_p_w1[im][2]; Indint_p[im][3] ])[pb] for im in indIc ]
      # third is from Indint_p, because the third is not saved into Indint_p_w1

      # 1D array of subindices of indIc for which the indices of the corr. nodes are not out of bounds, important for not passing bad values to indIJK_to_M
      ind_check = indIc[ [ prod( (ind_tmp.<=Nvec) .& (ind_tmp.>=1) ) for ind_tmp in Ind_tmp ].>0 ]

      # 1D array of m-indices of valid correction points; indices given to indIJK_to_M to get m-index
      (iv1_t[m])[ind_check] = [ indIJK_to_M[ind_tmp[1],ind_tmp[2],ind_tmp[3]] for ind_tmp in Ind_tmp[ind_check] ]

      for itmp in ind_check # indices sets which are not out of bounds
          η = (1-2*insidepoint((Int_zpt[itmp])[pb]))*norm(Permzpt-Int_zpt[itmp]) 
          (eta_vec[m])[itmp] = η
          α, β = ab_single_t[m, itmp]
          D0_now = D0(η)

          # weights for the different kernels
          w_K11_single_t[ m ][ itmp ] = w_k0_ptilde1([α;β]; lfun=(t->i0_PB_a(t,D0_now,Amat, Mmat, epsl_ratio)))*secondt_const
          w_K22_single_t[ m ][ itmp ] = w_k0_ptilde1([α;β]; lfun=(t->j0_PB_a(t,D0_now,Amat, Mmat, 1/epsl_ratio)))*(-secondt_const)
          w_K21_single_t[ m ][ itmp ] = w_k0_ptilde1([α;β]; lfun=(t->1.0./ψ0_PB_a(t,D0_now,Amat)))*(secondt_const)*0.5*kappa_val^2
          # w_K12_single_t[ m ][ itmp ] = secondt_const*kappa_val # weight is always one
          # technically useless array
      end
      
      # only keep the values corresponding to nodes inside the neighborhod
      p1tmp = (iv1_t[m]) .> 0
      iv1_t[m] = (iv1_t[m])[p1tmp]
      
      eta_vec[ m ] = (eta_vec[ m ])[p1tmp]     
      
      w_K11_single_t[ m ] = (w_K11_single_t[ m ])[p1tmp]
      w_K22_single_t[ m ] = (w_K22_single_t[ m ])[p1tmp]
      # w_K12_single_t[ m ] = (w_K12_single_t[ m ])[p1tmp]
      w_K21_single_t[ m ] = (w_K21_single_t[ m ])[p1tmp]

    end
  end
  println("Correction arrays filled - targets")
  # println("Only targets considered")

  println("--- Precomputations over ---")

  return M, Pg, PΓ, indIJK_to_M, indM_to_IJK, dsignes, iv1, w_K11_single, w_K22_single, w_K21_single, normal, iv1_t, w_K11_single_t, w_K22_single_t, w_K21_single_t, target_normal
end

#########################################################
#########################################################
# Function potentials, matrix-vector functions
#########################################################
#########################################################

tmprand = [ rand(3) for i=1:20];
tmprands = rand(20)
tmprandk = rand(20,2)

tmprandInc = [ ones(Int64,5) for i=1:20]
for j=1:20
    for i=1:4
        tmp = setdiff(1:20,tmprandInc[j][1:i])
        tmprandInc[j][i+1] = rand(tmp,1)[1]
    end
end

tmprandW = [ rand(5) for i=1:20 ]
tmprandD = [ rand(5) for i=1:20 ]
tmprandLim = [ rand(5) for i=1:20 ]

tmprandInc2 = [ tmprandInc[m][1:end-2] for m in 1:length(tmprandInc)]
tmprandW1 = [ rand(3) for i=1:20 ]

function K11_PB_Q1corr_target(α::Array{Float64,1}; source::Array{Array{Float64,1},1}=tmprand, 
  normal::Array{Array{Float64,1},1}=tmprand, salvare::Array{Float64,1}=tmprands,
  iv1::Array{Array{Int64,1},1}=tmprandInc2 ,
  wv1::Array{Array{Float64,1},1}=tmprandW1 , 
  h::Real=0.1,    
  targets::Array{Array{Float64,1},1}=tmprand,
  kappa_val::Real=1.0, theta_val::Real=1.0)

  Mt = length(targets)
  Mm = length(α)
  Aα = Array{Float64,1}(undef,Mt);

  @threads for m=1:Mt
      ind_1 = Array{Int64,1}(1:Mm);
      setdiff!( ind_1, iv1[m] )
      Aα[m] = sum( [ α[l] * salvare[l] * dGny_diff(targets[m], source[l], normal[l], kappa_val, theta_val) for l in ind_1 ] )*h^3
      Aα[m] += sum( α[iv1[m]] .* salvare[iv1[m]] .*  wv1[m] )*h^2 
  end
  return Aα
end
K11_PB_Q1corr_target(ones(20))
K11_PB_Q1corr_target(ones(20))
K11_PB_Q1corr_target(ones(20))

function K11_PB_Q1corr_DEBUG(α::Array{Float64,1}; source::Array{Array{Float64,1},1}=tmprand, 
  normal::Array{Array{Float64,1},1}=tmprand, salvare::Array{Float64,1}=tmprands,
  iv1::Array{Array{Int64,1},1}=tmprandInc2 ,
  wv1::Array{Array{Float64,1},1}=tmprandW1 , 
  h::Real=0.1,    
  targets::Array{Array{Float64,1},1}=tmprand,
  kappa_val::Real=1.0, theta_val::Real=1.0)

  Mt = length(targets)
  Mm = length(α)
  Aα = Array{Float64,1}(undef,Mt);

  for m=1:Mt
      ind_1 = Array{Int64,1}(1:Mm);
      setdiff!( ind_1, iv1[m] )
      for l in ind_1
        if isnan(dGny_diff(targets[m], source[l], normal[l], kappa_val, theta_val))
          println("dGdny  = $(dGny_diff(targets[m], source[l], normal[l], kappa_val, theta_val))")
          println("index  = target $(m), source $(l)")
          println("diff   = ",norm(targets[m].-source[l]))
          println("target = $(targets[m])")
          println("source = $(source[l])")
          println("normal = $(normal[l])")
        end
      end
      Aα[m] = sum( [ α[l] * salvare[l] * dGny_diff(targets[m], source[l], normal[l], kappa_val, theta_val) for l in ind_1 ] )*h^3
      Aα[m] += sum( α[iv1[m]] .* salvare[iv1[m]] .*  wv1[m] )*h^2 
  end
  return Aα
end
# K11_PB_Q1corr_DEBUG(ones(20))
# K11_PB_Q1corr_DEBUG(ones(20))
# K11_PB_Q1corr_DEBUG(ones(20))

function K22_PB_Q1corr_target(α::Array{Float64,1}; source::Array{Array{Float64,1},1}=tmprand, 
  salvare::Array{Float64,1}=tmprands,
  iv1::Array{Array{Int64,1},1}=tmprandInc2 ,
  wv1::Array{Array{Float64,1},1}=tmprandW1 , 
  h::Real=0.1,    
  targets::Array{Array{Float64,1},1}=tmprand,
  targetnormals::Array{Array{Float64,1},1}=tmprand,
  kappa_val::Real=1.0, theta_val::Real=1.0)

  Mt = length(targets)
  Mm = length(α)
  Aα = Array{Float64,1}(undef,Mt);

  @threads for m=1:Mt
      ind_1 = Array{Int64,1}(1:Mm);
      setdiff!( ind_1, iv1[m] )
      Aα[m] = sum( [ α[l] * salvare[l] * dGnx_diff(targets[m], source[l], targetnormals[m], kappa_val, theta_val) for l in ind_1 ] )*h^3
      Aα[m] += sum( α[iv1[m]] .* salvare[iv1[m]] .*  wv1[m] )*h^2 
  end
  return Aα
end
K22_PB_Q1corr_target(ones(20))
K22_PB_Q1corr_target(ones(20))
K22_PB_Q1corr_target(ones(20))

function K21_PB_Q1corr_target(α::Array{Float64,1}; source::Array{Array{Float64,1},1}=tmprand, 
  normal::Array{Array{Float64,1},1}=tmprand, salvare::Array{Float64,1}=tmprands,
  iv1::Array{Array{Int64,1},1}=tmprandInc2 ,
  wv1::Array{Array{Float64,1},1}=tmprandW1 , 
  h::Real=0.1,    
  targets::Array{Array{Float64,1},1}=tmprand,
  targetnormals::Array{Array{Float64,1},1}=tmprand,
  kappa_val::Real=1.0)

  Mt = length(targets)
  Mm = length(α)
  Aα = Array{Float64,1}(undef,Mt);

  @threads for m=1:Mt
      ind_1 = Array{Int64,1}(1:Mm);
      setdiff!( ind_1, iv1[m] )
      Aα[m] = sum( [ α[l] * salvare[l] * d2G_diff(targets[m], source[l], targetnormals[m], normal[l], kappa_val) for l in ind_1 ] )*h^3
      Aα[m] += sum( α[iv1[m]] .* salvare[iv1[m]] .*  wv1[m] )*h^2 
  end
  return Aα
end
K21_PB_Q1corr_target(ones(20))
K21_PB_Q1corr_target(ones(20))
K21_PB_Q1corr_target(ones(20))

function K12_PB_Q1corr_target(α::Array{Float64,1}; source::Array{Array{Float64,1},1}=tmprand, 
  salvare::Array{Float64,1}=tmprands,
  iv1::Array{Array{Int64,1},1}=tmprandInc2 ,
  #wv1::Array{Array{Float64,1},1}=tmprandW1 , 
  h::Real=0.1,    
  targets::Array{Array{Float64,1},1}=tmprand,
  kappa_val::Real=1.0)

  Mt = length(targets)
  Mm = length(α)
  Aα = Array{Float64,1}(undef,Mt);

  @threads for m=1:Mt
      ind_1 = Array{Int64,1}(1:Mm);
      setdiff!( ind_1, iv1[m] )
      Aα[m] = sum( [ α[l] * salvare[l] * G0Gk(targets[m], source[l], kappa_val) for l in ind_1 ] )*h^3
      Aα[m] += sum( α[iv1[m]] .* salvare[iv1[m]] )*(h^3)*(kappa_val/(4*pi))
  end
  return Aα
end
K12_PB_Q1corr_target(ones(20))
K12_PB_Q1corr_target(ones(20))
K12_PB_Q1corr_target(ones(20))

# 2D functions for testing on 2D planes
function K11_PB_Q1corr_target2d(α::Array{Float64,1}; source::Array{Array{Float64,1},1}=tmprand, 
  normal::Array{Array{Float64,1},1}=tmprand, salvare::Array{Float64,1}=tmprands,
  iv1::Array{Array{Int64,1},1}=tmprandInc2 ,
  wv1::Array{Array{Float64,1},1}=tmprandW1 , 
  h::Real=0.1,    
  targets::Array{Array{Float64,1},1}=tmprand,
  kappa_val::Real=1.0, theta_val::Real=1.0)

  Mt = length(targets)
  Mm = length(α)
  Aα = Array{Float64,1}(undef,Mt);

  for m=1:Mt
      ind_1 = Array{Int64,1}(1:Mm);
      setdiff!( ind_1, iv1[m] )
      Aα[m] = sum( [ α[l] * salvare[l] * dGny_diff(targets[m], source[l], normal[l], kappa_val, theta_val) for l in ind_1 ] )*h^2
      Aα[m] += sum( α[iv1[m]] .* salvare[iv1[m]] .*  wv1[m] )*h
  end
  return Aα
end
K11_PB_Q1corr_target2d(ones(20))
K11_PB_Q1corr_target2d(ones(20))
K11_PB_Q1corr_target2d(ones(20))

function K22_PB_Q1corr_target2d(α::Array{Float64,1}; source::Array{Array{Float64,1},1}=tmprand, 
  salvare::Array{Float64,1}=tmprands,
  iv1::Array{Array{Int64,1},1}=tmprandInc2 ,
  wv1::Array{Array{Float64,1},1}=tmprandW1 , 
  h::Real=0.1,    
  targets::Array{Array{Float64,1},1}=tmprand,
  targetnormals::Array{Array{Float64,1},1}=tmprand,
  kappa_val::Real=1.0, theta_val::Real=1.0)

  Mt = length(targets)
  Mm = length(α)
  Aα = Array{Float64,1}(undef,Mt);

  for m=1:Mt
      ind_1 = Array{Int64,1}(1:Mm);
      setdiff!( ind_1, iv1[m] )
      Aα[m] = sum( [ α[l] * salvare[l] * dGnx_diff(targets[m], source[l], targetnormals[m], kappa_val, theta_val) for l in ind_1 ] )*h^2
      Aα[m] += sum( α[iv1[m]] .* salvare[iv1[m]] .*  wv1[m] )*h
  end
  return Aα
end
K22_PB_Q1corr_target2d(ones(20))
K22_PB_Q1corr_target2d(ones(20))
K22_PB_Q1corr_target2d(ones(20))


################ IBIM potentials
function K11_PB_IBIM_target(α::Array{Float64,1}; source::Array{Array{Float64,1},1}=tmprand, normal::Array{Array{Float64,1},1}=tmprand, salvare::Array{Float64,1}=tmprands, targets::Array{Array{Float64,1},1}=tmprand, kappa_val::Real=1.0, theta_val::Real=1.0, tau::Real=0.1)
  Mm = length(targets)
  Aα = Array{Float64,1}(undef,Mm);
  @threads for i=1:Mm
    Aα[i] = sum( [ α[j] * salvare[j] * K11_PB(targets[i], source[j], normal[j], kappa_val, theta_val, tau) for j in 1:length(α) ] )
  end
  return Aα
end
K11_PB_IBIM_target(ones(20))
K11_PB_IBIM_target(ones(20))
K11_PB_IBIM_target(ones(20))

function K22_PB_IBIM_target(α::Array{Float64,1}; source::Array{Array{Float64,1},1}=tmprand, salvare::Array{Float64,1}=tmprands, targets::Array{Array{Float64,1},1}=tmprand, targetnormals::Array{Array{Float64,1},1}=tmprand, kappa_val::Real=1.0, theta_val::Real=1.0, tau::Real=0.1)
  Mm = length(targets)
  Aα = Array{Float64,1}(undef,Mm);
  @threads for i=1:Mm
    Aα[i] = sum( [ α[j] * salvare[j] * K22_PB(targets[i], source[j], targetnormals[i], kappa_val, theta_val, tau) for j in 1:length(α) ] )
  end
  return Aα
end
K22_PB_IBIM_target(ones(20))
K22_PB_IBIM_target(ones(20))
K22_PB_IBIM_target(ones(20))

function K21_PB_IBIM_target(α::Array{Float64,1}; source::Array{Array{Float64,1},1}=tmprand, normal::Array{Array{Float64,1},1}=tmprand, salvare::Array{Float64,1}=tmprands, targets::Array{Array{Float64,1},1}=tmprand, targetnormals::Array{Array{Float64,1},1}=tmprand, kappa_val::Real=1.0, tau::Real=0.1)
  Mm = length(targets)
  Aα = Array{Float64,1}(undef,Mm);
  @threads for i=1:Mm
    Aα[i] = sum( [ α[j] * salvare[j] * K21_PB(targets[i], source[j], targetnormals[i], normal[j], kappa_val, tau) for j in 1:length(α) ] )
  end
  return Aα
end
K21_PB_IBIM_target(ones(20))
K21_PB_IBIM_target(ones(20))
K21_PB_IBIM_target(ones(20))

function K12_PB_IBIM_target(α::Array{Float64,1}; source::Array{Array{Float64,1},1}=tmprand, salvare::Array{Float64,1}=tmprands, targets::Array{Array{Float64,1},1}=tmprand, kappa_val::Real=1.0, tau::Real=0.1)
  Mm = length(targets)
  Aα = Array{Float64,1}(undef,Mm);
  @threads for i=1:Mm
    Aα[i] = sum( [ α[j] * salvare[j] * K12_PB(targets[i], source[j], kappa_val, tau) for j in 1:length(α) ] )
  end
  return Aα
end
K12_PB_IBIM_target(ones(20))
K12_PB_IBIM_target(ones(20))
K12_PB_IBIM_target(ones(20))

#############################################3

#############################################3
# # functions for finite differences
#############################################3

#################################################################

vecval_f2_2nd_c = [1;-2;1]

function gradient_distance_2nd(Avec,Bvec,Cvec,h)
  return [dot(vecval2nd_c, Avec);dot(vecval2nd_c, Bvec);dot(vecval2nd_c, Cvec)]/h
end

global c1dervec = [1/12;-2/3;0;2/3;-1/12];
# firstder for 4th order 1st derivative (first method) and 4th order mixed 2nd derivatives (second method)
# the vectors need to be with 5 components, and the matrix needs to be 5x5

# first method
function firstder_2nd(v::Vector{T}) where {T<:Real}
  return dot([-1/2;0;1/2],v)
end
# second method
function firstder_2nd(V::Matrix{T}) where {T<:Real}
  return firstder_2nd([ firstder_2nd(V[:,i]) for i=1:size(V,2) ])
end

function hessian3D_2nd(A::Array{T,3}, dx, dy, dz) where {T<:Real}
  # I expect A of size 3x3x3 for 2nd order approximations
  Hmat = zeros(3,3);
  Hmat[1,1] = dot(vecval_f2_2nd_c,A[:,2,2])/(dx^2)
  Hmat[2,2] = dot(vecval_f2_2nd_c,A[2,:,2])/(dy^2)
  Hmat[3,3] = dot(vecval_f2_2nd_c,A[2,2,:])/(dz^2)

  Hmat[1,2] = firstder_2nd(A[:,:,2])/((dx)*(dy))
  Hmat[1,3] = firstder_2nd(A[:,2,:])/((dx)*(dz))
  Hmat[2,3] = firstder_2nd(A[2,:,:])/((dz)*(dy))
  Hmat[2,1] = Hmat[1,2]
  Hmat[3,1] = Hmat[1,3]
  Hmat[3,2] = Hmat[2,3]
  return Hmat
end

vecval_f2_2nd_f = [2;-5;4;-1]
vecval_f2_2nd_b = -vecval_f2_2nd_f[end:-1:1]
vecval2nd_c = [-1/2;0;1/2]
vecval1st_f = [-1;1]
vecval2nd_A = [-3 4 -1; 1 -4 3]/2
vecval3rd_A = [-11 18 -9 2; -2 9 -18 11]/6
vecval4th_A = [-25/12 4 -3 4/3 -1/4; 1/4 -4/3 3 -4 25/12]

vecval4th_c = [1/12;-2/3;0;2/3;-1/12]
vecval_f2_4th_c = [-1/12;4/3;-5/2;4/3;-1/12]

function get_indices_b(Vec::Array{Float64,1}) # returns 0:2 or -2:0
  Sp = dot(vecval2nd_c, Vec[3:5])
  Sm = dot(vecval2nd_c, Vec[1:3])
  vs = Sp-Sm
  if abs(vs)<1e-3
    return -1:1, vecval2nd_c
  elseif vs<0
    return 0:1:2, vecval2nd_A[1,:]
  else 
    return -2:1:0, vecval2nd_A[2,:]
  end
end

function numjac3D_AB_2nd!(i::Int64, j::Int64, k::Int64, Pg::Array{Float64}, A::Array{Float64,2}, B::Array{Float64,2}, C::Array{Float64,2})
  N=size(Pg,1);

  if (i<3)
    A[:,1]=Pg[i:i+4,j,k,1];
    A[:,2]=Pg[i:i+4,j,k,2];
    A[:,3]=Pg[i:i+4,j,k,3];
  elseif ((i+2)>N)
    A[1:5,1]=Pg[i-4:i,j,k,1];
    A[1:5,2]=Pg[i-4:i,j,k,2];
    A[1:5,3]=Pg[i-4:i,j,k,3];
  else
    A[1:5,1]=Pg[i-2:i+2,j,k,1];
    A[1:5,2]=Pg[i-2:i+2,j,k,2];
    A[1:5,3]=Pg[i-2:i+2,j,k,3];
  end

  if (j<3)
    B[1:5,1]=Pg[i,j:j+4,k,1];
    B[1:5,2]=Pg[i,j:j+4,k,2];
    B[1:5,3]=Pg[i,j:j+4,k,3];
  elseif ((j+2)>N)
    B[1:5,1]=Pg[i,j-4:j,k,1];
    B[1:5,2]=Pg[i,j-4:j,k,2];
    B[1:5,3]=Pg[i,j-4:j,k,3];
  else
    B[1:5,1]=Pg[i, j-2:j+2, k,1];
    B[1:5,2]=Pg[i, j-2:j+2, k,2];
    B[1:5,3]=Pg[i, j-2:j+2, k,3];
  end

  if (k<3)
    C[1:5,1]=Pg[i,j,k:k+4,1];
    C[1:5,2]=Pg[i,j,k:k+4,2];
    C[1:5,3]=Pg[i,j,k:k+4,3];
  elseif ((k+2)>N)
    C[1:5,1]=Pg[i,j,k-4:k,1];
    C[1:5,2]=Pg[i,j,k-4:k,2];
    C[1:5,3]=Pg[i,j,k-4:k,3];
  else
    C[1:5,1]=Pg[i, j, k-2:k+2, 1];
    C[1:5,2]=Pg[i, j, k-2:k+2, 2];
    C[1:5,3]=Pg[i, j, k-2:k+2, 3];
  end
end

function gennumjac3D_2nd(Nx,Ny,Nz,A,B,C,h,i,j,k)

  #A serve per lo stencil in i, B per lo stencil in j
  # A e B hanno dimensione 9x9
  dPg = zeros(3,3);
  cvec = [-25/12 4 -3 4/3 -0.25]/h;
  cvec2 = [0.25 -4/3 3 -4 25/12]/h;
  # cvec2 = -cvec[end:-1:1];
  cvec3 = [1.0/12 -2.0/3 0.0 2.0/3 -1.0/12]/h

  if (i<3)
    dPg[1,1] = (cvec*A[1:5,1])[1]
    dPg[1,2] = (cvec*B[1:5,1])[1]
    dPg[1,3] = (cvec*C[1:5,1])[1]
  elseif (i>(Nx-2))
    dPg[1,1] = dot(cvec2,A[1:5,1])
    dPg[1,2] = dot(cvec2,B[1:5,1])
    dPg[1,3] = dot(cvec2,C[1:5,1])
  else
    dPg[1,1] = dot(cvec3, A[1:5,1])
    dPg[1,2] = dot(cvec3, B[1:5,1])
    dPg[1,3] = dot(cvec3, C[1:5,1])
  end

  if (j<3)
    dPg[2,1] = (cvec*A[1:5,2])[1]
    dPg[2,2] = (cvec*B[1:5,2])[1]
    dPg[2,3] = (cvec*C[1:5,2])[1]
  elseif (j>(Ny-2))
    dPg[2,1] = dot(cvec2,A[1:5,2])
    dPg[2,2] = dot(cvec2,B[1:5,2])
    dPg[2,3] = dot(cvec2,C[1:5,2])
  else
    dPg[2,1] = dot(cvec3, A[1:5,2])
    dPg[2,2] = dot(cvec3, B[1:5,2])
    dPg[2,3] = dot(cvec3, C[1:5,2])
  end

  if (k<3)
    dPg[3,1] = (cvec*A[1:5,3])[1]
    dPg[3,2] = (cvec*B[1:5,3])[1]
    dPg[3,3] = (cvec*C[1:5,3])[1]
  elseif (k>(Nz-2))
    dPg[3,1] = dot(cvec2,A[1:5,3])
    dPg[3,2] = dot(cvec2,B[1:5,3])
    dPg[3,3] = dot(cvec2,C[1:5,3])
  else
    dPg[3,1] = dot(cvec3, A[1:5,3])
    dPg[3,2] = dot(cvec3, B[1:5,3])
    dPg[3,3] = dot(cvec3, C[1:5,3])
  end
  return dPg
end

#################################################################

vecvalm2 = [1/4;-4/3;3;-4;25/12]
vecvalm1 = [-1/12;1/2;-3/2;5/6;1/4]
vecval0 = [1/12;-2/3;0;2/3;-1/12]
vecvalp1 = -vecvalm1[end:-1:1]
vecvalp2 = -vecvalm2[end:-1:1]
global ValMatrix = zeros(5,5);
ValMatrix[1,:] = vecvalm2
ValMatrix[2,:] = vecvalm1
ValMatrix[3,:] = vecval0
ValMatrix[4,:] = vecvalp1
ValMatrix[5,:] = vecvalp2

vecval2nd_p = [-3;4;-1]/2
vecval2nd_m = -[-3;4;-1]/2

function gradient_distance_WENO(Avec,Bvec,Cvec,h)
  SAp = dot(vecval2nd_c,Avec[3:5])
  SAm = dot(vecval2nd_c,Avec[1:3])
  vSa = abs(SAp) <= abs(SAm)
  retA = (vSa)*dot(vecval2nd_p,Avec[3:5]) + (!vSa)*dot(vecval2nd_p,Avec[1:3])
  SBp = dot(vecval2nd_c,Bvec[3:5])
  SBm = dot(vecval2nd_c,Bvec[1:3])
  vSb = abs(SBp) <= abs(SBm)
  retB = (vSb)*dot(vecval2nd_p,Bvec[3:5]) + (!vSb)*dot(vecval2nd_p,Bvec[1:3])
  SCp = dot(vecval2nd_c,Cvec[3:5])
  SCm = dot(vecval2nd_c,Cvec[1:3])
  vSc = abs(SCp) <= abs(SCm)
  retC = (vSc)*dot(vecval2nd_p,Cvec[3:5]) + (!vSc)*dot(vecval2nd_p,Cvec[1:3])
  return [retA; retB; retC]/h
end

#############################################3

#############################################3
# # plotting functions
#############################################3

function plotting_comparison_dir(dir, nr, surface)
  hval = readdlm(dir*"/$(nr)/surf_PB_$(nr)_hvals.dat")[:,1]
  val_abs = readdlm(dir*"/$(nr)/surf_PB_$(nr)_val_abs.dat")[:,1]
  s = readdlm(dir*"/$(nr)/surf_PB_$(nr)_val_size.dat")[:,1]; s = Int.(s);
  surf_val = readdlm(dir*"/$(nr)/surf_PB_$(nr)_surf_val.dat")[:,1]
  sb = readdlm(dir*"/$(nr)/surf_PB_$(nr)_surf_val_size.dat")[:,1]; sb = Int.(sb);
  val_abs = reshape(val_abs, (s[1], s[2], s[3]))
  surf_val = reshape(surf_val, (sb[1], sb[2], sb[3]))
  detail = join( readdlm(dir*"/$(nr)/surf_PB_$(nr)_detail.dat") )
  plotting_comparison_1(hval, val_abs, surf_val, dir, nr, surface, detail)
end

function plotting_comparison_1(hval, val_abs, surf_val, dir, nr, surface, detail)
  nrun[1]+=1
  begin # ref values and all errors 
    refK11 = mean(val_abs[5,:,end])
    refK22 = mean(val_abs[6,:,end])
    refK21 = mean(val_abs[7,:,end])
    refK12 = mean(val_abs[8,:,end])

    refbK11 = mean(val_abs[1,:,end])
    refbK22 = mean(val_abs[2,:,end])
    refbK21 = mean(val_abs[3,:,end])
    refbK12 = mean(val_abs[4,:,end])

    err_K11_IBIM = abs.(mean(val_abs[1,:,1:end-1], dims=1).-refK11)/abs(refK11)
    err_K11_corr   = abs.(mean(val_abs[5,:,1:end-1], dims=1).-refK11)/abs(refK11)

    err_K22_IBIM = abs.(mean(val_abs[2,:,1:end-1], dims=1).-refK22)/abs(refK22)
    err_K22_corr   = abs.(mean(val_abs[6,:,1:end-1], dims=1).-refK22)/abs(refK22)
    
    err_K21_IBIM = abs.(mean(val_abs[3,:,1:end-1], dims=1).-refK21)/abs(refK21)
    err_K21_corr   = abs.(mean(val_abs[7,:,1:end-1], dims=1).-refK21)/abs(refK21)
    
    err_K12_IBIM = abs.(mean(val_abs[4,:,1:end-1], dims=1).-refK12)/abs(refK12)
    err_K12_corr   = abs.(mean(val_abs[8,:,1:end-1], dims=1).-refK12)/abs(refK12)

    errb_K11_IBIM = abs.(mean(val_abs[1,:,1:end-1], dims=1).-refbK11)/abs(refbK11)
    errb_K11_corr   = abs.(mean(val_abs[5,:,1:end-1], dims=1).-refbK11)/abs(refbK11)

    errb_K22_IBIM = abs.(mean(val_abs[2,:,1:end-1], dims=1).-refbK22)/abs(refbK22)
    errb_K22_corr   = abs.(mean(val_abs[6,:,1:end-1], dims=1).-refbK22)/abs(refbK22)
    
    errb_K21_IBIM = abs.(mean(val_abs[3,:,1:end-1], dims=1).-refbK21)/abs(refbK21)
    errb_K21_corr   = abs.(mean(val_abs[7,:,1:end-1], dims=1).-refbK21)/abs(refbK21)
    
    errb_K12_IBIM = abs.(mean(val_abs[4,:,1:end-1], dims=1).-refbK12)/abs(refbK12)
    errb_K12_corr   = abs.(mean(val_abs[8,:,1:end-1], dims=1).-refbK12)/abs(refbK12)

    diff_K11 = abs.(diff(mean(val_abs[5,:,:], dims=1), dims=2))
    diff_K22 = abs.(diff(mean(val_abs[6,:,:], dims=1), dims=2))
    diff_K21 = abs.(diff(mean(val_abs[7,:,:], dims=1), dims=2))
    diff_K12 = abs.(diff(mean(val_abs[8,:,:], dims=1), dims=2))

    s1 = diff(surf_val[1,:])
    s2 = diff(surf_val[2,:])
  end
  
  # ref corr
  figure(150, figsize=(8,6)); clf()
  subplot(1,2,1)
  loglog(hval[1:end-1], err_K11_IBIM[1,:], "-+b",label="K11")
  loglog(hval[1:end-1], err_K22_IBIM[1,:], "-^r",label="K22")
  loglog(hval[1:end-1], err_K21_IBIM[1,:], "-xm",label="K21")
  loglog(hval[1:end-1], err_K12_IBIM[1,:], "-vk",label="K12")
  title("Ref=CTR. Method: IBIM")
  subplot(1,2,2)
  loglog(hval[1:end-1], err_K11_corr[1,:], "-+b")
  loglog(hval[1:end-1], err_K22_corr[1,:], "-^r")
  loglog(hval[1:end-1], err_K21_corr[1,:], "-xm")
  loglog(hval[1:end-1], err_K12_corr[1,:], "-vk")
  title("Method: CTR")
  # subplot(2,2,3)
  # loglog(hval[1:end-1], errb_K11_IBIM[1,:], "-+b")
  # loglog(hval[1:end-1], errb_K22_IBIM[1,:], "-^r")
  # loglog(hval[1:end-1], errb_K21_IBIM[1,:], "-xm")
  # loglog(hval[1:end-1], errb_K12_IBIM[1,:], "-vk")
  # title("Ref=IBIM. Method: IBIM")
  # subplot(2,2,4)
  # loglog(hval[1:end-1], errb_K11_corr[1,:], "-+b")
  # loglog(hval[1:end-1], errb_K22_corr[1,:], "-^r")
  # loglog(hval[1:end-1], errb_K21_corr[1,:], "-xm")
  # loglog(hval[1:end-1], errb_K12_corr[1,:], "-vk")
  # title("Method: CTR")

  for i=1:2
      subplot(1,2,i)
      loglog(hval[1:end-1], 50*hval[1:end-1].^2, "--", label="O(h^2)")
      loglog(hval[1:end-1], 4*hval[1:end-1].^1.5, "--", label="O(h^{1.5})")
      loglog(hval[1:end-1], 2*hval[1:end-1].^2.5, "--", label="O(h^{2.5})")
      ylim([5e-7;1e-1])
  end
  subplot(1,2,1)
  legend(loc="best")
  tight_layout()
  savefig(dir*"/$(nr)/t$(nr)_$(nrun[1])_"*surface*"_"*detail*".png")

  figure(160, figsize=(8,6)); clf()
  subplot(1,2,1)
  loglog(hval[1:end-1], diff_K11[1,:], "-+b",label="K11")
  loglog(hval[1:end-1], diff_K22[1,:], "-^r",label="K22")
  loglog(hval[1:end-1], diff_K21[1,:], "-xm",label="K21")
  loglog(hval[1:end-1], diff_K12[1,:], "-vk",label="K12")
  loglog(hval, 5*hval.^1.5, "--", label="O(h^{1.5})")
  loglog(hval, 5*hval.^2, "--", label="O(h^2)")
  loglog(hval, 60*hval.^3, "--", label="O(h^3)")
  title("Kij. Diff. of consecutive values"); grid()
  legend(loc="best")  

  subplot(1,2,2)
  loglog(hval[1:end-1], s1, ".-k", label="ex. J")
  loglog(hval[1:end-1], s2, "-xr", label="2nd POS")
  loglog(hval[1:end-1], 5*hval[1:end-1].^2, "--", label="O(h^2)")
  loglog(hval[1:end-1], 60*hval[1:end-1].^3, "--", label="O(h^3)")
  loglog(hval[1:end-1], 2000*hval[1:end-1].^4, "--", label="O(h^4)")
  ylim([1e-7;3e-2])
  title("Surface area. Diff. of consecutive values")
  legend(loc="best"); grid()
  tight_layout()
  savefig(dir*"/$(nr)/t$(nr)_$(nrun[1])_surferr_"*surface*"_"*detail*".png")
end