
#########################
### rotations along x,y,z axes
Rx(a) = [1 0 0;0 cos(a) sin(a);0 -sin(a) cos(a)]
Ry(a) = [cos(a) 0 sin(a);0 1 0;-sin(a) 0 cos(a)]
Rz(a) = [cos(a) -sin(a) 0;sin(a) cos(a) 0;0 0 1]

function Pgamma_circle(x,R,x0)
    return (norm(x.-x0)>1e-3)*((x.-x0)*R/norm(x.-x0).+x0) + (norm(x.-x0)<=1e-3)*(x0 .+ R*x0/norm(x0))
end

function Pgamma_torus(z,Rv2,Rinv,R1,R2,xshift)
    zn = Rinv*(z.-xshift);
    phi = atan(zn[2],zn[1]);
    ctmp = R1*[cos(phi);sin(phi);0];
    return Rv2*((zn.-ctmp)/norm(zn.-ctmp)*R2.+ctmp).+xshift
end

function Pg_spheres(z, x0vec, Rvec, R3)
    delta = norm(x0vec[1].-x0vec[2])-Rvec[1]-Rvec[2];
    L1 = Rvec[1]+R3
    L2 = Rvec[2]+R3
    L3 = Rvec[1]+Rvec[2]+delta
    alpha = 0.5+(L1^2-L2^2)/(2L3^2)
    x1 = x0vec[1]
    x2 = x0vec[2]
    x3 = x0vec[1].+(alpha)*(x0vec[2].-x0vec[1])
    a = sqrt(abs(L1^2-alpha^2*L3^2))
    b = R3;
    X = x0vec[2].-x0vec[1]
    B = X/norm(X); 
    A = [0;1;0.]
    a1 = dot(A,B); a2 = norm(cross(A,B))
    G = [a1 -a2 0;a2 a1 0;0 0 1];
    u1 = A; u2 = B-a1*A; u2 /=norm(u2); u3 = cross(B,A);
    F = ([u1 u2 u3])
    Rottot = F*G*inv(F)
    Rot3 = transpose(Rottot)

    x4 = Rottot*[0;0;a].+x3
    cone_compA = (1-2*(norm(x2-x3)>norm(x1-x2)))
    sintC1 = a/norm(x4-x1); costC1 = norm(x1-x3)/norm(x4-x1); tant = sintC1/costC1
    cone_compB = (1-2*(norm(x1-x3)>norm(x1-x2)))
    sintC2 = a/norm(x4-x2); costC2 = norm(x2-x3)/norm(x4-x2); tant = sintC2/costC2
    
    p2 = Rx(-pi*0.5)*Rot3*(z.-x1)
    p3 = Rx(pi*0.5)*Rot3*(z.-x2)
    p2v = isInCone(sintC1, costC1, norm(x1-x3),p2, cone_compA)
    p3v = isInCone(sintC2, costC2, norm(x2-x3),p3, cone_compB)

    if ( (cone_compA*cone_compB>0) && (p2v || p3v) ) || ( (cone_compA<0) && (p3v && ~p2v) ) || ( (cone_compB<0) && (p2v && ~p3v) )
        RotX = Rottot*Rx(-pi*0.5)
        RotXinv = transpose(RotX)
        xshift = x3;
        return Pgamma_torus(z,RotX,RotXinv,a,b,xshift), 1
    else
        Pz1 = Pgamma_circle(z, Rvec[1], x0vec[1])
        Pz2 = Pgamma_circle(z, Rvec[2], x0vec[2])
        if norm(z.-x1)<Rvec[1]
            return Pz1, -1
        elseif norm(z.-x2)<Rvec[2]
            return Pz2, -1
        else
            return (norm(z-Pz1)<norm(z-Pz2))*Pz1 .+ (norm(z-Pz1)>=norm(z-Pz2))*Pz2, -1
        end
    end
end

function spheres_torus_test(x0vec, Rvec, R3)
  
    delta = norm(x0vec[1].-x0vec[2])-Rvec[1]-Rvec[2];
    if delta<0
      @warn("Spheres intersecting")
    end
    L1 = Rvec[1]+R3
    L2 = Rvec[2]+R3
    L3 = Rvec[1]+Rvec[2]+delta
    alpha = 0.5+(L1^2-L2^2)/(2L3^2)
    x3 = x0vec[1].+(alpha)*(x0vec[2].-x0vec[1])
  
    x1 = x0vec[1]
    x2 = x0vec[2]
    figure(199); clf()
    tvec = 0:0.02:1;
    nt = length(tvec);
    Sphere = [ [sin(tvec[i]*pi)*cos(2*pi*tvec[j]);sin(tvec[i]*pi)*sin(2*pi*tvec[j]);cos(tvec[i]*pi)] for i=1:nt for j=1:nt]
    xs = [x[1] for x in Sphere]
    ys = [x[2] for x in Sphere]
    zs = [x[3] for x in Sphere]
    plot3D(x1[1],x1[2],x1[3],"ok")
    scatter3D(x1[1].+Rvec[1]*xs, x1[2].+Rvec[1]*ys, x1[3].+Rvec[1]*zs, s=0.6, c="k")
    plot3D(x2[1],x2[2],x2[3],"or")
    scatter3D(x2[1].+Rvec[2]*xs, x2[2].+Rvec[2]*ys, x2[3].+Rvec[2]*zs, s=0.6, c="k")
    plot3D(x3[1],x3[2],x3[3],"ob")
    
  
    a = sqrt(abs(L1^2-alpha^2*L3^2))
    b = R3;
    Torus = [ [(a+b*cos(2*tvec[i]*pi))*sin(2*tvec[j]*pi);b*sin(2*pi*tvec[i]);(a+b*cos(2*tvec[i]*pi))*cos(2*tvec[j]*pi)] for i=1:nt for j=1:nt]
    # Torus = [ [b*sin(2*pi*tvec[i]);(a+b*cos(2*tvec[i]*pi))*cos(2*tvec[j]*pi);(a+b*cos(2*tvec[i]*pi))*sin(2*tvec[j]*pi)] for i=1:nt for j=1:nt]
    # Torus = [ [b*sin(2*pi*tvec[i]);(a+b*cos(2*tvec[i]*pi))*cos(2*tvec[j]*pi);(a+b*cos(2*tvec[i]*pi))*sin(2*pi*tvec[j])] for i=1:nt for j=1:nt]
  
    # figure(200);clf()
    # tmp1 = x1.+xshift
    # plot3D(tmp1[1],tmp1[2],tmp1[3],"ok")
    # plot3D([0;1],[0;0],[0;0],"-k")
    # plot3D([0;0],[0;1],[0;0],"-k")
    # plot3D([0;0],[0;0],[0;1],"-k")
    # # scatter3D(tmp1[1].+Rvec[1]*xt, tmp1[2].+Rvec[1]*yt, tmp1[3].+Rvec[1]*zt, s=1, c="k")
    # tmp2 = x2.+xshift
    # plot3D(tmp2[1],tmp2[2],tmp2[3],"or")
    # # scatter3D(tmp2[1].+Rvec[2]*xt, tmp2[2].+Rvec[2]*yt, tmp2[3].+Rvec[2]*zt, s=1, c="k")
  
    # X = tmp2
    # println(tmp2)
    # Y = cross(X,rand(3))
    # X = X/norm(X); Y = Y/norm(Y);
    # psi = asin(X[2]/sqrt(1-X[3]^2))
    # phi = asin(Y[3]/sqrt(1-X[3]^2))
    # theta = asin(-X[3])
    # Rot3 = Rx(-phi)*Ry(-theta)*Rz(-psi)
    X = x0vec[2].-x0vec[1]
    B = X/norm(X); 
    A = [0;1;0.]
    a1 = dot(A,B); a2 = norm(cross(A,B))
    G = [a1 -a2 0;a2 a1 0;0 0 1];
    u1 = A; u2 = B-a1*A; u2 /=norm(u2); u3 = cross(B,A);
    F = ([u1 u2 u3])
    Rottot = F*G*inv(F)
    Rot3 = transpose(Rottot)


    tmp3 = Rot3*(x2.-x1)
    # plot3D(tmp3[1],tmp3[2],tmp3[3],"sr")
  
    Torus2 = [ Rottot*(x).+x3 for x in Torus]
    x4 = Rottot*[a*sin(2*tvec[1]*pi);0;a*cos(2*tvec[1]*pi)].+x3
    sintC1 = a/norm(x4-x1); costC1 = norm(x1-x3)/norm(x4-x1); tant = sintC1/costC1
    figure(199)
    xt = [x[1] for x in Torus2]
    yt = [x[2] for x in Torus2]
    zt = [x[3] for x in Torus2]
    scatter3D(xt,yt,zt, s=1, c="b")

    hcone = norm(x1-x3)
    cone_compA = (1-2*(norm(x2-x3)>norm(x1-x2))) # 1 if cones have to be intersected, -1 if they need to be subtracted
    Cone = cone_compA*[ [tvec[i]*hcone*tant*cos(tvec[j]*2*pi);tvec[i]*hcone*tant*sin(tvec[j]*2*pi);tvec[i]*hcone] for i=1:nt for j=1:nt ]
    ConeA = [ Rottot*Rx(pi*0.5)*x.+x1 for x in Cone ]
    xcA = [x[1] for x in ConeA]
    ycA = [x[2] for x in ConeA]
    zcA = [x[3] for x in ConeA]
    scatter3D(xcA,ycA,zcA, s=1, c="g")
    C3 = [ Rottot*[a*sin(2*tvec[j]*pi);0;a*cos(2*tvec[j]*pi)].+x3 for j=1:nt ]
    xc3 = [x[1] for x in C3]
    yc3 = [x[2] for x in C3]
    zc3 = [x[3] for x in C3]
    scatter3D(xc3,yc3,zc3, s=1, c="r")

    hcone = norm(x2-x3)
    sintC2 = a/norm(x4-x2); costC2 = norm(x2-x3)/norm(x4-x2); tant = sintC2/costC2
    cone_compB = (1-2*(norm(x1-x3)>norm(x1-x2))) # 1 if cones have to be intersected, -1 if they need to be subtracted
    Cone = cone_compB * [ [tvec[i]*hcone*tant*cos(tvec[j]*2*pi);tvec[i]*hcone*tant*sin(tvec[j]*2*pi);tvec[i]*hcone] for i=1:nt for j=1:nt ]
    ConeB = [ Rottot*Rx(-pi*0.5)*x.+x2 for x in Cone ]
    xcB = [x[1] for x in ConeB]
    ycB = [x[2] for x in ConeB]
    zcB = [x[3] for x in ConeB]
    scatter3D(xcB,ycB,zcB, s=1, c="g")
    C3 = [ Rottot*[a*sin(2*tvec[j]*pi);0;a*cos(2*tvec[j]*pi)].+x3 for j=1:nt ]
    xc3 = [x[1] for x in C3]
    yc3 = [x[2] for x in C3]
    zc3 = [x[3] for x in C3]
    scatter3D(xc3,yc3,zc3, s=1, c="r")

    xlim([-1;1].+x3[1])
    ylim([-1;1].+x3[2])
    zlim([-1;1].+x3[3])

    figure(200); clf()
    C1 = [ Rot3*(x.-x1) for x in ConeA ]
    C2 = [ Rot3*(x.-x2).+[0;norm(x1-x2);0] for x in ConeB ]
    xC1 = [x[1] for x in C1]
    yC1 = [x[2] for x in C1]
    zC1 = [x[3] for x in C1]
    # scatter3D(xC1,yC1,zC1, s=1, c="b")
    scatter3D(xcA,ycA,zcA, s=1, c="b")
    xC2 = [x[1] for x in C2]
    yC2 = [x[2] for x in C2]
    zC2 = [x[3] for x in C2]
    # scatter3D(xC2,yC2,zC2, s=1, c="r")

    nab = 20000;
    # pv = [ rand(3).*[1;0.5;1].-[0.5;0;0.5] for i=1:nab ]
    pv = [ rand(3) for i=1:nab ]
    p2v = zeros(nab)
    p3v = zeros(nab)
    for i=1:nab
        p1b = pv[i].*[0.4;2;1.5].+[0.3;0.;0.]
        p1a = pv[i].*[0.3;1.6;1.3].+[0.22;0.;0.]
        p2 = Rx(-pi*0.5)*Rot3*(p1a.-x1)
        p3 = Rx(pi*0.5)*Rot3*(p1b.-x2)
        p2v[i] = (1-2*isInCone(sintC1, costC1, norm(x1-x3),p2, cone_compA))
        p3v[i] = (1-2*isInCone(sintC2, costC2, norm(x2-x3),p3, cone_compB))
    end
    indins2 = p2v .<= 0
    indout2 = p2v .> 0
    indins3 = p3v .<= 0
    indout3 = p3v .> 0
    xp = [ (x.*[0.3;1.6;1.3].+[0.22;0.;0.])[1] for x in pv ]
    yp = [ (x.*[0.3;1.6;1.3].+[0.22;0.;0.])[2] for x in pv ]
    zp = [ (x.*[0.3;1.6;1.3].+[0.22;0.;0.])[3] for x in pv ]
    scatter3D(xp[indins2], yp[indins2], zp[indins2], s=1.5, c="r")
    cb = scatter3D(xp[indout2], yp[indout2], zp[indout2], s=1.5, c="k")
    colorbar(cb)
    xlabel("x")
    ylabel("y")
    zlabel("z")

    figure(201); clf()
    C1 = [ Rot3*(x.-x1) for x in ConeA ]
    C2 = [ Rot3*(x.-x2).+[0;norm(x1-x2);0] for x in ConeB ]
    # xC1 = [x[1] for x in C1]
    # yC1 = [x[2] for x in C1]
    # zC1 = [x[3] for x in C1]
    # scatter3D(xC1,yC1,zC1, s=1, c="b")
    xC2 = [x[1] for x in C2]
    yC2 = [x[2] for x in C2]
    zC2 = [x[3] for x in C2]
    xp2 = [ (x.*[0.4;2;1.5].+[0.3;0.;0.])[1] for x in pv ]
    yp2 = [ (x.*[0.4;2;1.5].+[0.3;0.;0.])[2] for x in pv ]
    zp2 = [ (x.*[0.4;2;1.5].+[0.3;0.;0.])[3] for x in pv ]
    scatter3D(xcB,ycB,zcB, s=1, c="b")
    scatter3D(xp2[indins3], yp2[indins3], zp2[indins3], s=1.5, c="r")
    scatter3D(xp2[indout3], yp2[indout3], zp2[indout3], s=1.5, c="k")
    xlabel("x")
    ylabel("y")
    zlabel("z")

    pv = [ rand(3).*[0.7;2.5;2.5].+[0.;-0.5;-0.5] for i=1:nab ]
    
    
    figure(202); clf()
    PzVecA = [ Pg_spheres(pv[i], x0vec, Rvec, R3) for i=1:nab ]
    # return PzVecA 
    
    PzVec = [ (PzVecA[i])[1] for i=1:nab ]
    PzVecS = [ (PzVecA[i])[2] for i=1:nab ]
    xpz = [ x[1] for x in PzVec ]
    ypz = [ x[2] for x in PzVec ]
    zpz = [ x[3] for x in PzVec ]
    scatter3D(xpz,ypz,zpz, s=1, c=PzVecS )
    scatter3D(xt,yt,zt, s=0.5, c="b")

    scatter3D(x1[1].+Rvec[1]*xs, x1[2].+Rvec[1]*ys, x1[3].+Rvec[1]*zs, s=0.6, c="k")
    # plot3D(x2[1],x2[2],x2[3],"or")
    scatter3D(x2[1].+Rvec[2]*xs, x2[2].+Rvec[2]*ys, x2[3].+Rvec[2]*zs, s=0.6, c="k")
    # scatter3D(x1[1].+Rvec[1]*xs, x1[2].+Rvec[1]*ys, x1[3].+Rvec[1]*zs, s=0.6, c="k")
    # plot3D(x2[1],x2[2],x2[3],"or")
    # scatter3D(x2[1].+Rvec[2]*xs, x2[2].+Rvec[2]*ys, x2[3].+Rvec[2]*zs, s=0.6, c="k")
    xlabel("x")
    
    xlim([-1;1].+x3[1])
    ylim([-1;1].+x3[2])
    zlim([-1;1].+x3[3])
    ylabel("y")
    zlabel("z")
    
    figure(203); clf()
    xpz2 = [ x[1] for x in pv ]
    ypz2 = [ x[2] for x in pv ]
    zpz2 = [ x[3] for x in pv ]
    indins = PzVecS.>0
    scatter3D(xpz2[indins],ypz2[indins],zpz2[indins], s=1, c="r" )
    
    # scatter3D(x1[1].+Rvec[1]*xs, x1[2].+Rvec[1]*ys, x1[3].+Rvec[1]*zs, s=0.6, c="k")
    # plot3D(x2[1],x2[2],x2[3],"or")
    # scatter3D(x2[1].+Rvec[2]*xs, x2[2].+Rvec[2]*ys, x2[3].+Rvec[2]*zs, s=0.6, c="k")
    xlabel("x")
    ylabel("y")
    zlabel("z")
end

function isInCone(sint, cost, h, p, compl)
    return ( ((dot(p[1:2],p[1:2])*cost^2 - p[3]^2*sint^2) <=0) && (compl*p[3]>=0) && (compl*p[3]<=h) )
end

function spheres_torus_test_2d(x0vec, Rvec, R3)

    delta = norm(x0vec[1].-x0vec[2])-Rvec[1]-Rvec[2];
    if delta<0
        @warn("Circles intersecting")
    end
    L1 = Rvec[1]+R3
    L2 = Rvec[2]+R3
    L3 = Rvec[1]+Rvec[2]+delta
    alpha = 0.5+(L1^2-L2^2)/(2L3^2)
    h = sqrt(abs(L1^2-alpha^2*L3^2))

    x1 = x0vec[1]
    x2 = x0vec[2]
    figure(299); clf()
    tvec = 0:0.02:1;
    nt = length(tvec);
    Circle = [ [cos(tvec[i]*2*pi);sin(tvec[i]*2*pi)] for i=1:nt ]
    xc = [x[1] for x in Circle]
    yc = [x[2] for x in Circle]
    plot(x1[1],x1[2],"ok")
    plot(x1[1].+Rvec[1]*xc,x1[2].+Rvec[1]*yc,"-k")
    plot(x2[1],x2[2],"ob")
    plot(x2[1].+Rvec[2]*xc,x2[2].+Rvec[2]*yc,"-k")

    x3 = x0vec[1].+(alpha)*(x0vec[2].-x0vec[1])
    plot(x3[1],x3[2],"or")
    x4 = x3 .+ h*[0;1]
    plot(x4[1],x4[2],"sr")
    plot(x4[1].+R3*xc,x4[2].+R3*yc,"-k")

    xshift = -x0vec[1];
    
end

function euler_test()
    # V = [0.86;0.8;0.4]; 
    # V = [0.06737652098535762, 0.42183435657796853, 0.3423207711891334]
    V = [1.06737652098535762, 0.42183435657796853, 0.3423207711891334]
    # X = V/norm(V)
    # tmp1 = rand(3); #tmp1 = [0;0;1.]
    # Y = cross(X,tmp1); Y /= norm(Y)
    # Z = cross(X,Y)
    # # psi = asin(X[2]/sqrt(1-X[3]^2))
    # # phi = asin(Y[3]/sqrt(1-X[3]^2))
    # # theta = asin(-X[3])
    # # Rot1 = Rx(phi)
    # # Rot2 = Ry(theta)
    # # Rot3 = Rz(psi)
    # alpha = acos(-Z[2]/sqrt(1-Z[3]^2))
    # beta = acos(Z[3])
    # gamma = acos(Y[3]/sqrt(1-Z[3]^2))
    # # Rot1 = Rx(alpha)
    # # Rot2 = Rz(beta)
    # # Rot3 = Rx(gamma)
    # Rot1 = Rz(alpha)
    # Rot2 = Rx(beta)
    # Rot3 = Rz(gamma)
    # Rottot = Rot1*Rot2*Rot3 #Rot3*Rot2*
    # Rotinv = transpose(Rottot)
    B = V/norm(V); 
    A = [0;1;0.]
    a1 = dot(A,B); a2 = norm(cross(A,B))
    G = [a1 -a2 0;a2 a1 0;0 0 1];
    u1 = A; u2 = B-a1*A; u2 /=norm(u2); u3 = cross(B,A);
    F = ([u1 u2 u3])
    # println(u1)
    # println(F)
    Rottot = F*G*inv(F)
    Rotinv = transpose(Rottot)

    figure(4); clf()
    plot3D(0,0,0,"ok")
    plot3D([0;1],[0;0],[0;0],"-k")
    plot3D([0;0],[0;1],[0;0],"-k")
    plot3D([0;0],[0;0],[0;1],"-k")
    plot3D(V[1],V[2],V[3],"xr")
    xlim([-1.5;1.5])
    ylim([-1.5;1.5])
    zlim([-1.5;1.5])
    
    p = [1;0;0]#*norm(X)
    plot3D(p[1],p[2],p[3],"ob")
    p2 = Rottot*p
    p2 *= norm(V)
    plot3D(p2[1],p2[2],p2[3],"+b")

    V2 = Rotinv*V
    # V2 = transpose(Rx(-alpha))*transpose(Rz(-beta))*transpose(Rx(-gamma))*V
    plot3D(V2[1],V2[2],V2[3],"sr")    
end

function torus_proj(z)
    x1 = [0.7588462711677448; 0.05775529878652241; 0.33318634166802674]
    x2 = [0.9604612655979043; 0.31100867895918594; 0.3439301850417964]
    xvec = [x1,x2]
    Rvec = [0.5;0.4]
    R3 = 0.2;
   
    delta = norm(x1-x2)-Rvec[1]-Rvec[2]
    L1 = Rvec[1]+R3
    L2 = Rvec[2]+R3
    L3 = Rvec[1]+Rvec[2]+delta
    alpha = 0.5+(L1^2-L2^2)/(2L3^2)
    x3 = xvec[1].+(alpha)*(xvec[2].-xvec[1])

    a = sqrt(L1^2-alpha^2*L3^2); b = R3
    X = xvec[2].-xvec[1]
    B = X/norm(X); 
    A = [0;1;0.]
    a1 = dot(A,B); a2 = norm(cross(A,B))
    G = [a1 -a2 0;a2 a1 0;0 0 1];
    u1 = A; u2 = B-a1*A; u2 /=norm(u2); u3 = cross(B,A);
    F = ([u1 u2 u3])
    Rot3inv = F*G*inv(F)
    Rot3 = transpose(Rot3inv)

    v = xvec[1]; w = xvec[2]; 
    vb = Rot3inv*v; wb = Rot3inv*w; zb = Rot3inv*z;
    cvec = [wb[1]-vb[1] wb[3]-vb[3];zb[1]-vb[1] zb[3]-vb[3]]\[-wb[2]+vb[2]; -zb[2]+vb[2]]
    local f(x) = cvec[1]*a*sin(x)+cvec[2]*a*cos(x)-vb[2]-cvec[1]*vb[1]-cvec[2]*vb[3]
    local fp(x) = cvec[1]*a*cos(x)-cvec[2]*a*sin(x)
    phi1 = mynewton(2pi*0.15,f,fp)
    phi2 = mynewton(2pi*0.85,f,fp)
    # figure(8); clf()
    # plot([0;2*pi],[0;0],"-k")
    # tvec = 0:0.02:1;
    # plot(tvec*2*pi,f.(tvec*2*pi),"-r")
    # plot(phi1,0,"ob")
    # plot(phi2,0,"or")
    # ub = 
end

function mynewton(x0,f,fp;TOL::Real=1e-9)
    x = x0;
    d = 1;
    maxIt = 10;
    nit = 0;
    while abs(d)>TOL && nit<maxIt
        d = f(x)/fp(x);
        x -= d;
        nit += 1;
    end
    return x
end

function cone(theta,h)
    figure(2); clf()
    tvec = 0:0.1:1;
    nt = length(tvec)
    Cone = [ [tvec[i]*h*tan(theta)*cos(tvec[j]*2*pi);tvec[i]*h*tan(theta)*sin(tvec[j]*2*pi);tvec[i]*h] for i=1:nt for j=1:nt ]
    xc = [x[1] for x in Cone]
    yc = [x[2] for x in Cone]
    zc = [x[3] for x in Cone]

    plot3D(0,0,0,"ok")
    plot3D([0;1],[0;0],[0;0],"-k")
    plot3D([0;0],[0;1],[0;0],"-k")
    plot3D([0;0],[0;0],[0;1],"-k")
    
    scatter3D(xc, yc, zc, s=1, c="b")
end

function f_para(a, x, y)
    n = length(x)
    # z = SharedArray{Float64,1}(n)
    z = zeros(n)
    @threads for i=1:n
        z[i] = sum( [ a[i,j]*x[j] for j=1:n ] )+y[i]
    end
    return z
end

function f_single(a, x, y)
    n = length(x)
    z = zeros(n)
    for i=1:n
        z[i] = sum( [ a[i,j]*x[j] for j=1:n ] )+y[i]
    end
    return z
end