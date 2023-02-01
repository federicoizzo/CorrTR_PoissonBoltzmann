Rx(a) = [1 0 0;0 cos(a) sin(a);0 -sin(a) cos(a)]
Ry(a) = [cos(a) 0 sin(a);0 1 0;-sin(a) 0 cos(a)]
Rz(a) = [cos(a) -sin(a) 0;sin(a) cos(a) 0;0 0 1]
xshift = [0.05475547095598521; 0.06864792402110276; 0.03502726366462485]
# xshift = zeros(3)

a,b,c = 0.24402412255508432, 0.7454097947651017, 2.219760487439292
# a,b,c = 1.99487, 2.54097947651017, 4.219760487439292

# a,b,c = 0., 0., 0.
Rot = Rz(c)*Ry(b)*Rx(a);
Rinv = Rx(-a)*Ry(-b)*Rz(-c);

function toroid1(N1::Int64, N2::Int64; ifplot::Bool=false)
    
    # xvec2 = 0:0.01:1;
    xvec1 = LinRange(0,1, N1+1); xvec1 = xvec1[1:end-1]
    xvec2 = LinRange(0,1, N2+1); xvec2 = xvec2[1:end-1]
    # xvec2 = 0:0.0005:1;
    # xvec1 = [0.2]

    tvec = xvec1*2*pi;
    pvec = xvec2*2*pi;

    Nt = length(tvec); Np = length(pvec);

    Arot1(t) = [cos(2*t) sin(2*t);-sin(2*t) cos(2*t)];
    Arot2(t) = [cos(2*t) -sin(2*t);sin(2*t) cos(2*t)];
    # Pvec(t,a,b) = Arot2(t)*[a 0;0 b]*Arot1(t)*[cos(t); sin(t)];
    # Pvec(t,a,b) = Arot2(t)*[a 0;0 b]*[cos(t); sin(t)];

    RZfun(t,p,a,b,R) = [R;0] .+ Arot2(t)*[a 0;0 b]*Arot1(t)*[cos(p);sin(p)];
    # Rot = 
    

    # X = zeros(Nt*Np,3);
    X = Array{Array{Float64,1},1}(undef,Nt*Np);
    @threads for i=1:Nt
        for j=1:Np
            R =0.8; a = 0.4; b = 0.3;
            tmp = RZfun(tvec[i],pvec[j],a,b,R);
            # X[j+(i-1)*Np,:] = [(R + (Pvec(tvec[i],a,b))[1])*cos(pvec[j]); (R + (Pvec(tvec[i],a,b))[1])*sin(pvec[j]); (Pvec(tvec[i],a,b))[2] ];
            X[j+(i-1)*Np] = xshift .+ Rot * [tmp[1]*cos(tvec[i]); tmp[1]*sin(tvec[i]); tmp[2]]
        end
    end

    if ifplot
        figure(1); clf()
        scatter3D([x[1] for x in X], [x[2] for x in X], [x[3] for x in X], s=0.2)
        # scatter3D(X[:,1], X[:,2], X[:,3], s=0.2)
        xlim([-3;3])
        ylim([-3;3])
        zlim([-2;2])
        xlabel("x")
        ylabel("y")
        zlabel("z")
    end

    return X
end

if false # generating 8Mil points on surface
    X = toroid1(4000, 2000, ifplot=false)
    open(newdir*"/deformtorus.txt", "w") do io
        writedlm(io,X)
    end
    nx = length(X);
    dx = zeros(nx);
end

function Pgamma_deft(z)
    j=1; minv = 100.;
    @threads for i=1:nx
        dv = norm(z.-X[i])
        if dv<minv
            minv = dv
            j = i;
        end
    end
    return X[j]
end

function Pgamma_deft2(z)
    @threads for i=1:nx
        dx[i] = norm(z.-X[i])
    end
    j = argmin(dx)
    return X[j]
end

function Pgamma_deft3(z)
    dx = [norm(z.-x) for x in X]
    j = argmin(dx)
    return X[j]
end

# toroid1(100, 50, ifplot=true)
if false # checking fastest CPM
    M = 100;
    V = [rand(3).*[4;4;2].-[2;2;1] for m=1:M ];
    @time begin
        for m = 1:M
            x = V[m];
            y = Pgamma_deft(x)
            # plot3D([x[1];y[1]], [x[2];y[2]], [x[3];y[3]], "-r")
            # plot3D([y[1]], [y[2]], [y[3]], ".k")
            # plot3D([x[1]], [x[2]], [x[3]], "or")
        end
    end
    @time begin
        for m = 1:M
            x = V[m];
            y = Pgamma_deft2(x)
        end
    end
    @time begin
        for m = 1:M
            x = V[m];
            y = Pgamma_deft3(x)
        end
    end
end


Arot1(t) = [cos(2*t) sin(2*t);-sin(2*t) cos(2*t)];
Arot2(t) = [cos(2*t) -sin(2*t);sin(2*t) cos(2*t)];
Arot1b(cost,sint) = [cost^2-sint^2 2*sint*cost;-2*sint*cost cost^2-sint^2];
Arot2b(cost,sint) = [cost^2-sint^2 -2*sint*cost;2*sint*cost cost^2-sint^2];
RZfun(t,p,a,b,R) = [R;0] .+ Arot2(t)*[a 0;0 b]*Arot1(t)*[cos(p);sin(p)];
RZfun2(cost,sint,p,a,b,R) = [R;0] .+ Arot2b(cost,sint)*[a 0;0 b]*Arot1b(cost,sint)*[cos(p);sin(p)];

toroid1(50, 50, ifplot=true)

function Pgamma_deft4(z)
    zn = Rinv*(z.-xshift);
    phi = atan(zn[2],zn[1]);
    pvec = (0:0.001:1)*2*pi;
    Np = length(pvec);
    dv = zeros(Np)
    X = Array{Array{Float64,1},1}(undef,Np);
    @threads for j=1:Np
        R = 1.5; a = 0.5; b = 0.3;
        R = 0.8; a = 0.4; b = 0.3;
        tmp = RZfun2(cos(phi),sin(phi),pvec[j],a,b,R);
        X[j] = [tmp[1]*cos(phi); tmp[1]*sin(phi); tmp[2]]
        dv[j] = norm( zn - X[j])
    end
    return xshift .+ Rot*(X[argmin(dv)])
end

if false # checking fastest CPM
    M = 10000;
    V = [rand(3).*[4;4;2].-[2;2;1] for m=1:M ];
    # theta = rand()*2*pi;
    # phi = rand()*2*pi;
    # V = [ rand()*3*[sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)] for m=1:M ];
    # pvec = (0:0.001:1)*2*pi;
    # Np = length(pvec);
    # # dv = zeros(Np)
    # X = Array{Array{Float64,1},1}(undef,Np);
    # @threads for j=1:Np
    #     R = 1.5; a = 0.5; b = 0.3;
    #     tmp = RZfun2(cos(theta),sin(theta),pvec[j],a,b,R);
    #     X[j] = [tmp[1]*cos(theta); tmp[1]*sin(theta); tmp[2]]
    #     # dv[j] = norm( z - X[j])
    #     scatter3D(X[j][1], X[j][2], X[j][3], s=0.5, c="b")
    # end

    @time begin
        for m = 1:M
            x = V[m];
            y = Pgamma_deft4(x)
            # plot3D([x[1];y[1]], [x[2];y[2]], [x[3];y[3]], "-r")
            # plot3D([y[1]], [y[2]], [y[3]], ".k")
            # plot3D([x[1]], [x[2]], [x[3]], "or")
        end
    end
end

function insidepoint_deft(z)
    z = Rinv*(z.-xshift);
    phi = atan(z[2],z[1]);
    Rv = sqrt(z[1]^2+z[2]^2)
    Zv = z[3];
    R0 = 1.5; a = 0.5; b = 0.3;
    v = [a 0;0 b]*Arot2b(cos(phi),sin(phi))*[1/a 0;0 1/b]*Arot1b(cos(phi),sin(phi))*[Rv-R0;Zv];
    return 2*(v[1]^2/a^2+v[2]^2/b^2-1>0)-1
end

if false # checking fastest CPM
    toroid1(50, 50, ifplot=true)

    M = 40;
    # V = [rand(3).*[4;4;2].-[2;2;1] for m=1:M ];
    theta = rand()*2*pi;
    theta = 0.56*pi
    phi = rand()*2*pi;
    V = [ rand()*3*[sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)] for m=1:M ];
    # pvec = (0:0.001:1)*2*pi;
    # Np = length(pvec);
    # # dv = zeros(Np)
    # Xx = Array{Array{Float64,1},1}(undef,Np);
    # @threads for j=1:Np
    #     R = 1.5; a = 0.5; b = 0.3;
    #     tmp = RZfun2(cos(theta),sin(theta),pvec[j],a,b,R);
    #     Xx[j] = [tmp[1]*cos(theta); tmp[1]*sin(theta); tmp[2]]
    #     # dv[j] = norm( z - Xx[j])
    #     plot3D(Xx[j][1], Xx[j][2], Xx[j][3], ".b")
    # end

    @time begin
        for m = 1:M
            x = V[m];
            y = Pgamma_deft4(x)
            s = insidepoint_deft(x)
            colr = "b"
            if s>0
                colr = "r"
            end
            # colr = (s>0)*"r" + (s<0)*"b"
            plot3D([x[1];y[1]], [x[2];y[2]], [x[3];y[3]], "-r")
            # plot3D([y[1]], [y[2]], [y[3]], ".k")
            plot3D([x[1]], [x[2]], [x[3]], "o"*colr)
        end
    end
end

if ~false
    M = 3000;
    # V = [rand(3).*[-0.4;-0.25;0.5].+[0.2;0.1;0] for m=1:M ];
    # V = [rand(3).*0.5.*[1;1;1].-0.2*[1;1;1] for m=1:M ];
    V = [rand(3).*[2;2;2].-[1;1;1] for m=1:M ];
    W = [zeros(3) for i=1:M]
    for i=1:M
        # x = V[i]
        W[i] = Pgamma_deft4(V[i])
        # plot3D([y[1]], [y[2]], [y[3]], ".k")
    end
    figure(2); clf()
    scatter3D( [x[1] for x in W], [x[2] for x in W], [x[3] for x in W], s = 1 )
    # xlim([-1.2;1.2])
    # ylim([-1.2;1.2])
    # zlim([-0.5;0.5])
    # xlim([0;1.2])
    # ylim([0;1.2])
    # zlim([-0.5;0.5])
    xlim([-0.5;0])
    ylim([-0.5;0])
    zlim([0;0.5])
    # for i=1:M
    #     x = V[i]; y = W[i];
    #     plot3D([x[1];y[1]], [x[2];y[2]], [x[3];y[3]], "-") 
    # end
end

# toroid1(200, 100, ifplot=true)
