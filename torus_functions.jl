# rotation matrices along the three axes
Rx(a) = [1 0 0;0 cos(a) sin(a);0 -sin(a) cos(a)]
Ry(a) = [cos(a) 0 sin(a);0 1 0;-sin(a) 0 cos(a)]
Rz(a) = [cos(a) -sin(a) 0;sin(a) cos(a) 0;0 0 1]

# rotation parameters I usually use
a,b,c = 1.99487, 2.54097947651017, 4.219760487439292
# center of torus, origin is fine
xshift = zeros(3)

# rotation
Qrot = Rz(c)*Ry(b)*Rx(a);
# inverse rotation
Qinv = Rx(-a)*Ry(-b)*Rz(-c);

# R1 R2 radii of torus
R1 = 0.7; R2 = 0.2;

function Pgamma_torus(z,Qrot,Qinv,R1,R2,xshift; perm::Array{Int64,1}=[1;2;3], perm2::Array{Int64,1}=[1;2;3])
    zn = z[perm]
    zn = Qinv*(z.-xshift);
    phi = atan(zn[2],zn[1]);
    ctmp = R1*[cos(phi);sin(phi);0];
    return (Qrot*((zn.-ctmp)/norm(zn.-ctmp)*R2.+ctmp).+xshift)[perm2]
end

function insidepoint_torus(z,Qinv,R1,R2,xshift)
    zn = Qinv*(z.-xshift);
    phi = atan(zn[2],zn[1]);
    ctmp = R1*[cos(phi);sin(phi);0];
    return norm(ctmp-zn)<R2
end

function signed_distance(z,Qrot,Qinv,R1,R2,xshift)
    Pz = Pgamma_torus(z,Qrot,Qinv,R1,R2,xshift)
    zsign = insidepoint_torus(z,Qinv,R1,R2,xshift)
    return (1-2*zsign)*norm(z-Pz)
end
