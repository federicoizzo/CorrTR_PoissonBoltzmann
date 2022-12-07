include("corrected_TR_sub.jl")
include("corrected_TR_fun.jl")


if false # weight approximation test. 5 correct digits
  k1, k2 = rand(2)*4 .- 2
  Mmat = [k1 0;0 k2];
  η = rand()*0.2 -0.1
  D0mat = [1/(1-η*k1) 0; 0 1/(1-η*k2)];
  vtmp = rand(3); 
  vtmp2 = rand(3); 
  vtmp /= norm(vtmp)
  vtmp3 = cross(vtmp,vtmp2); vtmp3 /= norm(vtmp3)
  vtmp2 = cross(vtmp3,vtmp); vtmp2 /= norm(vtmp2)
  # Amat =  [-0.229457 0.973319; 0.973319 0.229457];
  Amat =  [vtmp[1] vtmp[2]; vtmp2[1] vtmp2[2]]
  epsl_ratio = 80*rand()
  f1(t) = i0_PB_a(t, D0mat, Amat, Mmat, epsl_ratio)
  f2(t) = j0_PB_a(t, D0mat, Amat, Mmat, 1/epsl_ratio)
  f3(t) = 1.0./ψ0_PB_a(t, D0mat, Amat)

  nab = 1000;
  w_a = zeros(3,nab)
  w_e = zeros(3,nab)
  abv = [rand(2).-0.5 for i=1:nab]

  @time for i=1:nab
    α, β = abv[i]
    w_a[1,i] = w_k0_ptilde1( [α;β]; lfun=f1 )
    w_a[2,i] = w_k0_ptilde1( [α;β]; lfun=f2 )
    w_a[3,i] = w_k0_ptilde1( [α;β]; lfun=f3 )
  end
  @time for i=1:nab
    α, β = abv[i]
    w_e[1,i] = weights2D_k0_single_TOL(α,β; lfun=f1)
    w_e[2,i] = weights2D_k0_single_TOL(α,β; lfun=f2)
    w_e[3,i] = weights2D_k0_single_TOL(α,β; lfun=f3)
  end
  println(norm(w_a[1,:].-w_e[1,:],Inf)," ",norm(w_a[1,:].-w_e[1,:])/sqrt(nab))
  println(norm(w_a[2,:].-w_e[2,:],Inf)," ",norm(w_a[2,:].-w_e[2,:])/sqrt(nab))
  println(norm(w_a[3,:].-w_e[3,:],Inf)," ",norm(w_a[3,:].-w_e[3,:])/sqrt(nab))
end

# dir = newdir
# nr = 4
# dir_ref = path*"/2022/2022-03/2022-03-15"
# nr_ref = 29
# surface = "torus"
# plotting_comparison_dir_ref(dir, nr, surface, dir_ref, nr_ref)

if false # tests for FD stencils
  f(x) = exp(-3x)+0.1*sin(4x)
  fp(x) = -3exp(-3x)+0.4*cos(4x)
  fpp(x) = 9exp(-3x)-1.6sin(4x)
  x0 = 0.124;
  h = [1.];
  rval = fp(x0)
  r2val = fpp(x0)
  rvec = zeros(10,9)
  hvec = zeros(10)
  for i=1:10
    h[1] *= 0.5
    hvec[i] = h[1]
    grid = x0.+(-2:3)*h[1]
    fgrid = f.(grid)
    rvec[i,1] = dot(vecval2nd_A[1,:],fgrid[3:5])/h[1]
    rvec[i,2] = dot(vecval2nd_A[2,:],fgrid[1:3])/h[1]
    rvec[i,3] = dot([-1;0;1]/2,fgrid[2:4])/h[1]
    rvec[i,4] = dot([-1/3;-1/2;1;-1/6],fgrid[2:5])/h[1]
    rvec[i,5] = dot([1/12;-2/3;0;2/3;-1/12],fgrid[1:5])/h[1]
    rvec[i,6] = dot(vecval_f2_4th_c,fgrid[1:5])/(h[1]^2)
    rvec[i,7] = dot(vecval_f2_2nd_c,fgrid[2:4])/(h[1]^2)
    rvec[i,8] = dot(vecval_f2_2nd_f,fgrid[3:6])/(h[1]^2)
    rvec[i,9] = dot(vecval_f2_p,fgrid[3:5])/(h[1]^2)
    # rvec[i,9] = dot(vecval_f2_2nd_b,fgrid[1:3])/(h[1]^2)
  end

  figure(2); clf()
  loglog(hvec,hvec.^2,"--")
  loglog(hvec,hvec.^3,"--")
  loglog(hvec,hvec.^4,"--")
  loglog(hvec,abs.(rvec[:,1].-rval),"-+")
  # loglog(hvec,abs.(rvec[:,2].-rval),"-x")
  # loglog(hvec,abs.(rvec[:,3].-rval),"-o")
  # loglog(hvec,abs.(rvec[:,4].-rval),"-s")
  # loglog(hvec,abs.(rvec[:,5].-rval),"-^")
  loglog(hvec,abs.(rvec[:,6].-r2val),"-v")
  loglog(hvec,abs.(rvec[:,7].-r2val),"->")
  loglog(hvec,abs.(rvec[:,8].-r2val),"-+")
  loglog(hvec,abs.(rvec[:,9].-r2val),"-x")
end

if false  # ε~2h test on the spheres 
  nrunf!(nrun)

  Ttarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_t.dat")[:,1]
  Ptarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_p.dat")[:,1]
  
  TPtargets = [ [Ttarget[i];Ptarget[i];(0.5*Ttarget[i]-0.12)] for i=1:length(Ttarget) ]

  Ntargets = length(Ttarget)
  println("Number of targets = $(Ntargets)")

  x = TPtargets;
  open(newdir_nrun*"/torus_err_const_$(nrun[1])_trg_xyz.dat","w") do io
    writedlm(io,x)
  end

  # Nvec = 64:10:404
  # Nvec = 64:10:304
  # Nvec = 64:10:264
  # Nvec = 64:8:400
  Nvec = 64:8:300
  # Nvec = 64:8:248
  # Nvec = 64:10:224
  
  # fε(x) = (0.30*x.^0.5)
  # fε(x) = (2*x.^0.8)
  # fε(x) = 4x
  fε(x) = 2x
  # detail = "epsl03h05"
  # detail = "epsl2h08"
  # detail = "epsl4h"
  detail = "epsl2h"

  nmax = length(Nvec)
  hval = zeros(nmax)
  mval_abs = zeros(8,nmax)
  mval_err = zeros(4,nmax)
  surf_val= zeros(2,nmax)
  
  val_abs = zeros(8,Ntargets,nmax)
  val_err = zeros(4,Ntargets,nmax)

  shift=[0.026834;0.01346254;-0.06293856]

  @time for i=1:nmax
      println("\nRun $i/$nmax")
      local h = 0.1
      ε = 0.1
      @time hval[i], val_abs[:,:,i], val_err[:,:,i], mval_abs[:,i], mval_err[:,i],  surf_val[:,i] = PB_gen_shape_WENO( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i)
      GC.gc()
  end
  x = detail;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_detail.dat","w") do io
    writedlm(io,x)
  end
  x = hval;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_hvals.dat","w") do io
    writedlm(io,x)
  end
  x = val_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs.dat","w") do io
    writedlm(io,x)
  end
  x = val_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_abs);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_size.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_err);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_size.dat","w") do io
    writedlm(io,x)
  end
  x = mval_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs_means.dat","w") do io
    writedlm(io,x)
  end
  x = mval_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_means.dat","w") do io
    writedlm(io,x)
  end
  x = surf_val;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val.dat","w") do io
    writedlm(io,x)
  end
  x = size(surf_val);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val_size.dat","w") do io
    writedlm(io,x)
  end
end
# plotting_comparison_1(hval, val_abs, surf_val, newdir, nrun[1], "spheres", detail)
# plotting_comparison_1(hval, val_abs, surf_val, newdir, 2, "spheres", detail)
if false  # ε~4h test on the spheres
  nrunf!(nrun)

  Ttarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_t.dat")[:,1]
  Ptarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_p.dat")[:,1]
  
  TPtargets = [ [Ttarget[i];Ptarget[i];(0.5*Ttarget[i]-0.12)] for i=1:length(Ttarget) ]

  Ntargets = length(Ttarget)
  println("Number of targets = $(Ntargets)")

  x = TPtargets;
  open(newdir_nrun*"/torus_err_const_$(nrun[1])_trg_xyz.dat","w") do io
    writedlm(io,x)
  end

  # Nvec = 64:10:404
  # Nvec = 64:10:304
  # Nvec = 64:10:264
  # Nvec = 64:8:400
  Nvec = 64:8:300
  # Nvec = 64:8:248
  # Nvec = 64:10:224
  
  # fε(x) = (0.30*x.^0.5)
  # fε(x) = (2*x.^0.8)
  fε(x) = 4x
  # fε(x) = 2x
  # detail = "epsl03h05"
  # detail = "epsl2h08"
  detail = "epsl4h"
  # detail = "epsl2h"

  nmax = length(Nvec)
  hval = zeros(nmax)
  mval_abs = zeros(8,nmax)
  mval_err = zeros(4,nmax)
  surf_val= zeros(8,nmax)
  
  val_abs = zeros(8,Ntargets,nmax)
  val_err = zeros(4,Ntargets,nmax)

  shift=[0.026834;0.01346254;-0.06293856]

  @time for i=1:nmax
      println("\nRun $i/$nmax")
      local h = 0.1
      ε = 0.1
      @time hval[i], val_abs[:,:,i], val_err[:,:,i], mval_abs[:,i], mval_err[:,i],  surf_val[:,i] = PB_gen_shape_WENO( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i)
      GC.gc()
  end
  x = detail;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_detail.dat","w") do io
    writedlm(io,x)
  end
  x = hval;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_hvals.dat","w") do io
    writedlm(io,x)
  end
  x = val_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs.dat","w") do io
    writedlm(io,x)
  end
  x = val_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_abs);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_size.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_err);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_size.dat","w") do io
    writedlm(io,x)
  end
  x = mval_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs_means.dat","w") do io
    writedlm(io,x)
  end
  x = mval_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_means.dat","w") do io
    writedlm(io,x)
  end
  x = surf_val;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val.dat","w") do io
    writedlm(io,x)
  end
end
# plotting_comparison_1(hval, val_abs, surf_val, newdir, nrun[1], "spheres", detail)
if false  # ε~h^{0.8} test on the spheres
  nrunf!(nrun)

  Ttarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_t.dat")[:,1]
  Ptarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_p.dat")[:,1]
  
  TPtargets = [ [Ttarget[i];Ptarget[i];(0.5*Ttarget[i]-0.12)] for i=1:length(Ttarget) ]

  Ntargets = length(Ttarget)
  println("Number of targets = $(Ntargets)")

  x = TPtargets;
  open(newdir_nrun*"/torus_err_const_$(nrun[1])_trg_xyz.dat","w") do io
    writedlm(io,x)
  end

  # Nvec = 64:10:404
  # Nvec = 64:10:304
  # Nvec = 64:10:264
  # Nvec = 64:8:400
  Nvec = 64:8:300
  # Nvec = 64:8:248
  # Nvec = 64:10:224
  
  # fε(x) = (0.30*x.^0.5)
  fε(x) = (2*x.^0.8)
  # fε(x) = 4x
  # fε(x) = 2x
  # detail = "epsl03h05"
  detail = "epsl2h08"
  # detail = "epsl4h"
  # detail = "epsl2h"

  nmax = length(Nvec)
  hval = zeros(nmax)
  mval_abs = zeros(8,nmax)
  mval_err = zeros(4,nmax)
  surf_val= zeros(8,nmax)
  
  val_abs = zeros(8,Ntargets,nmax)
  val_err = zeros(4,Ntargets,nmax)

  shift=[0.026834;0.01346254;-0.06293856]

  @time for i=1:nmax
      println("\nRun $i/$nmax")
      local h = 0.1
      ε = 0.1
      @time hval[i], val_abs[:,:,i], val_err[:,:,i], mval_abs[:,i], mval_err[:,i],  surf_val[:,i] = PB_gen_shape_WENO( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i)
      GC.gc()
  end
  x = detail;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_detail.dat","w") do io
    writedlm(io,x)
  end
  x = hval;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_hvals.dat","w") do io
    writedlm(io,x)
  end
  x = val_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs.dat","w") do io
    writedlm(io,x)
  end
  x = val_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_abs);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_size.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_err);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_size.dat","w") do io
    writedlm(io,x)
  end
  x = mval_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs_means.dat","w") do io
    writedlm(io,x)
  end
  x = mval_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_means.dat","w") do io
    writedlm(io,x)
  end
  x = surf_val;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val.dat","w") do io
    writedlm(io,x)
  end
end
# plotting_comparison_1(hval, val_abs, surf_val, newdir, nrun[1], "spheres", detail)
if false  # ε~h^{0.5} test on the spheres
  nrunf!(nrun)

  Ttarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_t.dat")[:,1]
  Ptarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_p.dat")[:,1]
  
  TPtargets = [ [Ttarget[i];Ptarget[i];(0.5*Ttarget[i]-0.12)] for i=1:length(Ttarget) ]

  Ntargets = length(Ttarget)
  println("Number of targets = $(Ntargets)")

  x = TPtargets;
  open(newdir_nrun*"/torus_err_const_$(nrun[1])_trg_xyz.dat","w") do io
    writedlm(io,x)
  end

  # Nvec = 64:10:404
  # Nvec = 64:10:304
  # Nvec = 64:10:264
  # Nvec = 64:8:400
  Nvec = 64:8:300
  # Nvec = 64:8:248
  # Nvec = 64:10:224
  
  fε(x) = (0.30*x.^0.5)
  # fε(x) = (2*x.^0.8)
  # fε(x) = 4x
  # fε(x) = 2x
  detail = "epsl03h05"
  # detail = "epsl2h08"
  # detail = "epsl4h"
  # detail = "epsl2h"

  nmax = length(Nvec)
  hval = zeros(nmax)
  mval_abs = zeros(8,nmax)
  mval_err = zeros(4,nmax)
  surf_val= zeros(8,nmax)
  
  val_abs = zeros(8,Ntargets,nmax)
  val_err = zeros(4,Ntargets,nmax)

  shift=[0.026834;0.01346254;-0.06293856]

  @time for i=1:nmax
      println("\nRun $i/$nmax")
      local h = 0.1
      ε = 0.1
      @time hval[i], val_abs[:,:,i], val_err[:,:,i], mval_abs[:,i], mval_err[:,i],  surf_val[:,i] = PB_gen_shape_WENO( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i)
      GC.gc()
  end
  x = detail;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_detail.dat","w") do io
    writedlm(io,x)
  end
  x = hval;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_hvals.dat","w") do io
    writedlm(io,x)
  end
  x = val_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs.dat","w") do io
    writedlm(io,x)
  end
  x = val_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_abs);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_size.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_err);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_size.dat","w") do io
    writedlm(io,x)
  end
  x = mval_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs_means.dat","w") do io
    writedlm(io,x)
  end
  x = mval_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_means.dat","w") do io
    writedlm(io,x)
  end
  x = surf_val;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val.dat","w") do io
    writedlm(io,x)
  end
end
# plotting_comparison_1(hval, val_abs, surf_val, newdir, nrun[1], "spheres", detail)

# solving system on sphere
if false  # ε~5h test on the surfaces
  nrunf!(nrun)

  Ttarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_t.dat")[:,1]
  Ptarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_p.dat")[:,1]
  
  TPtargets = [ [Ttarget[i];Ptarget[i];(0.5*Ttarget[i]-0.12)] for i=1:length(Ttarget) ]

  Ntargets = length(Ttarget)
  println("Number of targets = $(Ntargets)")

  x = TPtargets;
  open(newdir_nrun*"/torus_err_const_$(nrun[1])_trg_xyz.dat","w") do io
    writedlm(io,x)
  end

  # Nvec = 64:10:404
  # Nvec = 64:10:304
  # Nvec = 64:10:264
  # Nvec = 64:8:400
  # Nvec = 64:8:300
  # Nvec = [40;64;80]
  # Nvec = [23]
  Nvec = [20;40;80]
  # Nvec = [25;50;100]
  # Nvec = 64:8:248
  # Nvec = 64:10:224
  
  # fε(x) = (0.30*x.^0.5)
  # fε(x) = (2*x.^0.8)
  # fε(x) = 4x
  fε(x) = 5x
  # detail = "epsl03h05"
  # detail = "epsl2h08"
  # detail = "epsl4h"
  detail = "epsl5h"

  nmax = length(Nvec)
  errvec = zeros(nmax,2,3)
  hval = zeros(nmax)
  mval_abs = zeros(8,nmax)
  mval_err = zeros(4,nmax)
  surf_val= zeros(2,nmax)
  
  val_abs = zeros(8,Ntargets,nmax)
  val_err = zeros(4,Ntargets,nmax)

  shift=[0.026834;0.01346254;-0.06293856]

  @time for i=1:nmax
      println("\nRun $i/$nmax")
      local h = 0.1
      ε = 0.1
      # @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-0.75, Zlim2=0.75, zcenter=zeros(3) )

      @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-0.75, Zlim2=0.75, zcenter=zeros(3) )
      # @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-1.5, Zlim2=1.0, zcenter=zeros(3) )
      GC.gc()
  end
  x = detail;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_detail.dat","w") do io
    writedlm(io,x)
  end
  x = hval;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_hvals.dat","w") do io
    writedlm(io,x)
  end
  x = errvec;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_errvec.dat","w") do io
    writedlm(io,x)
  end
  x = val_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs.dat","w") do io
    writedlm(io,x)
  end
  x = val_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_abs);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_size.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_err);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_size.dat","w") do io
    writedlm(io,x)
  end
  x = mval_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs_means.dat","w") do io
    writedlm(io,x)
  end
  x = mval_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_means.dat","w") do io
    writedlm(io,x)
  end
  x = surf_val;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val.dat","w") do io
    writedlm(io,x)
  end
  x = size(surf_val);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val_size.dat","w") do io
    writedlm(io,x)
  end
end
if false  # ε~2h test on the surfaces
  nrunf!(nrun)

  Ttarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_t.dat")[:,1]
  Ptarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_p.dat")[:,1]
  
  TPtargets = [ [Ttarget[i];Ptarget[i];(0.5*Ttarget[i]-0.12)] for i=1:length(Ttarget) ]

  Ntargets = length(Ttarget)
  println("Number of targets = $(Ntargets)")

  x = TPtargets;
  open(newdir_nrun*"/torus_err_const_$(nrun[1])_trg_xyz.dat","w") do io
    writedlm(io,x)
  end

  # Nvec = 64:10:404
  # Nvec = 64:10:304
  # Nvec = 64:10:264
  # Nvec = 64:8:400
  # Nvec = 64:8:300
  # Nvec = [40;64;80]
  # Nvec = [23]
  Nvec = [20;40;80]
  # Nvec = 64:8:248
  # Nvec = 64:10:224
  
  # fε(x) = (0.30*x.^0.5)
  # fε(x) = (2*x.^0.8)
  # fε(x) = 4x
  fε(x) = 2x
  # detail = "epsl03h05"
  # detail = "epsl2h08"
  # detail = "epsl4h"
  detail = "epsl2h"

  nmax = length(Nvec)
  errvec = zeros(nmax,2,3)
  hval = zeros(nmax)
  mval_abs = zeros(8,nmax)
  mval_err = zeros(4,nmax)
  surf_val= zeros(2,nmax)
  
  val_abs = zeros(8,Ntargets,nmax)
  val_err = zeros(4,Ntargets,nmax)

  shift=[0.026834;0.01346254;-0.06293856]

  @time for i=1:nmax
      println("\nRun $i/$nmax")
      local h = 0.1
      ε = 0.1
      # @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-0.75, Zlim2=0.75, zcenter=zeros(3) )
      @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-1.5, Zlim2=1.0, zcenter=zeros(3) )
      GC.gc()
  end
  x = detail;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_detail.dat","w") do io
    writedlm(io,x)
  end
  x = hval;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_hvals.dat","w") do io
    writedlm(io,x)
  end
  x = errvec;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_errvec.dat","w") do io
    writedlm(io,x)
  end
  x = val_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs.dat","w") do io
    writedlm(io,x)
  end
  x = val_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_abs);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_size.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_err);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_size.dat","w") do io
    writedlm(io,x)
  end
  x = mval_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs_means.dat","w") do io
    writedlm(io,x)
  end
  x = mval_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_means.dat","w") do io
    writedlm(io,x)
  end
  x = surf_val;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val.dat","w") do io
    writedlm(io,x)
  end
  x = size(surf_val);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val_size.dat","w") do io
    writedlm(io,x)
  end
end

# solving system on spheres

if false  # ε~2h test on the surfaces
  nrunf!(nrun)

  Ttarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_t.dat")[:,1]
  Ptarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_p.dat")[:,1]
  
  TPtargets = [ [Ttarget[i];Ptarget[i];(0.5*Ttarget[i]-0.12)] for i=1:length(Ttarget) ]

  Ntargets = length(Ttarget)
  println("Number of targets = $(Ntargets)")

  x = TPtargets;
  open(newdir_nrun*"/torus_err_const_$(nrun[1])_trg_xyz.dat","w") do io
    writedlm(io,x)
  end

  # Nvec = 64:10:404
  # Nvec = 64:10:304
  # Nvec = 64:10:264
  # Nvec = 64:8:400
  # Nvec = 64:8:300
  # Nvec = [40;64;80]
  # Nvec = [23]
  Nvec = [20;40;80]
  # Hvec = 0.075*(2.0.^(0:-1:-2))
  # Hvec = 0.1*(2.0.^(0:-1:-2))
  # Hvec = 0.115*(2.0.^(0:-1:-2))

  Hvec = 0.075*(2.0.^(0:-1:-2))
  # Nvec = [25;50;100]
  # Nvec = 64:8:248
  # Nvec = 64:10:224
  
  # fε(x) = (0.30*x.^0.5)
  # fε(x) = (2*x.^0.8)
  # fε(x) = 4x
  fε(x) = 2x
  # detail = "epsl03h05"
  # detail = "epsl2h08"
  # detail = "epsl4h"
  detail = "epsl2h"

  nmax = length(Nvec)
  errvec = zeros(nmax,2,3)
  hval = zeros(nmax)
  mval_abs = zeros(8,nmax)
  mval_err = zeros(4,nmax)
  surf_val= zeros(2,nmax)
  
  val_abs = zeros(8,Ntargets,nmax)
  val_err = zeros(4,Ntargets,nmax)

  shift=[0.026834;0.01346254;-0.06293856]

  @time for i=1:nmax
      println("\nRun $i/$nmax")
      local h = 0.1
      ε = 0.1
      # @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-0.75, Zlim2=0.75, zcenter=zeros(3) )
      # @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Hvec[i], Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-.75, Zlim2=.75, zcenter=zeros(3), hflag=true )
      @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Hvec[i], Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-1.2, Zlim2=1.2, zcenter=zeros(3), hflag=true )
      GC.gc()
  end
  x = detail;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_detail.dat","w") do io
    writedlm(io,x)
  end
  x = hval;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_hvals.dat","w") do io
    writedlm(io,x)
  end
  x = errvec;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_errvec.dat","w") do io
    writedlm(io,x)
  end
  x = val_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs.dat","w") do io
    writedlm(io,x)
  end
  x = val_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_abs);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_size.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_err);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_size.dat","w") do io
    writedlm(io,x)
  end
  x = mval_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs_means.dat","w") do io
    writedlm(io,x)
  end
  x = mval_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_means.dat","w") do io
    writedlm(io,x)
  end
  x = surf_val;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val.dat","w") do io
    writedlm(io,x)
  end
  x = size(surf_val);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val_size.dat","w") do io
    writedlm(io,x)
  end
end
if false  # ε~4h test on the surfaces
  nrunf!(nrun)

  Ttarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_t.dat")[:,1]
  Ptarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_p.dat")[:,1]
  
  TPtargets = [ [Ttarget[i];Ptarget[i];(0.5*Ttarget[i]-0.12)] for i=1:length(Ttarget) ]

  Ntargets = length(Ttarget)
  println("Number of targets = $(Ntargets)")

  x = TPtargets;
  open(newdir_nrun*"/torus_err_const_$(nrun[1])_trg_xyz.dat","w") do io
    writedlm(io,x)
  end

  # Nvec = 64:10:404
  # Nvec = 64:10:304
  # Nvec = 64:10:264
  # Nvec = 64:8:400
  # Nvec = 64:8:300
  # Nvec = [40;64;80]
  # Nvec = [23]
  Nvec = [20;40;80]
  # Hvec = 0.075*(2.0.^(0:-1:-2))
  # Hvec = 0.1*(2.0.^(0:-1:-2))
  # Hvec = 0.1*(2.0.^(0:-1:-1))
  Hvec = 0.075*(2.0.^(0:-1:-2))
  # Nvec = [25;50;100]
  # Nvec = 64:8:248
  # Nvec = 64:10:224
  
  # fε(x) = (0.30*x.^0.5)
  # fε(x) = (2*x.^0.8)
  fε(x) = 4x
  # fε(x) = 2x
  # detail = "epsl03h05"
  # detail = "epsl2h08"
  detail = "epsl4h"
  # detail = "epsl2h"

  nmax = length(Nvec)
  errvec = zeros(nmax,2,3)
  hval = zeros(nmax)
  mval_abs = zeros(8,nmax)
  mval_err = zeros(4,nmax)
  surf_val= zeros(2,nmax)
  
  val_abs = zeros(8,Ntargets,nmax)
  val_err = zeros(4,Ntargets,nmax)

  shift=[0.026834;0.01346254;-0.06293856]

  @time for i=1:nmax
      println("\nRun $i/$nmax")
      local h = 0.1
      ε = 0.1
      # @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-0.75, Zlim2=0.75, zcenter=zeros(3) )
      # @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Hvec[i], Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-.5, Zlim2=.5, zcenter=zeros(3), hflag=true )
      @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Hvec[i], Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-1., Zlim2=0.5, zcenter=[0.;0.25;0.4], hflag=true )
      GC.gc()
  end
  x = detail;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_detail.dat","w") do io
    writedlm(io,x)
  end
  x = hval;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_hvals.dat","w") do io
    writedlm(io,x)
  end
  x = errvec;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_errvec.dat","w") do io
    writedlm(io,x)
  end
  x = val_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs.dat","w") do io
    writedlm(io,x)
  end
  x = val_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_abs);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_size.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_err);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_size.dat","w") do io
    writedlm(io,x)
  end
  x = mval_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs_means.dat","w") do io
    writedlm(io,x)
  end
  x = mval_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_means.dat","w") do io
    writedlm(io,x)
  end
  x = surf_val;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val.dat","w") do io
    writedlm(io,x)
  end
  x = size(surf_val);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val_size.dat","w") do io
    writedlm(io,x)
  end
end


if false  # J=1, ε~4h test on the surfaces
  nrunf!(nrun)

  Ttarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_t.dat")[:,1]
  Ptarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_p.dat")[:,1]
  
  TPtargets = [ [Ttarget[i];Ptarget[i];(0.5*Ttarget[i]-0.12)] for i=1:length(Ttarget) ]

  Ntargets = length(Ttarget)
  println("Number of targets = $(Ntargets)")

  x = TPtargets;
  open(newdir_nrun*"/torus_err_const_$(nrun[1])_trg_xyz.dat","w") do io
    writedlm(io,x)
  end

  # Nvec = 64:10:404
  # Nvec = 64:10:304
  # Nvec = 64:10:264
  # Nvec = 64:8:400
  # Nvec = 64:8:300
  # Nvec = [40;64;80]
  # Nvec = [23]
  Nvec = [20;40;80]
  Hvec = 0.075*(2.0.^(0:-1:-2))
  # Nvec = [25;50;100]
  # Nvec = 64:8:248
  # Nvec = 64:10:224
  
  # fε(x) = (0.30*x.^0.5)
  # fε(x) = (2*x.^0.8)
  fε(x) = 4x
  # fε(x) = 2x
  # detail = "epsl03h05"
  # detail = "epsl2h08"
  detail = "epsl4h"
  # detail = "epsl2h"

  nmax = length(Nvec)
  errvec = zeros(nmax,2,3)
  hval = zeros(nmax)
  mval_abs = zeros(8,nmax)
  mval_err = zeros(4,nmax)
  surf_val= zeros(2,nmax)
  
  val_abs = zeros(8,Ntargets,nmax)
  val_err = zeros(4,Ntargets,nmax)

  shift=[0.026834;0.01346254;-0.06293856]

  @time for i=1:nmax
      println("\nRun $i/$nmax")
      local h = 0.1
      ε = 0.1
      # @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-0.75, Zlim2=0.75, zcenter=zeros(3) )
      @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Hvec[i], Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-.75, Zlim2=.75, zcenter=zeros(3), hflag=true, Jflag=true )
      GC.gc()
  end
  x = detail;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_detail.dat","w") do io
    writedlm(io,x)
  end
  x = hval;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_hvals.dat","w") do io
    writedlm(io,x)
  end
  x = errvec;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_errvec.dat","w") do io
    writedlm(io,x)
  end
  x = val_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs.dat","w") do io
    writedlm(io,x)
  end
  x = val_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_abs);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_size.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_err);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_size.dat","w") do io
    writedlm(io,x)
  end
  x = mval_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs_means.dat","w") do io
    writedlm(io,x)
  end
  x = mval_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_means.dat","w") do io
    writedlm(io,x)
  end
  x = surf_val;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val.dat","w") do io
    writedlm(io,x)
  end
  x = size(surf_val);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val_size.dat","w") do io
    writedlm(io,x)
  end
end
if false  # J=1, ε~2h test on the surfaces
  nrunf!(nrun)

  Ttarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_t.dat")[:,1]
  Ptarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_p.dat")[:,1]
  
  TPtargets = [ [Ttarget[i];Ptarget[i];(0.5*Ttarget[i]-0.12)] for i=1:length(Ttarget) ]

  Ntargets = length(Ttarget)
  println("Number of targets = $(Ntargets)")

  x = TPtargets;
  open(newdir_nrun*"/torus_err_const_$(nrun[1])_trg_xyz.dat","w") do io
    writedlm(io,x)
  end

  # Nvec = 64:10:404
  # Nvec = 64:10:304
  # Nvec = 64:10:264
  # Nvec = 64:8:400
  # Nvec = 64:8:300
  # Nvec = [40;64;80]
  # Nvec = [23]
  Nvec = [20;40;80]
  Hvec = 0.075*(2.0.^(0:-1:-2))
  # Nvec = [25;50;100]
  # Nvec = 64:8:248
  # Nvec = 64:10:224
  
  # fε(x) = (0.30*x.^0.5)
  # fε(x) = (2*x.^0.8)
  # fε(x) = 4x
  fε(x) = 2x
  # detail = "epsl03h05"
  # detail = "epsl2h08"
  # detail = "epsl4h"
  detail = "epsl2h"

  nmax = length(Nvec)
  errvec = zeros(nmax,2,3)
  hval = zeros(nmax)
  mval_abs = zeros(8,nmax)
  mval_err = zeros(4,nmax)
  surf_val= zeros(2,nmax)
  
  val_abs = zeros(8,Ntargets,nmax)
  val_err = zeros(4,Ntargets,nmax)

  shift=[0.026834;0.01346254;-0.06293856]

  @time for i=1:nmax
      println("\nRun $i/$nmax")
      local h = 0.1
      ε = 0.1
      # @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-0.75, Zlim2=0.75, zcenter=zeros(3) )
      @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Hvec[i], Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-.75, Zlim2=.75, zcenter=zeros(3), hflag=true, Jflag=true )
      GC.gc()
  end
  x = detail;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_detail.dat","w") do io
    writedlm(io,x)
  end
  x = hval;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_hvals.dat","w") do io
    writedlm(io,x)
  end
  x = errvec;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_errvec.dat","w") do io
    writedlm(io,x)
  end
  x = val_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs.dat","w") do io
    writedlm(io,x)
  end
  x = val_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_abs);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_size.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_err);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_size.dat","w") do io
    writedlm(io,x)
  end
  x = mval_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs_means.dat","w") do io
    writedlm(io,x)
  end
  x = mval_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_means.dat","w") do io
    writedlm(io,x)
  end
  x = surf_val;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val.dat","w") do io
    writedlm(io,x)
  end
  x = size(surf_val);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val_size.dat","w") do io
    writedlm(io,x)
  end
end

# test for single plane, checking curvatures
if false
  nrunf!(nrun)

  Ttarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_t.dat")[:,1]
  Ptarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_p.dat")[:,1]
  
  TPtargets = [ [Ttarget[i];Ptarget[i];(0.5*Ttarget[i]-0.12)] for i=1:length(Ttarget) ]

  Ntargets = length(Ttarget)
  println("Number of targets = $(Ntargets)")

  x = TPtargets;
  open(newdir_nrun*"/torus_err_const_$(nrun[1])_trg_xyz.dat","w") do io
    writedlm(io,x)
  end

  # Nvec = 64:10:404
  # Nvec = 64:10:304
  # Nvec = 64:10:264
  # Nvec = 64:8:400
  # Nvec = 64:8:300
  # Nvec = [40;64;80]
  # Nvec = [23]
  # Nvec = [20;40;80]
  Nvec = [80]
  # Hvec = 0.075*(2.0.^(0:-1:-2))
  # Hvec = 0.1*(2.0.^(0:-1:-2))
  # Hvec = 0.1*(2.0.^(0:-1:-1))
  Hvec = 0.04*(2.0.^(0:-1:-1))
  # Nvec = [25;50;100]
  # Nvec = 64:8:248
  # Nvec = 64:10:224
  
  # fε(x) = (0.30*x.^0.5)
  # fε(x) = (2*x.^0.8)
  fε(x) = 4x
  # fε(x) = 2x
  # detail = "epsl03h05"
  # detail = "epsl2h08"
  detail = "epsl4h"
  # detail = "epsl2h"

  nmax = length(Nvec)
  errvec = zeros(nmax,2,3)
  hval = zeros(nmax)
  mval_abs = zeros(8,nmax)
  mval_err = zeros(4,nmax)
  surf_val= zeros(2,nmax)
  
  val_abs = zeros(8,Ntargets,nmax)
  val_err = zeros(4,Ntargets,nmax)

  shift=[0.026834;0.01346254;-0.06293856]

  @time for i=1:nmax
      println("\nRun $i/$nmax")
      local h = 0.1
      ε = 0.1
      # @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-0.75, Zlim2=0.75, zcenter=zeros(3) )
      # @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Hvec[i], Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-.5, Zlim2=.5, zcenter=zeros(3), hflag=true )
      @time _ = PB_gen_shape_plane( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==1), count=i)#, Zlim1=-1., Zlim2=0.5, zcenter=[0.;0.25;0.4], hflag=true )
      GC.gc()
  end
  # x = detail;
  # open(newdir_nrun*"/surf_PB_$(nrun[1])_detail.dat","w") do io
  #   writedlm(io,x)
  # end
  # x = hval;
  # open(newdir_nrun*"/surf_PB_$(nrun[1])_hvals.dat","w") do io
  #   writedlm(io,x)
  # end
  # x = errvec;
  # open(newdir_nrun*"/surf_PB_$(nrun[1])_errvec.dat","w") do io
  #   writedlm(io,x)
  # end
  # x = val_abs;
  # open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs.dat","w") do io
  #   writedlm(io,x)
  # end
  # x = val_err;
  # open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err.dat","w") do io
  #   writedlm(io,x)
  # end
  # x = size(val_abs);
  # open(newdir_nrun*"/surf_PB_$(nrun[1])_val_size.dat","w") do io
  #   writedlm(io,x)
  # end
  # x = size(val_err);
  # open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_size.dat","w") do io
  #   writedlm(io,x)
  # end
  # x = mval_abs;
  # open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs_means.dat","w") do io
  #   writedlm(io,x)
  # end
  # x = mval_err;
  # open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_means.dat","w") do io
  #   writedlm(io,x)
  # end
  # x = surf_val;
  # open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val.dat","w") do io
  #   writedlm(io,x)
  # end
  # x = size(surf_val);
  # open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val_size.dat","w") do io
  #   writedlm(io,x)
  # end
end

if false # tests checking number of bad curvature points as h->0
  nrunf!(nrun)

  Ttarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_t.dat")[:,1]
  Ptarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_p.dat")[:,1]
  
  TPtargets = [ [Ttarget[i];Ptarget[i];(0.5*Ttarget[i]-0.12)] for i=1:length(Ttarget) ]

  Ntargets = length(Ttarget)
  println("Number of targets = $(Ntargets)")

  x = TPtargets;
  open(newdir_nrun*"/torus_err_const_$(nrun[1])_trg_xyz.dat","w") do io
    writedlm(io,x)
  end

  Nvec = 20:10:400
  # Nvec = 64:10:304
  # Nvec = 64:10:264
  # Nvec = 64:8:400
  # Nvec = 64:8:300
  # Nvec = [40;64;80]
  # Nvec = [23]
  # Nvec = [20;40;80]
  # Nvec = [80]
  # Hvec = 0.075*(2.0.^(0:-1:-2))
  # Hvec = 0.1*(2.0.^(0:-1:-2))
  # Hvec = 0.1*(2.0.^(0:-1:-1))
  Hvec = 0.04*(2.0.^(0:-1:-1))
  # Nvec = [25;50;100]
  # Nvec = 64:8:248
  # Nvec = 64:10:224
  
  # fε(x) = (0.30*x.^0.5)
  # fε(x) = (2*x.^0.8)
  # fε(x) = 0.25
  fε(x) = 4x
  # fε(x) = 2x
  # detail = "epsl03h05"
  # detail = "epsl2h08"
  detail = "epsl4h"
  # detail = "epsl2h"

  nmax = length(Nvec)
  errvec = zeros(nmax,2,3)
  hval = zeros(nmax)
  mval_abs = zeros(8,nmax)
  mval_err = zeros(4,nmax)
  surf_val= zeros(2,nmax)
  
  val_abs = zeros(8,Ntargets,nmax)
  val_err = zeros(4,Ntargets,nmax)

  shift=[0.026834;0.01346254;-0.06293856]
  Mvec = zeros(nmax,2)
  pvec = zeros(nmax)

  @time for i=1:nmax
      println("\nRun $i/$nmax")
      local h = 0.1
      ε = 0.1
      # @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-0.75, Zlim2=0.75, zcenter=zeros(3) )
      # @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Hvec[i], Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-.5, Zlim2=.5, zcenter=zeros(3), hflag=true )
      @time hval[i], _, _, _, _, _, Mvec[i,:], pvec[i] = PB_gen_shape_bad_kappas( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i)#, Zlim1=-1., Zlim2=0.5, zcenter=[0.;0.25;0.4], hflag=true )
      GC.gc()
  end
end
if false
  Mvec = Mvec[hval.>0,:]
  pvec = pvec[hval.>0]
  hval = hval[hval.>0]  
  x = hval;
  open(newdir_nrun*"/surf_PB_badpoints_$(nrun[1])_hvals.dat","w") do io
    writedlm(io,x)
  end
  x = Mvec;
  open(newdir_nrun*"/surf_PB_badpoints_$(nrun[1])_Mvec.dat","w") do io
    writedlm(io,x)
  end
  x = pvec;
  open(newdir_nrun*"/surf_PB_badpoints_$(nrun[1])_pvec.dat","w") do io
    writedlm(io,x)
  end

  figure(90); 
  # clf()
  subplot(1,2,1)
  # semilogx(hval, pvec, label="epsl=0.25")
  semilogx(hval, pvec, label="epsl=4h")
  legend(loc="best")
  ylabel("% of bad points")
  xlabel(L"h")
  title("Semilogx")
  subplot(1,2,2)
  # loglog(hval, pvec, label="epsl=0.25")
  loglog(hval, pvec, label="epsl=4h")
  legend(loc="best")
  ylabel("% of bad points")
  xlabel(L"h")
  title("Loglog")
  savefig(newdir_nrun*"/surf_PB_badpoints.png")
end
# epsl = 0.25, 2022-07-06/38
# epsl = 4h, 2022-07-06/1
# 
if false
  figure(90); clf()
  subplot(1,2,1)
  loglog(hval1, pvec1, "-", label=L"\varepsilon=0.25")
  loglog(hval2, pvec2, "-.", label=L"\varepsilon=4h")
  loglog(hval2, 100*hval2, "--", label=L"O(h)")
  legend(loc="best")
  # ylabel("% of bad points")
  xlabel(L"h")
  title("Percentage of bad points")
  subplot(1,2,2)
  loglog(hval1, Mvec1[:,1], "-", label=L"\varepsilon=0.25")
  loglog(hval2, Mvec2[:,1], "-.", label=L"\varepsilon=4h")
  loglog(hval2, 6hval2.^(-3), ":", label=L"O(h^{-3})")
  loglog(hval2, 100hval2.^(-2), "--", label=L"O(h^{-2})")
  legend(loc="best")
  # ylabel("% of bad points")
  xlabel(L"h")
  title("Interior nodes")
  suptitle("Two conjoined spheres")
  savefig(newdir_nrun*"/surf_PB_badpoints_plot.png")
  savefig(newdir_nrun*"/surf_PB_badpoints_plot.eps")
end

# tests for molecular surface area
if false
  Xlim, Ylim, Zlim, Nx, h, dsignes = read_protein("2022"; scaling=0)
end
if false # tests with finest grid
  nrunf!(nrun)
  # reference values
  farray = [(x->1.0), (x->x)];
  carray = Array{Float64,1}(1.4:-0.1:0.1);
  # carray = [1.0;0.5;0.25];
  N = length(carray)
  surf_ref = zeros(N,2);
  for i=1:N
    println("\nRun $i/$N")
    surf_ref[i,:] = PB_gen_shape_molecule_area(Nx; Xlim=Xlim, Ylim=Ylim, Zlim=Zlim, Nx=Nx, h=h, dsignes=dsignes, fε = farray[1], cf = carray[i], plotting_surface=(i==2))
  end
end

if ~false # test 256^3, PB system
  farray = [(x->x)]
  carray = [4]
  Xlim, Ylim, Zlim, Nx, h, dsignes = read_protein("2022"; scaling=0)
  Gpol_IBIM, Gpol_corr, surf_vals = PB_gen_shape_molecule_system(Nx; Xlim=Xlim, Ylim=Ylim, Zlim=Zlim, Nx=Nx, h=h, dsignes=dsignes, fε = farray[1], cf = carray[1])
end


if false # tests LTC molecular surface. Computing area
  vecx = [2;1;0]
  # farray = [(x->1.0), (x->x)];
  farray = [(x->1.0)];
  # carray = [5;4;3;2;1];
  carray = Array{Float64,1}(6:-0.2:0.5);
  N = length(carray)
  surf_vals = zeros(2,N,3)
  for i=1:3
    nrunf!(nrun)
    Xlim, Ylim, Zlim, Nx, h, dsignes = read_protein("2022"; scaling=vecx[i])
    for j=1:N
      println("\nRun $j/$N")
      surf_vals[:,j,i] = PB_gen_shape_molecule_area(Nx; Xlim=Xlim, Ylim=Ylim, Zlim=Zlim, Nx=Nx, h=h, dsignes=dsignes, fε = farray[1], cf = carray[j], count=i)
    end
  end
end

if false # plotting results from LTC area
  hvec = [0.993576;0.496788;0.248394]
  nrun[1]+=1
  surf_vals = readdlm(path*"/2022/2022-07/2022-07-05/moleculePB_tight_surf_vals_16.dat")
  s = Int.(readdlm(path*"/2022/2022-07/2022-07-05/moleculePB_tight_surf_vals_size_16.dat"))
  surf_vals = reshape(surf_vals, (s[1],s[2],s[3]))
  hvec = readdlm(path*"/2022/2022-07/2022-07-05/moleculePB_tight_hvec_16.dat")[:,1]
  carray = Array{Float64,1}(6:-0.2:0.5);
  # M = size(surf_ref,1)
  N = size(surf_vals,2)
  # errvec = zeros(M,N,)
  figure(8); clf()
  for i=1:2
    subplot(1,2,i)
    ylim([1e-4;1e-0])
    loglog(hvec[1:2],2e-2hvec[1:2],"--k",label="h")
    loglog(hvec[1:2],2e-2hvec[1:2].^2,"--b",label="h^2")
    loglog(hvec[1:2],2e-2hvec[1:2].^3,"--r",label="h^3")
    loglog(hvec[1:2],2e-2hvec[1:2].^4,"--m",label="h^4")
    suptitle("LTC Smooth molecule. Surface area rel. diff. consec. values")
  end
  for i=1:N
    vec1 = abs.(diff(surf_vals[1,i,:]))./abs.(surf_vals[1,i,end-1:end])
    vec2 = abs.(diff(surf_vals[2,i,:]))./abs.(surf_vals[2,i,end-1:end])
    subplot(1,2,1)
    # loglog(hvec[1:2], vec1, "-o", label="vals $i, J=1")
    loglog(hvec[1:2], vec1, "-o", label="epsl $(carray[i]), J=1")
    subplot(1,2,2)
    loglog(hvec[1:2], vec2, "-x", label="epsl $(carray[i]), J∼h^2")
  end
  subplot(1,2,1)
  legend(loc="best")
  subplot(1,2,2)
  legend(loc="best")
  savefig(newdir*"/LTC_protein_fixed_epsl.png")
end

if false
  hvec = [0.993576;0.496788;0.248394]
  nrun[1]+=1
  x = surf_ref
  open(newdir*"/moleculePB_tight_surf_ref_$(nrun[1]).dat","w") do io
    writedlm(io,x)
  end
  x = surf_vals
  open(newdir*"/moleculePB_tight_surf_vals_$(nrun[1]).dat","w") do io
    writedlm(io,x)
  end
  x = size(surf_ref)
  open(newdir*"/moleculePB_tight_surf_ref_size_$(nrun[1]).dat","w") do io
    writedlm(io,x)
  end
  x = size(surf_vals)
  open(newdir*"/moleculePB_tight_surf_vals_size_$(nrun[1]).dat","w") do io
    writedlm(io,x)
  end
  x = hvec
  open(newdir*"/moleculePB_tight_hvec_$(nrun[1]).dat","w") do io
    writedlm(io,x)
  end
  M = size(surf_ref,1)
  N = size(surf_vals,2)
  # errvec = zeros(M,N,)
  figure(8); clf()
  figure(9); clf()
  for i=1:N
    for j=1:M
      vec = abs.(surf_vals[2,i,:].-surf_ref[j,2])/surf_ref[j,2]
      figure(8)
      subplot(4,4,j)
      loglog(hvec, vec, "-o", label="vals $i, ref $j")
      vec = abs.(surf_vals[1,i,:].-surf_ref[j,1])/surf_ref[j,1]
      ylim([1e-4;1e-0])
      figure(9)
      subplot(4,4,j)
      loglog(hvec, vec, "-o", label="vals $i, ref $j")
      ylim([1e-4;1e-0])
    end
  end
  for j=1:M
    figure(8)
    subplot(4,4,j)
    loglog(hvec,hvec,"--k")
    loglog(hvec,hvec.^2,"--b")
    loglog(hvec,hvec.^3,"--r")
    loglog(hvec,hvec.^4,"--m")
    legend(loc="best")
    suptitle("Second order Jac")
    figure(9)
    subplot(4,4,j)
    loglog(hvec,hvec,"--k")
    loglog(hvec,hvec.^2,"--b")
    loglog(hvec,hvec.^3,"--r")
    loglog(hvec,hvec.^4,"--m")
    legend(loc="best")
    suptitle("J=1")
  end
end

# solving system on torus
if false  # ε~2h test on the surfaces
  nrunf!(nrun)

  Ttarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_t.dat")[:,1]
  Ptarget = readdlm(path*"/2022/2022-03/2022-03-21/14/surface_tp_14_trg_p.dat")[:,1]
  
  TPtargets = [ [Ttarget[i];Ptarget[i];(0.5*Ttarget[i]-0.12)] for i=1:length(Ttarget) ]

  Ntargets = length(Ttarget)
  println("Number of targets = $(Ntargets)")

  x = TPtargets;
  open(newdir_nrun*"/torus_err_const_$(nrun[1])_trg_xyz.dat","w") do io
    writedlm(io,x)
  end

  # Nvec = 64:10:404
  # Nvec = 64:10:304
  # Nvec = 64:10:264
  # Nvec = 64:8:400
  # Nvec = 64:8:300
  # Nvec = [40;64;80]
  # Nvec = [23]
  Nvec = [20;40;80]
  # Hvec = 0.1*(0.5.^(0:3));
  Hvec = 0.25*(0.5.^(0:3));
  # Nvec = [25;50;100]
  # Nvec = 64:8:248
  # Nvec = 64:10:224
  
  # fε(x) = (0.30*x.^0.5)
  fε(x) = 0.9*x.^0.5
  # fε(x) = (2*x.^0.8)
  # fε(x) = 4x
  # fε(x) = 4x
  # detail = "epsl03h05"
  # detail = "epsl2h08"
  # detail = "epsl4h"
  detail = "torus_epsl5h"

  nmax = length(Nvec)
  errvec = zeros(nmax,2,3)
  hval = zeros(nmax)
  mval_abs = zeros(8,nmax)
  mval_err = zeros(4,nmax)
  surf_val= zeros(2,nmax)
  
  val_abs = zeros(8,Ntargets,nmax)
  val_err = zeros(4,Ntargets,nmax)

  shift=[0.026834;0.01346254;-0.06293856]

  @time for i=1:nmax
      println("\nRun $i/$nmax")
      local h = 0.1
      ε = 0.1
      # @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-0.75, Zlim2=0.75, zcenter=zeros(3) )

      @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Hvec[i], Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-1.75, Zlim2=1.75, zcenter=zeros(3), hflag=true )
      # @time hval[i], val_abs[:,:,i], mval_abs[:,i], surf_val[:,i], errvec[i,:,:] = PB_gen_shape_system( Nvec[i]; surfTargetXYZ=TPtargets, shift=shift, epslI=1.0, epslE=80.0, kappa_val=0.1257, fε=fε, plotting_surface=(i==0), count=i, Zlim1=-1.5, Zlim2=1.0, zcenter=zeros(3) )
      GC.gc()
  end
  x = detail;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_detail.dat","w") do io
    writedlm(io,x)
  end
  x = hval;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_hvals.dat","w") do io
    writedlm(io,x)
  end
  x = errvec;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_errvec.dat","w") do io
    writedlm(io,x)
  end
  x = val_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs.dat","w") do io
    writedlm(io,x)
  end
  x = val_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_abs);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_size.dat","w") do io
    writedlm(io,x)
  end
  x = size(val_err);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_size.dat","w") do io
    writedlm(io,x)
  end
  x = mval_abs;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_abs_means.dat","w") do io
    writedlm(io,x)
  end
  x = mval_err;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_val_err_means.dat","w") do io
    writedlm(io,x)
  end
  x = surf_val;
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val.dat","w") do io
    writedlm(io,x)
  end
  x = size(surf_val);
  open(newdir_nrun*"/surf_PB_$(nrun[1])_surf_val_size.dat","w") do io
    writedlm(io,x)
  end
end