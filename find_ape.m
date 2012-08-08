fname = 'E:\Work\instability\ROMS\si_part\edge\3D\run-ape\ocean_his.nc';

 rho = ncread(fname,'rho',[1 1 1 1],[Inf Inf Inf 1]);

rgrid = roms_get_grid(fname,fname,0,1);

xpsi = rgrid.x_psi(1,:);
  zw = rgrid.z_w(:,1,1);
  dy = mean(diff(rgrid.y_psi(:,1)));

ape = calc_ape(xpsi,zw,dy,rho(:,1,:))