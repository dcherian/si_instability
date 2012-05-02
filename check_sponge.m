% run in directory iwth ocean_his.nc file
fname = 'ocean_his.nc';
roms_grid = roms_get_grid(fname,fname,0,1);

xr = roms_grid.x_rho;


%%

USER = [10000 10];

fac = USER(2);

visc2 = 10;
ng = 1;
jstart = 1; istart = 1;
jend = size(xr,2); iend = size(xr,1);

dx = xr(1,2)-xr(1,1);

for j=jstart:jend
   for i=istart:iend
           if (xr(i,j) < USER(1))
             cff=(USER(1)-xr(i,j))*fac*visc2(ng)/USER(1);
             visc2_r(i,j)=cff;
             visc2_p(i,j)=cff;
           else
             if (xr(i,j) > (xr(1,jend)-USER(1)))
                 cff=(xr(i,j)+USER(1)-xr(1,jend))*fac*visc2(ng)/USER(1);
                 visc2_r(i,j)=cff;
                 visc2_p(i,j)=cff;
             else
                 visc2_r(i,j)=0.0;
                 visc2_p(i,j)=0.0;3
             end
           end
   end
end

plot(visc2_r(5,:))
