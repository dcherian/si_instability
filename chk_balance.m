
[vars,atts,dims] = ncdfread('state.0000000000.t001.nc')
ph  = ncread('phiHyd.0000000000.t001.nc','phiHyd');
phl = ncread('phiHydLow.0000000000.t001.nc','phiHydLow');

%%

timeph = 1;
timefv = timeph;

nz = length(vars.Z);
ymid = ceil(length(vars.Y)/2);
dx = max(diff(vars.X,1,1));
v = vars.V;


dedx = diff(vars.Eta(1:end-1,ymid,timefv),1,1)./dx;
dpdx = diff(ph(1:end-1,:,:,timeph),1,1)./dx;
fv   = 1E-4*v(2:end-1,:,:,timefv);

dpfv = squeeze(dpdx(:,ymid,:,1) + fv(:,ymid,:));
dpfveta = squeeze(dpdx(:,ymid,:) + fv(:,ymid,:)) + 9.81*repmat(dedx(:,ymid),1,nz);
