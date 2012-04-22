fname = 'ocean_his.nc';

vinfo = ncinfo(fname,'temp');
s = vinfo.Size;
mids = ceil(s/2);
f0 = 1e-4;
g = 9.81;
rho0 = ncread('ocean_his.nc','rho0');
TCOEF = ncread('ocean_his.nc','Tcoef');

temp_field = squeeze(double(ncread('ocean_his.nc','temp',[1 1 1 1], [Inf 1 Inf Inf])));
%rho_field = squeeze(double(ncread('ocean_his.nc','rho',[1 1 1 1], [Inf 1 Inf Inf])));
grid = roms_get_grid('ocean_his.nc','ocean_his.nc',0,1);
time = ncread('ocean_his.nc','ocean_time')./86400;
temp = double(squeeze(ncread(fname,'temp',[1 mids(2) mids(3) 1],[Inf 1 1 Inf])));

dist = (bsxfun(@minus,temp,temp(:,1)));

xmat = ncread('ocean_his.nc','x_rho');
dx = (diff(xmat,1));
dx =  mean(dx(:));

grad = diff(temp,1,1)./dx;
grad = grad(2:end-1,:);

%%
figure;

for i=1:size(grad,2)
    subplot(2,1,1)
    plot(dist(i,:))
    ylim([min(dist(:)) max(dist(:))]);
    
    subplot(2,1,2)
    plot(grad(i,:));
    ylim([min(grad(:)) max(grad(:))]);
    pause(0.3);
end

%%
grad1 = bsxfun(@minus,grad,median(grad));
i = 80;

ipeak = zeros(s(end),5);

for i=1:s(end)
    clf
    [peaks,ind] = findpeaks(smooth(grad(:,i).^2,10),'npeaks',5,'sort','descend');
    plot(grad(:,i));
    hold on
    if length(ind) >=3, linex(ind(1:3)); end
    title(num2str(i))
    pause(0.1)
    
    ipeak(i,1:length(ind)) = xmat(ind,1);
end


%% Find region where symmetric rolls exist

u = double(ncread('ocean_his.nc','u',[1 1 1 1],[Inf Inf Inf Inf]));

j = 50;

% need to do better averaging
um   = mean(u(:,j,:,:),1);
up = squeeze(bsxfun(@minus,u(:,j,:,:),um));

%%
pcolor(avg1(xmat(:,1))/1000,time,squeeze(up(:,15,:))'); shading interp
linex(mean(avg1(xmat(:,1))/1000));

%%
ind = 20;
contourf(grid.x_rho(1,:)',grid.z_r(:,1,1),temp_field(:,:,i)',20);
linex(ind);
shading flat
hold on
contour(grid.x_rho(1,:)',grid.z_r(:,1,1),temp_field(:,:,1)',10,'w')

%% Compare PV and crappy region

[pv,xpv,ypv,zpv] = roms_pv('ocean_his.nc',[1 1]);
temp = ncread('ocean_his.nc','temp',[1 1 1 5],[Inf Inf Inf 1]);
a = char('contour(xpv,zpv,squeeze(pv(:,1,:,1))'',''w'')');
mod_movie('ocean_his.nc','u',[5 5],{},'y','mid','pcolor; shading interp;');
hold on
contour(xpv./1000,zpv,squeeze(pv(:,1,:,1))',40,'w');

