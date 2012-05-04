%% choose run

fname  = 'ocean_his.nc';
pvname = 'ocean_pv.nc';
enname = 'energy-avg-x.mat';
lname  = 'length_scales_u.mat';

volume = {};

% First calculate diagnostics if not already done so
if ~exist(enname,'file'), roms_energy(fname,[],4,1,0); end % average along x
if ~exist(pvname,'file'), roms_pv(fname,[]); end
if ~exist(lname,'file'), roms_length_scales(fname,'u',[],volume,0); end

% grid and misc stuff
roms_grid = roms_get_grid(fname,fname,0,1);
misc = roms_load_misc(fname);
time = ncread('ocean_his.nc','ocean_time');
f0 = mean(misc.f(:));
g = 9.81;

dx = mean(mean(diff(roms_grid.x_rho,1,2)));
dy = mean(mean(diff(roms_grid.y_rho,1,1)));
dz = squeeze(diff(roms_grid.z_r(:,1,1)));

xmid = ceil(size(roms_grid.x_rho,2)/2);
ymid = ceil(size(roms_grid.y_rho,1)/2);
zmid = ceil(size(roms_grid.z_r  ,1)/2);

% load data
load energy-avg-x.mat
load length_scales_u.mat
%load('pv.mat','tpv');

temp = squeeze(double(ncread(fname,'temp',[1 1 1 1],[Inf Inf Inf Inf])));
v = squeeze(double(ncread(fname,'v',[1 1 1 1],[Inf Inf Inf Inf])));
u = squeeze(double(ncread(fname,'u',[1 1 1 1],[Inf Inf Inf Inf])));
temp_mid = squeeze(double(ncread(fname,'temp',[1 ymid 1 1],[Inf 1 Inf 1])));
    
%% Energy and Growth rates

% Symmetric Instability
% calculate growth rate from linear stability theory for initial state
yind = ymid;
% M2 N2
M2 = g*misc.Tcoef*avg1(squeeze(diff(temp(:,yind,:,:),1,1)),2)./dx;
N2 = g*misc.Tcoef*avg1(bsxfun(@rdivide,diff(squeeze(temp(:,yind,:,:)),1,2),permute(dz,[3 1 2])),1);

M2mid = squeeze(M2(xmid,zmid,1));
N2mid = squeeze(N2(xmid,zmid,1));

m = 2*pi/100;
sigma0 = sqrt(M2mid^2./N2mid - f0^2); %- m^2*misc.nl_visc2*((M2mid./N2mid)^2 + 1);

% compare against energy calculation of growth rate
load energy-avg-x
figure
plot(time_A./86400,A*86400);
liney(sigma0*86400)
legend('Growth rate from energy curve','2D most unstable growth rate','Location','Best');
ylabel('Growth Rate (d^{-1})');
xlabel('Time (days)');

% Baroclinic instability - use to validate general picture of what is happening.

%% PV fluxes at edge of initally unstable region

% find initial symmetrically unstable region
pvi = squeeze(double(ncread('ocean_pv.nc','pv',[1 1 16 1],[Inf 1 1 1])));
ind = find(pvi < 0);    
volume = {'x' ind(1) ind(end)};

% calculate fluxes
qfluxl = roms_pv_flux(fname,pvname,[],'x','5000'); % left
qfluxr = roms_pv_flux(fname,pvname,[],'x','45000'); % right

figure;
plot(tpv./86400,qfluxl,'b');
hold on
plot(tpv./86400,qfluxr,'r');
legend('Left','Right','Location','best');

%% find angle of isotherms + most unstable SI mode growth rate + calculate Ri
clear slope mean_slope sigma M2grad N2grad

yind = ymid;
s = size(temp);

sigma = sqrt(fillnan((M2.^2)./N2 - f0^2,Inf)) * 86400;
slope = fillnan(-M2./N2,Inf);

% Velocity gradients
vz = avg1(bsxfun(@rdivide,squeeze(diff(v(:,yind,:,:),1,3)),permute(dz,[3 1 2])),1);
uz = bsxfun(@rdivide,squeeze(diff(u(:,yind,:,:),1,3)),permute(dz,[3 1 2]));

Ri = addnan(abs(N2./(uz.^2 + vz.^2)),0.25);
%animate(Ri,'pcolor;shading flat');

for i=1:s(end)   
    M2i = M2(ind,:,i);
    N2i = N2(ind,:,i);
    
    % kludge to find angle of max gradient - why is this needed
    M2grad(i) = nanmedian(addnan(M2i(:),-1e-6));
    N2grad(i) = nanmedian(-addnan(-N2i(:),-1e-6));
    
    % calculate isotherm slopes
    mean_slope(i) = -M2grad(i)./N2grad(i); %atan(-N2grad(i)/M2grad(i))*180/pi;
end

plot(time./86400,mean_slope)
xlabel('Time (days)');
ylabel('Slope');
title('(Mean) Rotation of isotherms');

% OLD
%     mean_slope(i) = atand((-mean(M2(:))./mean(N2(:))));
%     
%     slope(:,:,i) =  atand((-M2(2:end,2:end-1)./N2(2:end-1,2:end)));

%% FFT
clear freq cn
for i=2:size(u,4)
    [freq(:,i), ~,cn(:,i),~] = dcpgram(squeeze(u(100,:,27,i)),312.5,0);
end

pgram = abs(cn).^2;
x = [min(freq(:)) max(freq(:))];
y = [min(pgram(:)) max(pgram(:))];

for i=2:size(u,4)
    loglog(freq(:,i),pgram(:,i),'LineWidth',1.5);
    xlim(x);
    ylim(y);
    grid on;
    pause;
end
    
%% Find width of region with instabilities

% Find region where symmetric rolls exist

% problems 
%   - gradient in PV means right side has lesser growth rate.
%   - should I calculate growth rates for each region separately and then calculate width?

% calculate theoretical wavelengths and growth rates

clear all

%dir = 'E:\Work\instability\ROMS\si_part\edge\2D\';
%dirs = {'run01','run02','run03','run04','run04_2','run05','run06','run07'};
%dirs = {'run02'};
%plotx = [01 02 03 04 04.2 05 06 07];
fname = 'ocean_his.nc';

plot_flag = 1;
redo_en = 0;

%for ii=1:length(dirs)
%    cd([dir dirs{ii}]);
    u = double(ncread(fname,'u',[1 1 1 1],[Inf Inf Inf Inf]));

    roms_grid = roms_get_grid(fname,fname,0,1);
    time = ncread(fname,'ocean_time');
    j = 3;

    if ~exist('energy-avg-x.mat','file') || redo_en == 1, roms_energy(fname,[],{},4,1,0); end

    load energy-avg-x.mat

    % I AM SMOOTHING THE GROWTH RATE CURVE HERE
    [peaks,locs] = findpeaks(conv(A,[1 1]/2,'valid'));
    %tAmax = time_A(find(A == max(A)));
    tAmax = time_A(locs(1));
     tind = find_approx(time,tAmax,1);

    % need to do better averaging
    um   = mean(u(:,j,:,:),1);
    up = squeeze(bsxfun(@minus,u(:,j,:,:),um));
    data = up(:,20,tind);
    thresh = 0.5;

    xu = roms_grid.x_u(1,:);
    dx = mean(diff(xu));
    nx = length(xu);
    range = {[1:floor(nx/2)], [floor(nx/2)+1:nx]};

    [wave,period,scale,coi] = wavelet(data,dx,1);
    enwave = abs(wave).^2;

    % plot energy 'spectrum'
    if plot_flag
        h1 = figure;
        pcolor(xu/1000,log(period),enwave); shading interp
        hold on; revz
        plot(xu/1000,log(coi),'k','LineWidth',1.5);
        xlabel('X (km)'); ylabel('log(Wavenumber)');

        % plot u field
        mod_movie(fname,'u',[tind tind],{},'y','mid','pcolor;shading interp;liney(-50)');
        hold on;
        h2 = gcf;
    end

    for i=1:2
          clear m;
          rr = range{i};
         dat = data(rr);
        xdat =   xu(rr);

        % find maxima & calculate width
           m = trapz(enwave(:,rr),1);
         ind = find(m > thresh*max(m));    
        w(i) = (ind(end)-ind(1))*dx;

        % calculate wavelength
             m = enwave(:,rr);
        maxind = find(m == max(m(:)));
         [a,b] = ind2sub(size(m),maxind);
          l(i) = period(a);

        % put on 
        if plot_flag
            figure(h1)
            linex((rr(1)+ind(1))*dx/1000,''); linex((rr(1)+ind(end))*dx/1000,'');
            liney(log(period(a)));

            figure(h2)
            linex((rr(1)+ind(1))*dx/1000,''); linex((rr(1)+ind(end))*dx/1000,'');
            liney(log(period(a)));
        end
    end

    if plot_flag
        figure(h1); title(['Width = ' num2str(w./1000) ' km']);
        figure(h2); title(['Width = ' num2str(w./1000) ' km']);
    end
    w
    l
    % save data
 %   wsi(ii,:) = w;
%    lsi(ii,:) = l;
%end

%% compare with 2D run
dir = 'E:\Work\instability\ROMS\si_part\edge\';
dir2D = [dir 'run150_lsi_2d1.5\'];
dir3D = [dir 'run150_lsi_redo2\'];

file2d = [dir2D 'ocean_his.nc'];
file3d = [dir3D 'ocean_his.nc'];

% add some slab stuff
% %% movies
% u2d = ncread(file2D,'u');
% u3d = ncread(file3D,'u');
% 
% diffu = bsxfun(@minus,u3d,u2d);
% animate(diffu);

% compare diagnostics

h1 = figure;

volume = {'x' '12000' '32000'};

% 2d
cd(dir2D);
if ~exist('length_scales_u.mat')
    roms_length_scales(file2D,'u',[],volume);
end
if ~exist('energy-avg-x.mat')
    roms_energy(file3D,[],4,1,0);
end

load length_scales_u.mat
L2d = L;
t2d = time_L;figure(h1)
subplot(212)
plot(t2d/86400,L2d,'r','LineWidth',1.5); hold on
ylabel('Length Scales (m)');

load energy-avg-x.mat
A2d = A;
tA2d = time_A;
figure(h1)
subplot(211)
plot(tA2d/86400,A2d*86400,'r','LineWidth',1.5); hold on

% 3d
cd(dir3D);
if ~exist('length_scales_u.mat')
    roms_length_scales(file3d,'u',[],volume);
end
if ~exist('energy-avg-x.mat')
    roms_energy(file3d,[],4,1,0);
end

load length_scales_u.mat
L3d = L;
t3d = time_L;
figure(h1);
subplot(212)
plot(t3d/86400,L3d,'k','LineWidth',1.5); hold on
xlabel('Time (days)');

load energy-avg-x.mat
A3d = A;
tA3d = time_A;
figure(h1)
subplot(211)
plot(tA3d/86400,A3d*86400,'k','LineWidth',1.5); hold on
legend('2D','3D','Location','Best');
ylabel('Growth Rate (d^{-1})');

%% old scaling
N2 = 1e-5;
vz = 4e-3;
v0 = 0.2;
f  = 1e-4;
nu = 1e-4; % 3e-2

M2 = f*vz;
Ri = N2/vz^2;

sigma = f*(1/Ri - 1)^0.5 * 86400
sigma0 = sqrt(M2^2./N2 - f^2) * 86400
L_stone = 2*pi / (pi/sqrt(1-Ri) * f/v0)
%L_raf = 2*pi/sqrt(f/nu*sqrt(1/Ri - 1)/(1 + N2^2/M2^2))

% observed - 2.5km in x - works for mid level vleocity v0 = 0.2 m/s

%% Really Old Stuff

% length scales
% % 
% % % interpolate to standard z grid and find lz
% % dx = 208.33;
% % dy = 375;
% % u = ncread('ocean_his.nc','u');
% % for i=1:size(u,4)
% %     lx(i) = length_scale(u(:,:,15:32,i),1,dx);
% %     ly(i) = length_scale(u(:,:,15:32,i),2,dy);
% %     %lz(i) = length_scale(u(:,:,:,i),3,dz);
% % end
% % 
% % figure;
% % hold on
% % plot(lx,'r');
% % plot(ly,'g');
% % %plot(lz,'b');
% % legend('lx','ly','lz');
% 
% %%
% fname = 'ocean_his.nc';
% grid = roms_get_grid(fname,fname,0,1);
% 
% %% Plot v movie
% roms_diffmovie(fname,'v',[1 Inf],'y',3);
% 
% %% rho movie
% roms_diffmovie(fname,'u',[1 Inf],'y',3);
% 
% %% PV movie
% pv = roms_pv('ocean_his.nc',[1 Inf]);
% animate(pv(:,3,:,:));
% %%
% pv1 = pv > 0;
% animate(pv1(:,3,:,:));
% 
% %% integrate PV
% dz = permute(diff(grid.z_u,1,1),[3 2 1]);
% int_pv1 = (pv.*repmat(dz(:,1:end-1,:),[1 1 1 size(pv,4)]));
% int_pv = zeros(size(pv,4),1);
% for i=1:size(pv,4)
%     temp = int_pv1(:,:,:,i);
%     int_pv(i) = sum(temp(:));
% end
% 
% figure;
% plot(int_pv);
% ylabel('PV');
% xlabel('time');
% 
% %% Integrate energy
% 
% u = ncread(fname,'u');
% v = ncread(fname,'v');
% rho = ncread(fname,'rho');
% 
% KE1 = 0.5*rho(1:end-1,1:end-1,:,:).*(u(:,1:end-1,:,:).*u(:,1:end-1,:,:) + v(1:end-1,:,:,:).*v(1:end-1,:,:,:));
% PE1 = 9.81*rho.*repmat(permute(grid.z_r,[3 2 1]),[1 1 1 size(rho,4)]);
% 
% int_KE = zeros(size(KE1,4),1);
% int_PE = zeros(size(PE1,4),1);
% for i=1:size(pv,4)
%     temp = KE1(:,:,:,i);
%     int_KE(i) = sum(temp(:));
%     temp = PE1(:,:,:,i);
%     int_PE(i) = sum(temp(:));
% end
% 
% figure;
% plot(int_PE,'b'); hold on;
% plot(int_KE,'r');
% %plot(int_PE + int_KE,'k');
% ylabel('Energy');
% xlabel('time');
% legend('PE','KE','TE');
% 
% %% OLD CODE
% % [vars,atts,dims] = ncdfread(fname);
% % grid = roms_get_grid(fname,fname,0,1);
% % 
% % %% Plot v movie
% % 
% % yindex = 10;
% % s = size(vars.v);
% % dv = zeros([s(1) s(3) s(4)]);
% % for i=1:s(end)
% %     dv(:,:,i) = squeeze(vars.v(:,yindex,:,i)-vars.v(:,yindex,:,1));
% % end
% % 
% % animate(grid.x_v(1,:),grid.z_v(:,1,1),dv);
% % 
% % %% Plot temp movie
% % 
% % yindex = 10;
% % s = size(vars.temp);
% % dv = zeros([s(1) s(3) s(4)]);
% % for i=1:s(end)
% %     dv(:,:,i) = squeeze(vars.temp(:,yindex,:,i)-vars.temp(:,yindex,:,1));
% % end
% % 
% % animate(grid.x_rho(1,:),grid.z_r(:,1,1),vars.rho(:,yindex,:,:));