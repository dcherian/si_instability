% Things to do
%   - calculate growth rates for each half separately
%   - better averaging for mean
%   - calculate width at a bunch of depths and average
%   - negative N2!!!

% growth rate corresponding to particular wavelength

clear all

dir = 'E:\Work\instability\ROMS\si_part\edge\2D\';
%dirs = {'run01','run02','run03','run04','run05','run06','run07'};
dirs = {'run01','run02','run03','run04','run05','run06','run07','run08','run09','run10_2','run10_3','run10_4','run10_5','run11','run12_2','run13_2','run14_2','run15', ...
        'run16','run17','run18','run19','run21','run23','run24'}; % 4 and 7 are outliers (PVmin / PVmid > 1.1 - greater wavelength too)
runx = [01 02 03 04 05 06 07 08 09 10.2 10.3 10.4 10.5 11 12.2 13.2 14.2 15 16 17 18 19 21 23 24];
fname = 'ocean_his.nc';
volume = {};

ennames = {'energy-left.mat','energy-right.mat'};

plot_flag = 0;
redo_en = 0;
redo_pv = 0;

%dirs = {'run10_2','run10_3','run10_4','run10_5'}; runx = [10.2 10.3 10.4 10.5];
dir = 'E:\Work\instability\ROMS\si_part\edge\3D\'; dirs = {'run01-2D','run01','run01-bfric-1','run01-bfric-2','run01-bfric-3'}; runx=[0 0.5 1 2 3];
%dirs = {'run24'}; runx = [24];plot_flag = 1;

thresh = 0.5; % threshold for energy

for ii=1:length(dirs)
    cd([dir dirs{ii}]);
    
    % get grid details    
    roms_grid = roms_get_grid(fname,fname,0,1);
    misc = roms_load_misc(fname);
    time = ncread(fname,'ocean_time');
    f0 = mean(misc.f(:));
    g = 9.81;
    bfric(ii) = misc.rdrg;

    dx = mean(mean(diff(roms_grid.x_rho,1,2)));
    dy = mean(mean(diff(roms_grid.y_rho,1,1)));
    dz = squeeze(diff(roms_grid.z_r(:,1,1)));
    
    xu = roms_grid.x_u(1,:);
    xv = roms_grid.x_v(1,:);
    nx = length(xu);

    xmid = ceil(size(roms_grid.x_rho,2)/2);
    ymid = ceil(size(roms_grid.y_rho,1)/2);
    zmid = ceil(size(roms_grid.z_r  ,1)/2);

    yind = ymid;
    
    % read data
    temp = squeeze(double(ncread(fname,'temp',[1 ymid 1 1],[Inf 1 Inf Inf])));
       v = squeeze(double(ncread(fname,'v',[1 ymid 1 1],[Inf 1 Inf Inf])));
       u = squeeze(double(ncread(fname,'u',[1 ymid 1 1],[Inf 1 Inf Inf])));
    %temp_mid = squeeze(double(ncread(fname,'temp',[1 ymid 1 1],[Inf 1 Inf 1])));
    
    if ~exist('ocean_pv.nc','file') || redo_pv == 1, roms_pv(fname,[1 1]); end   
    pv = squeeze(ncread('ocean_pv.nc','pv',[1 ymid 1 1],[Inf 1 Inf 1]));
    xpv = ncread('ocean_pv.nc','x_pv');
    
    %%%%% find different pv & v regions
    [npvl cpvl cpvr npvr] = find_region(xpv,pv);
    [nvl  cvl   cvr  nvr] = find_region(xv,squeeze(v(:,:,1)));
     % x-length scale for velocity gradient (Rossby number calculations)
     Lvx = [cvl-nvl nvr-cvr];
     range_vx = {[find(xv == nvl):find(xv == cvl)], [find(xv == cvr):find(xv == nvr)]};
     range_nx = {[1:floor(nx/2)], [floor(nx/2)+1:nx]};
    
    %%%%%%%%%%%%%%%%%%% first left, then right
    for i=1:2
        rr = range_vx{i};
        volume = {'x' rr(1) rr(end)};
        % calculate / load energy diagnostics for current region rr
        if ~exist(ennames{i},'file') || redo_en == 1
            roms_energy(fname,[],volume,1,1,'growthrate_u'); 
            system(['move energy-avg-x.mat ' ennames{i}]);
        end
        load(ennames{i});
        
        %%%%%%%%%%%%%%%% find time of max growth rate
        % I AM SMOOTHING THE GROWTH RATE CURVE HERE
        %[peaks,locs] = findpeaks(conv(A./max(A),[1 1]/2,'valid'),'threshold',0.01);
        [peaks,locs] = extrema(conv(A./max(A),[1 1]/2,'valid'));
        Amax(ii,i) = peaks(1);
        %if locs(2) < locs(1) && peaks(2)/peaks(1) >= 0.8
        %    tAmax = time_A(locs(2));
        %else
            tAmax(ii,i) = time_A(locs(1));
        %end
        tind = find_approx(time,tAmax(ii,i),1);
        % verify
        if plot_flag
            figure;
            plot(time_A./86400, A*86400);
            linex(tAmax(ii,i)./86400);
            title(['Growth Rate' ennames{i}]);
        end

        %%%%%%%%%%%%%% Wavelet stuff
        % NEED TO DO BETTER AVERAGING? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % do this calculation at a bunch of depths and average????? <-- better estimate of the width
        um   = mean(u(:,:,:),1);
        up = squeeze(bsxfun(@minus,u(:,:,:),um));
        data = up(:,zmid,tind);

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

        %%%%% figure out width / wavelength of edge effects region
          clear m;
          rr = range_nx{i}; % search for width in appropriate half of the domain
         dat = data(rr);
        xdat =   xu(rr);

        %%%%% find maxima & calculate width
           m = trapz(enwave(:,rr),1);
         ind = find(m > thresh*max(m));    
        w(i) = (ind(end)-ind(1))*dx;
        edgeloc = [(rr(1)+ind(1))*dx/1000 (rr(1)+ind(end))*dx/1000];

        %%%%% calculate x-wavelength
             m = enwave(:,rr);
        maxind = find(m == max(m(:)));
         [a,b] = ind2sub(size(m),maxind);
          lx(i) = period(a);
          
        %%%%% calculate z-wavelength/wavenumber
        % first interpolate to regular grid
        clear m
        data = up(ceil((ind(1)+ind(end))/2),:,tind);
        deltaz = 1; % m
        zin = [roms_grid.z_u(1,1,1):deltaz:roms_grid.z_u(end,1,1)]; % interpolate to this grid.
        datain = interp1(roms_grid.z_u(:,1,1),data,zin);
        lz(ii,i) = length_scale(datain,2,deltaz)*4; % in m
        
        %%%%% find initial mean pv in the edge effects region.
        xpvrr = xpv(rr(1):rr(end)-1); % this is now the current HALF of the domain
        xpvle = find_approx(xpvrr,edgeloc(1)*1000,1);
        xpvre = find_approx(xpvrr,edgeloc(2)*1000,1);
        pvedge = pv(xpvle:xpvre,:);
        meanpv = mean(pvedge(:)); %max(pvedge(:)./pv(xmid,zmid));
        
        plotpv(ii,i) = meanpv;
        
        %%%%% Calculate INITIAL gradients in the current HALF
        tindex = 1; 
        M2 = g*misc.Tcoef*avg1(                diff(squeeze(temp(rr(1):rr(end),:,tindex)),1,1),2)./dx;
        N2 = g*misc.Tcoef*avg1(bsxfun(@rdivide,diff(squeeze(temp(rr(1):rr(end),:,tindex)),1,2),permute(dz,[3 1 2])),1);
        
        %%%% find mean slope in EDGE effects region ||||||  THIS IS NOISY.
         slope = M2(xpvle:xpvre,:)./N2(xpvle:xpvre,:);
         plotsl(ii,i) = median(slope(:));
        
        %%%%% find Ri in EDGE effects region
         Ri = N2(xpvle:xpvre,:)./(M2(xpvle:xpvre,:).^2) * f0^2;
        Ri0 = median(Ri(Ri>0));
        plotri(ii,i) = Ri0;
        
        %%%%% Calculate Rossby number
        v0 = max(abs(v(:,:,:,1)));
        v0 = max(v0(:));
        plotro(ii,i) = v0/f0/Lvx(i);
        
        %%%%% calculate theoretical wavelength
        % first velocity scale
        xvl = find_approx(roms_grid.x_v(1,:),edgeloc(1)*1000,1);
        xvr = find_approx(roms_grid.x_v(1,:),edgeloc(2)*1000,1);
        %v0 = abs(mean(v(xvl:xvr,ymid,zmid,tind)));
        l0(i) = 2*pi / (pi/sqrt(1-Ri0) * f0/v0);
        
        %%%%% Theoretical and inferred growth rate
        Amax(ii,i) = peaks(1);
          A0(ii,i) = sqrt(1/Ri0-1)*f0*86400;
        
        %%%%% decorrelation length scale
        ldc(ii,i) = length_scale(v(xvl:xvr,:,tind),1,dx)*4;
        
%         % check pv 
%         figure
%         plot(xpv/1000,pv(:,zmid));
%         linex([npvl npvr cpvl cpvr],' ');
%         linex([(rr(1)+ind(1))*dx/1000 (rr(1)+ind(end))*dx/1000],' ','r');
%         pause

        % put on 
        if plot_flag
            figure(h1) % on energy plot
            linex((rr(1)+ind(1))*dx/1000,''); linex((rr(1)+ind(end))*dx/1000,'');
            liney(log(period(a)));

            figure(h2) % on u plot
            linex(edgeloc,''); % mark width
            linex([npvl npvr cpvl cpvr]/1000,' ','w');% mark negative + constant pv region
            
            figure(h1); title(['Width = ' num2str(w./1000) ' km']);
            figure(h2); title(['Width = ' num2str(w./1000) ' km']);
            if i==1, pause; end
        end
    end
        
    % save data
     wsi(ii,:) = w;
     lsi(ii,:) = lx;
    lsi0(ii,:) = l0; % theoretical
end

%% plots
plotx = (abs(plotri));
titlex = 'Ri';

figure
subplot(311)
plot(runx,lsi(:,1)/1000,'r*'); hold on
plot(runx,lsi(:,2)/1000,'b*');
plot(runx,lsi0(:,1)/1000,'ro');
plot(runx,lsi0(:,2)/1000,'bo');
%plot(runx,ldc(:,1)/1000,'rx');
%plot(runx,ldc(:,2)/1000,'bx');
legend('Left','Right','Location','Best');
ylabel('Wavelength (km)');
title('Runs in order');
subplot(312)
plot(runx,wsi(:,1)/1000,'r*'); hold on
plot(runx,wsi(:,2)/1000,'b*'); hold on
ylabel('Width (km)');
subplot(313)
plot(runx,wsi(:,1)./lsi(:,1),'r*'); hold on
plot(runx,wsi(:,2)./lsi(:,2),'b*'); hold on
ylabel('No. of wavelengths');
%ylim([1 3]);

figure
subplot(311)
plot(plotx(:,1),lsi(:,1)/1000,'r*'); hold on
plot(plotx(:,2),lsi(:,2)/1000,'b*');
plot(plotx(:,1),lsi0(:,1)/1000,'ro');
plot(plotx(:,2),lsi0(:,2)/1000,'bo');
%plot(plotx(:,1),ldc(:,1)/1000,'rx');
%plot(plotx(:,2),ldc(:,2)/1000,'bx');
legend('Left','Right','Location','Best');
ylabel('Wavelength (km)');
title('Richardson Number');
subplot(312)
plot(plotx(:,1),wsi(:,1)/1000,'r*'); hold on
plot(plotx(:,2),wsi(:,2)/1000,'b*'); hold on
ylabel('Width (km)');
subplot(313)
plot(plotx(:,1),wsi(:,1)./lsi(:,1),'r*'); hold on
plot(plotx(:,2),wsi(:,2)./lsi(:,2),'b*'); hold on
ylabel('No. of wavelengths');
xlabel(titlex);
%ylim([1 3]);

plotx = (abs(plotro));
titlex = 'Ro';

figure
subplot(311)
plot(plotx(:,1),lsi(:,1)/1000,'r*'); hold on
plot(plotx(:,2),lsi(:,2)/1000,'b*');
plot(plotx(:,1),lsi0(:,1)/1000,'ro');
plot(plotx(:,2),lsi0(:,2)/1000,'bo');
%plot(plotx(:,1),ldc(:,1)/1000,'rx');
%plot(plotx(:,2),ldc(:,2)/1000,'bx');
legend('Left','Right','Location','Best');
ylabel('Wavelength (km)');
title('Rossby Number');
subplot(312)
plot(plotx(:,1),wsi(:,1)/1000,'r*'); hold on
plot(plotx(:,2),wsi(:,2)/1000,'b*'); hold on
ylabel('Width (km)');
subplot(313)
plot(plotx(:,1),wsi(:,1)./lsi(:,1),'r*'); hold on
plot(plotx(:,2),wsi(:,2)./lsi(:,2),'b*'); hold on
ylabel('No. of wavelengths');
xlabel(titlex);

figure;
subplot(211)
plot(bfric,Amax(:,1),'r*');
plot(bfric,Amax(:,2),'b*');
ylabel('Max growth rate A_{max} d^{-1}');
subplot(212)
plot(bfric,tAmax(:,1)/86400,'r*');
plot(bfric,tAmax(:,2)/86400,'b*');
ylabel('Time of A_{max} (days)');
xlabel('Bottom Friction');

% print table
disp('     lsi_l     lsi_r     wsi_l     wsi_r     ri_l      ri_r      ro_l       ro_r      run');
disp( [lsi wsi plotri*1e4 plotro*1e4 runx'*1e4]/1e4)

%% 3d

F = TriScatteredInterp(plotro(:,1),plotri(:,1),wsi(:,1));
[rol,ril] = meshgrid(plotro(:,1),plotri(:,1));
winterp = F(rol,ril);

figure
plot3(plotro(:,1),plotri(:,1),wsi(:,1)./lsi(:,1),'r*'); hold on
plot3(plotro(:,2),plotri(:,1),wsi(:,2)./lsi(:,2),  'b*');
xlabel('Ro'); ylabel('Ri'),zlabel('wsi');
