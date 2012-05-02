clear all

dir = 'E:\Work\instability\ROMS\si_part\edge\2D\';
dirs = {'run01','run02','run03','run04','run05','run06','run07'};
%dirs = {'run02'};
runx = [01 02 03 04 04.2 05 06 07];
fname = 'ocean_his.nc';

plot_flag = 0;
redo_en = 0;
redo_pv = 0;

for ii=1:length(dirs)
    cd([dir dirs{ii}]);
    
    % read data    
    roms_grid = roms_get_grid(fname,fname,0,1);
    misc = roms_load_misc(fname);
    time = ncread(fname,'ocean_time');
    f0 = mean(misc.f(:));
    g = 9.81;

    dx = mean(mean(diff(roms_grid.x_rho,1,2)));
    dy = mean(mean(diff(roms_grid.y_rho,1,1)));
    dz = squeeze(diff(roms_grid.z_r(:,1,1)));

    xmid = ceil(size(roms_grid.x_rho,2)/2);
    ymid = ceil(size(roms_grid.y_rho,1)/2);
    zmid = ceil(size(roms_grid.z_r  ,1)/2);
    
    j = ymid;
    
    temp = squeeze(double(ncread(fname,'temp',[1 1 1 1],[Inf Inf Inf Inf])));
    v = squeeze(double(ncread(fname,'v',[1 1 1 1],[Inf Inf Inf Inf])));
    u = squeeze(double(ncread(fname,'u',[1 1 1 1],[Inf Inf Inf Inf])));
    temp_mid = squeeze(double(ncread(fname,'temp',[1 ymid 1 1],[Inf 1 Inf 1])));

    if ~exist('energy-avg-x.mat','file') || redo_en == 1, roms_energy(fname,[],4,1,0); end
    if ~exist('ocean_pv.nc','file') || redo_pv == 1, roms_pv(fname,[1 1]); end

    load energy-avg-x.mat

    % figure out width / wavelength of edge effects region
    % I AM SMOOTHING THE GROWTH RATE CURVE HERE
    [peaks,locs] = findpeaks(conv(A,[1 1]/2,'valid'));
    %tAmax = time_A(find(A == max(A)));
    tAmax = time_A(locs(1));
     tind = find_approx(time,tAmax,1);

    % NEED TO DO BETTER AVERAGING? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

    % first left, then right
    for i=1:2
          clear m;
          rr = range{i};
         dat = data(rr);
        xdat =   xu(rr);

        % find maxima & calculate width
           m = trapz(enwave(:,rr),1);
         ind = find(m > thresh*max(m));    
        w(i) = (ind(end)-ind(1))*dx;
        edgeloc = [(rr(1)+ind(1))*dx/1000 (rr(1)+ind(end))*dx/1000];

        % calculate wavelength
             m = enwave(:,rr);
        maxind = find(m == max(m(:)));
         [a,b] = ind2sub(size(m),maxind);
          l(i) = period(a);
          
        % find negative pv region
         pv = squeeze(ncread('ocean_pv.nc','pv',[1 ymid 1 1],[Inf 1 Inf 1]));
        xpv = ncread('ocean_pv.nc','x_pv');
        
        pvind = find(pv(:,zmid) < 0);
         npvl = xpv(pvind(1))/1000; % negative pv left
         npvr = xpv(pvind(end))/1000; % and right
        pvind = find(pv(:,zmid) ==  pv(xmid,zmid));
         cpvl = xpv(pvind(1))/1000; % constant pv left
         cpvr = xpv(pvind(end))/1000; % and right
         
        % find initial mean pv in the edge effects region.
        xpvl = find_approx(xpv,edgeloc(1)*1000,1);
        xpvr = find_approx(xpv,edgeloc(2)*1000,1);
        pvedge = pv(xpvl:xpvr,:);
        meanpv = max(pvedge(:));
        
        plotpv(ii,i) = meanpv;
        
        % find 
        
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
            linex([npvl npvr cpvl cpvr],' ','w');% mark negative + constant pv region
        end
    end

    if plot_flag
        figure(h1); title(['Width = ' num2str(w./1000) ' km']);
        figure(h2); title(['Width = ' num2str(w./1000) ' km']);
    end
    
    % save data
    wsi(ii,:) = w;
    lsi(ii,:) = l;
end

% plots
plotx = abs(plotpv);
figure
% subplot(311)
% plot(runx,lsi(:,1)/1000,'r*'); hold on
% plot(runx,lsi(:,2)/1000,'b*');
% legend('Left','Right');
% ylabel('Wavelength (km)');
% subplot(312)
% plot(runx,wsi(:,1)/1000,'r*'); hold on
% plot(runx,wsi(:,2)/1000,'b*'); hold on
% ylabel('Width (km)');
% subplot(313)
% plot(runx,wsi(:,1)./lsi(:,1),'r*'); hold on
% plot(runx,wsi(:,2)./lsi(:,2),'b*'); hold on
% ylabel('No. of wavelengths');
% ylim([1 3]);

subplot(311)
plot(plotx(:,1),lsi(:,1)/1000,'r*'); hold on
plot(plotx(:,2),lsi(:,2)/1000,'b*');
legend('Left','Right');
ylabel('Wavelength (km)');
subplot(312)
plot(plotx(:,1),wsi(:,1)/1000,'r*'); hold on
plot(plotx(:,2),wsi(:,2)/1000,'b*'); hold on
ylabel('Width (km)');
subplot(313)
plot(plotx(:,1),wsi(:,1)./lsi(:,1),'r*'); hold on
plot(plotx(:,2),wsi(:,2)./lsi(:,2),'b*'); hold on
ylabel('No. of wavelengths');
ylim([1 3]);
