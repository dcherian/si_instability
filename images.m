%% set variables
imdir = 'E:\Work\instability\graphs\talk\';

%% Ken's run cross section - u
dir = 'E:\Work\instability\ken102\';
file = [dir 'STABB102_his_0001.nc'];
mod_movie(file,'u',[10 10],{},'y','mid','pcolor;shading interp;nocaxis;fancy_cmap');

%% Ken's run cross section - pv
dir = 'E:\Work\instability\ken102\';
file = [dir 'ocean_pv.nc'];
file2 = [dir 'STABB102_his_0001.nc'];
handles = mod_movie(file2,'temp',[1 1],{},'y','mid','shading flat; fancy_cmap; topresent');

v = ncread(file2,'v',[1 1 1 1], [Inf Inf Inf 1]);
roms_grid = roms_get_grid(file2,file2,0,1);
hold on
[c2,h2] = contour(roms_grid.x_v(1,:)/1000,roms_grid.z_v(:,1,1),squeeze(v(:,ceil(size(v,2)/2),:))','b','LineWidth',1.5);
clabel(c2,h2,[-0.35:0.1:0],'FontSize',16,'Color','w','LabelSpacing',72*6,'Rotation',0);

pv = ncread(file,'pv',[1 1 1 1],[Inf Inf Inf 1]);
xpv = ncread(file,'x_pv');
zpv = ncread(file,'z_pv');
pvmid = squeeze(pv(:,size(pv,2)/2,:,1));
pvpos = pvmid; pvpos(pvmid<0)=0;
pvneg = pvmid; pvneg(pvmid>0) = 0;
hold on
[c,h] = contour(xpv/1000,zpv,pvneg',10,'k','LineWidth',1.5);
set(h,'LineWidth',1.5);

set(gca,'YTick',[-100:25:0]);
title('Temperature (shaded)','FontSize',24);
set(handles.h_cbar,'YTick',[10:1:18]);
set(get(handles.h_cbar,'title'),'String', 'Celsius');
set(get(handles.h_cbar,'Title'),'FontSize',20)

export_fig([imdir 'ken-cross.png']);

%% Initial conditions cross section


 dir = ['E:\Work\instability\ROMS\si_part\edge\3D\run01\'];
file = [dir 'ocean_pv.nc'];
file2 = [dir 'ocean_his.nc'];
clear handles

tindex = 1;

temp = ncread(file2,'temp',[1 240 1 tindex], [Inf 1 Inf 1]);
roms_grid = roms_get_grid(file2,file2,0,1);
pv = ncread(file,'pv',[1 1 1 1],[Inf Inf Inf 1]);
xpv = ncread(file,'x_pv');
zpv = ncread(file,'z_pv');
xmid = ceil(size(pv,1)/2);

pvmid = squeeze(pv(:,size(pv,2)/2,ceil(size(pv,3)/2),1));
pv0l = xpv(find_approx(pvmid(1:xmid),0))/1000;
pv0r = xpv(find_approx(pvmid(xmid:end),0) + xmid-1)/1000;

pvnl = 20.88;%xpv(find_approx(pvmid(1:xmid),pvmid(xmid),1))/1000;
pvnr = 49.3;%xpv(find(pvmid(xmid:end) == pvmid(xmid),1,'last') + xmid-1)/1000;
pvpos = pvmid; pvpos(pvmid<0)=0;
pvneg = pvmid; pvneg(pvmid>0) = 0;

% plot u
hfig2 = figure;
set(hfig2,'Units','pixels');
set(hfig2,'Position',[0 0 1600 900]);
handles.sub1 = subplot(211);
[c2,handles.contourf] = contourf(roms_grid.x_rho(1,:)/1000,roms_grid.z_r(:,1,1),squeeze(temp)',30); shading flat;
fancy_map = flipud(cbrewer('div', 'RdYlGn', 32));
colormap(fancy_map);
handles.h_cbar = colorbar

handles.lines = linex([8 70-8],'','w');
set(handles.lines,'LineStyle','-')
set(handles.lines,'LineWidth',2)

title('Temperature','FontSize',24);
set(get(handles.h_cbar,'title'),'String', 'Celsius');
set(get(handles.h_cbar,'Title'),'FontSize',20)
set(handles.h_cbar,'YTick',[27:3:42]);

set(gca,'YTick',[-50:10:0]);
ylabel('Z (m)', 'FontSize',20);
xtick = get(handles.sub1,'XTick');
xtick(3) = str2num(num2str(pv0l,'%.2f'));
xtick(6) = str2num(num2str(pv0r,'%.2f')); 
xtick(end+1) = str2num(num2str(pvnl,'%.2f'));
xtick(end+1) = str2num(num2str(pvnr,'%.2f'));
xtick = sort(xtick);
set(handles.sub1,'XTick',xtick);
    
% Now do PV , V

figure(hfig2)
handles.sub2 = subplot(212)
hold on

[x,y,z,~,~,~] = roms_var_grid(file2,'v');
v = double(squeeze(ncread(file2,'v',[1 ceil(length(y)/2) length(z) 1],[Inf 1 1 1])));

[ax,h1,h2] = plotyy(x/1000,v,xpv/1000,fillnan(pvpos,0));
legend('v','PV','Location','SouthWest');
set(h1,'LineWidth',1.5)
set(h2,'LineWidth',1.5);
axes(ax(2))
hline = line(xpv/1000,fillnan(pvneg,0))
set(hline,'LineWidth',1.5)
set(hline,'LineStyle','--')
set(hline,'Color','k');

% Move subplots into position
set(handles.sub1,'Units','pixels')
set(handles.sub2,'Units','pixels')

pos1 = get(handles.sub1,'Position');
set(handles.sub1,'Position',pos1 + [0 -100 0 100])
set(handles.sub2,'Position', pos1 + [0 -400 0 -100]);

% mimic the axis above
xtick = get(handles.sub1,'XTick');
xlim = get(handles.sub1,'XLim');
set(ax(1),'XTick',xtick);
set(ax(2),'XTick',[]);
set(ax(1),'XLim',xlim);
set(ax(2),'XLim',xlim);

set(ax(1),'FontSize',20);
set(ax(2),'FontSize',20);
set(ax(1),'YLim',[-0.5 0]);
set(ax(2),'YLim',[-2 3]*1e-10);
set(ax(2),'ytick',[-1 0 1]*1e-10);
set(ax(1),'YTick',[min(v) 0]);

set(get(ax(1),'ylabel'),'String','v (m/s)');
set(get(ax(2),'ylabel'),'String','PV (m^{-1} s^{-1})');
set(get(ax(1),'ylabel'),'FontSize',20);
set(get(ax(2),'ylabel'),'FontSize',20);

xlabel('X (km)', 'FontSize', 20);
set(handles.sub1,'FontSize', 20);
set(handles.sub2,'FontSize', 20);

export_fig('-painters',[imdir 'initial-cross.png']);

%% Initial conditions cross section - OLD

dir = ['E:\Work\instability\ROMS\si_part\edge\3D\run01\'];
file = [dir 'ocean_pv.nc'];
file2 = [dir 'ocean_his.nc'];
clear handles

temp = ncread(file2,'temp',[1 1 1 1], [Inf Inf Inf 1]);
roms_grid = roms_get_grid(file2,file2,0,1);

hfig2 = figure;
set(hfig2,'Units','pixels');
set(hfig2,'Position',[0 0 1600 900]);
handles.sub1 = subplot(211)
[c2,handles.contourf] = contourf(roms_grid.x_rho(1,:)/1000,roms_grid.z_r(:,1,1),squeeze(temp(:,ceil(size(temp,2)/2),:))',30); shading flat;
fancy_map = flipud(cbrewer('div', 'RdYlGn', 32));
colormap(fancy_map);
title('Temperature','FontSize',24);
handles.h_cbar = colorbar;
set(get(handles.h_cbar,'title'),'String', 'Celsius');
set(get(handles.h_cbar,'Title'),'FontSize',20)
set(gca,'YTick',[-45 -30 -15]);
set(handles.h_cbar,'YTick',[28:2:40]);
ylabel('Z (m)', 'FontSize',20);
handles.lines = linex([8 70-8],'','w');
set(handles.lines,'LineStyle','-')
set(handles.lines,'LineWidth',2)

pv = ncread(file,'pv',[1 1 1 1],[Inf Inf Inf 1]);
xpv = ncread(file,'x_pv');
zpv = ncread(file,'z_pv');

pvmid = squeeze(pv(:,size(pv,2)/2,ceil(size(pv,3)/2),1));
pv0l = xpv(find_approx(pvmid(1:xmid),0))/1000;
pv0r = xpv(find_approx(pvmid(xmid:end),0) + xmid-1)/1000;
xmid = ceil(size(pvmid,1)/2);

pvnl = xpv(find_approx(pvmid(1:xmid),pvmid(xmid),1))/1000;
pvnr = xpv(find(pvmid(xmid:end) == pvmid(xmid),1,'last') + xmid-1)/1000;
pvpos = pvmid; pvpos(pvmid<0)=0;
pvneg = pvmid; pvneg(pvmid>0) = 0;

figure(hfig2)
handles.sub2 = subplot(212)
hold on
plot(xpv/1000,fillnan(pvpos,0),'Color',[0 0.5 0],'LineWidth',1.5);
plot(xpv/1000,fillnan(pvneg,0),'k--','LineWidth',1.5);
%linex([pv0l pv0r],'');
set(handles.sub1,'Units','pixels')
set(handles.sub2,'Units','pixels')

pos1 = get(handles.sub1,'Position');
set(handles.sub1,'Position',pos1 + [0 -100 0 100])
set(handles.sub2,'Position', pos1 + [0 -400 0 -100]);
title('PV (m^{-1} s^{-1})','FontSize',24);
ylim([-2e-10 2e-10]);
xtick = get(gca,'XTick');
xtick(3) = str2num(num2str(pv0l,'%.2f'));
xtick(6) = str2num(num2str(pv0r,'%.2f'));
xtick(end+1) = str2num(num2str(pvnl,'%.2f'));
xtick(end+1) = str2num(num2str(pvnr,'%.2f'));
xtick = sort(xtick);
set(gca,'XTick',xtick);

xlabel('X (km)', 'FontSize', 20);
set(handles.sub1,'FontSize', 20);
set(handles.sub2,'FontSize', 20);
xtick = get(handles.sub1,'XTick');
xtick(3) = str2num(num2str(pv0l,'%.2f'));
xtick(6) = str2num(num2str(pv0r,'%.2f')); 
xtick(end+1) = str2num(num2str(pvnl,'%.2f'));
xtick(end+1) = str2num(num2str(pvnr,'%.2f'));
xtick = sort(xtick);
set(handles.sub1,'XTick',xtick);

set(get(handles.h_cbar,'title'),'String', 'Celsius');
set(get(handles.h_cbar,'Title'),'FontSize',20)

export_fig('-painters',[imdir 'initial-cross.png']);

%% Velocity vectors over temp. 2D run

dir = 'E:\Work\instability\ROMS\si_part\edge\3D\run01-2D\';

file = [dir 'ocean_his.nc'];
tindex = 10;

%for tindex = 10:10
    u = squeeze(double(ncread(file,'u',[1 3 1 tindex],[Inf 1 Inf 1])));
    w = squeeze(double(ncread(file,'w',[1 3 1 tindex],[Inf 1 Inf 1])));

    w = avg1(w(2:end-1,:),2);
    u = avg1(u,1);

    [xr,yr,zr,~,~,~] = roms_var_grid(file,'rho');
    xr = xr(2:end-1);

    range = 205:220;
    handles = mod_movie(file,'temp',[tindex tindex],{'x' range(1) range(end)+2},'y','mid','shading flat; fancy_cmap');

    figure(gcf);
    set(gcf,'renderer','zbuffer');
    hold on
    
    temp = ncread(file,'temp',[1 1 1 1],[Inf Inf Inf 1]);

    [xmat,zmat] = meshgrid(xr(range),zr);
    quiver(xmat/1000,zmat,u(range,:)',w(range,:)'*1000,0,'Color','k','MarkerSize',0.5);
    contour(xmat/1000,zmat,squeeze(temp(range,3,:))',10,'k-','LineWidth',1.5);
    axis tight
    %pause
%end

set(gcf,'Position',[0 0 1600 900])
 %export_fig('-nocrop',[imdir 'quiver-uw.png']);
 
 %% velocity vectors y-z

%% temperature animation map - 2D

cd('E:\Work\instability\ROMS\si_part\edge\3D\run01-2D\');
 mod_movie('ocean_his.nc','temp',[],{'x' '30000' '40000'},'y','mid','fancy_cmap;movieman;');
 
 %% tilted rolls movie - temp
 
 cd('E:\Work\instability\ROMS\si_part\edge\3D\run01\');
 
 mod_movie('ocean_his.nc','temp',[],{'x' '8000' '62000'},'z','mid','fancy_cmap;movieman;');
 system(['move mm_output.avi ' imdir  '\3D-temp.mpg']);
 
 %% u- z- movie
 cd('E:\Work\instability\ROMS\si_part\edge\3D\run01\');
 
 mod_movie('ocean_his.nc','u',[2 17],{'x' '8000' '62000'},'z','mid','fancy_cmap;movieman;caxis([-0.05 0.05])');
 system(['move mm_output.avi ' imdir  '\3D-u-z.mpg']);
 
 %% u - snapshot
 cd('E:\Work\instability\ROMS\si_part\edge\3D\run01\');
 
 mod_movie('ocean_his.nc','u',[14 14],{'x' '8000' '62000'},'z','mid','fancy_cmap;topresent');
 set(gcf,'renderer','zbuffer');
 export_fig([imdir '/u-snapshot.png']);
 %system(['move mm_output.avi ' imdir  '\3D-u.mpg']);
 
 %% u - snapshot + initial PV
 
 dir = ['E:\Work\instability\ROMS\si_part\edge\3D\run01\'];
file = [dir 'ocean_pv.nc'];
file2 = [dir 'ocean_his.nc'];
clear handles

tindex = 14;

temp = ncread(file2,'u',[1 1 10 tindex], [Inf Inf 1 1]);
roms_grid = roms_get_grid(file2,file2,0,1);
pv = ncread(file,'pv',[1 1 1 1],[Inf Inf Inf 1]);
xpv = ncread(file,'x_pv');
zpv = ncread(file,'z_pv');
xmid = ceil(size(pv,1)/2);

pvmid = squeeze(pv(:,size(pv,2)/2,ceil(size(pv,3)/2),1));
pv0l = xpv(find_approx(pvmid(1:xmid),0))/1000;
pv0r = xpv(find_approx(pvmid(xmid:end),0) + xmid-1)/1000;

pvnl = 20.88;%xpv(find_approx(pvmid(1:xmid),pvmid(xmid),1))/1000;
pvnr = 49.3;%xpv(find(pvmid(xmid:end) == pvmid(xmid),1,'last') + xmid-1)/1000;
pvpos = pvmid; pvpos(pvmid<0)=0;
pvneg = pvmid; pvneg(pvmid>0) = 0;

% plot u
hfig2 = figure;
set(hfig2,'Units','pixels');
set(hfig2,'Position',[0 0 1600 900]);
handles.sub1 = subplot(211);
[c2,handles.contourf] = contourf(roms_grid.x_u(1,:)/1000,roms_grid.y_u(:,1,1)/1000,squeeze(temp(:,:))',30); shading flat;
fancy_map = flipud(cbrewer('div', 'RdYlGn', 32));
colormap(fancy_map);
handles.h_cbar = colorbar;

title('u','FontSize',24);

set(get(handles.h_cbar,'title'),'String', 'm/s');
set(get(handles.h_cbar,'Title'),'FontSize',20)
set(handles.h_cbar,'YTick',[-0.05:0.01:0.05]);

set(gca,'YTick',[0:25:150]);
ylabel('Y (km)', 'FontSize',20);
xtick = get(handles.sub1,'XTick');
xtick(3) = str2num(num2str(pv0l,'%.2f'));
xtick(6) = str2num(num2str(pv0r,'%.2f')); 
xtick(end+1) = str2num(num2str(pvnl,'%.2f'));
xtick(end+1) = str2num(num2str(pvnr,'%.2f'));
xtick = sort(xtick);
set(handles.sub1,'XTick',xtick);
    
% Now do PV , V

figure(hfig2)
handles.sub2 = subplot(212)
hold on

[x,y,z,~,~,~] = roms_var_grid(file2,'v');
v = double(squeeze(ncread(file2,'v',[1 ceil(length(y)/2) length(z) 1],[Inf 1 1 1])));

[ax,h1,h2] = plotyy(x/1000,v,xpv/1000,fillnan(pvpos,0));
legend('v','PV','Location','SouthWest');
set(h1,'LineWidth',1.5)
set(h2,'LineWidth',1.5);
axes(ax(2))
hline = line(xpv/1000,fillnan(pvneg,0))
set(hline,'LineWidth',1.5)
set(hline,'LineStyle','--')
set(hline,'Color','k');

% Move subplots into position
set(handles.sub1,'Units','pixels')
set(handles.sub2,'Units','pixels')

pos1 = get(handles.sub1,'Position');
set(handles.sub1,'Position',pos1 + [0 -100 0 100])
set(handles.sub2,'Position', pos1 + [0 -400 0 -100]);

% mimic the axis above
xtick = get(handles.sub1,'XTick');
xlim = get(handles.sub1,'XLim');
set(ax(1),'XTick',xtick);
set(ax(2),'XTick',[]);
set(ax(1),'XLim',xlim);
set(ax(2),'XLim',xlim);

set(ax(1),'FontSize',20);
set(ax(2),'FontSize',20);
set(ax(1),'YLim',[-0.5 0]);
set(ax(2),'YLim',[-2 3]*1e-10);
set(ax(2),'ytick',[-1 0 1]*1e-10);
set(ax(1),'YTick',[min(v) 0]);

set(get(ax(1),'ylabel'),'String','v (m/s)');
set(get(ax(2),'ylabel'),'String','PV (m^{-1} s^{-1})');
set(get(ax(1),'ylabel'),'FontSize',20);
set(get(ax(2),'ylabel'),'FontSize',20);

xlabel('X (km)', 'FontSize', 20);
set(handles.sub1,'FontSize', 20);
set(handles.sub2,'FontSize', 20);

export_fig('-painters',[imdir 'u-snapshot-pv-v.png']);

%% u - snapshot + Ri,Ro

 dir = ['E:\Work\instability\ROMS\si_part\edge\3D\run01\'];
file = [dir 'ocean_pv.nc'];
file2 = [dir 'ocean_his.nc'];
clear handles

tindex = 9; % time of max. growth rate

temp = ncread(file2,'u',[1 240 1 tindex], [Inf 1 Inf 1]);
roms_grid = roms_get_grid(file2,file2,0,1);
pv = ncread(file,'pv',[1 1 1 1],[Inf Inf Inf 1]);
xpv = ncread(file,'x_pv');
zpv = ncread(file,'z_pv');
xmid = ceil(size(pv,1)/2);

pvmid = squeeze(pv(:,size(pv,2)/2,ceil(size(pv,3)/2),1));
pv0l = xpv(find_approx(pvmid(1:xmid),0))/1000;
pv0r = xpv(find_approx(pvmid(xmid:end),0) + xmid-1)/1000;

pvnl = 20.88;%xpv(find_approx(pvmid(1:xmid),pvmid(xmid),1))/1000;
pvnr = 49.3;%xpv(find(pvmid(xmid:end) == pvmid(xmid),1,'last') + xmid-1)/1000;
pvpos = pvmid; pvpos(pvmid<0)=0;
pvneg = pvmid; pvneg(pvmid>0) = 0;

% plot u
hfig2 = figure;
set(hfig2,'Units','pixels');
set(hfig2,'Position',[0 0 1600 900]);
handles.sub1 = subplot(211);
[c2,handles.contourf] = contourf(roms_grid.x_u(1,:)/1000,roms_grid.z_u(:,1,1),squeeze(temp(:,:))',30); shading flat;
fancy_map = flipud(cbrewer('div', 'RdYlGn', 32));
colormap(fancy_map);
handles.h_cbar = colorbar;

title('u','FontSize',24);

set(get(handles.h_cbar,'title'),'String', 'm/s');
set(get(handles.h_cbar,'Title'),'FontSize',20)
set(handles.h_cbar,'YTick',[-0.05:0.01:0.05]);

set(gca,'YTick',[-50:10:0]);
ylabel('Z (m)', 'FontSize',20);
xtick = get(handles.sub1,'XTick');
xtick(3) = str2num(num2str(pv0l,'%.2f'));
xtick(6) = str2num(num2str(pv0r,'%.2f')); 
xtick(end+1) = str2num(num2str(pvnl,'%.2f'));
xtick(end+1) = str2num(num2str(pvnr,'%.2f'));
xtick = sort(xtick);
set(handles.sub1,'XTick',xtick);

% Now do Ri, Ro

figure(hfig2)
handles.sub2 = subplot(212)
hold on

[x,y,z,~,~,~] = roms_var_grid(file2,'v');

zindex = 10;

dx = x(2)-x(1); 
f0 = 1e-4;
N2 = 1e-5; %%%%%%%%%%%%%%%%%%%%%%%%% HARD CODED VALUE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = double(squeeze(ncread(file2,'v',[1 ceil(length(y)/2) 1 1],[Inf 1 Inf 1])));
Ro = (diff(v(:,zindex),1,1)./dx./f0);
vz = bsxfun(@rdivide,diff(v,1,2),diff(z)');
Ri = N2./vz.^2;

plotRi = addnan(Ri(:,zindex),2);
sigma = f0^2 * (avg1(1./plotRi) - 1 -Ro);

[ax,h1,h2] = plotyy(avg1(x)/1000,Ro,x/1000,plotRi);
legend('Ro','Ri','Location','SouthEast');
set(h1,'LineWidth',1.5)
set(h2,'LineWidth',1.5);
% axes(ax(2))
% hline = line(xpv/1000,fillnan(pvneg,0))
% set(hline,'LineWidth',1.5)
% set(hline,'LineStyle','--')
% set(hline,'Color','k');

% Move subplots into position
set(handles.sub1,'Units','pixels')
set(handles.sub2,'Units','pixels')

pos1 = get(handles.sub1,'Position');
set(handles.sub1,'Position',pos1 + [0 -100 0 100])
set(handles.sub2,'Position', pos1 + [0 -400 0 -100]);

% mimic the axis above
xtick = get(handles.sub1,'XTick');
xlim = get(handles.sub1,'XLim');
set(ax(1),'XTick',xtick);
set(ax(2),'XTick',[]);
set(ax(1),'XLim',xlim);
set(ax(2),'XLim',xlim);

set(ax(1),'FontSize',20);
set(ax(2),'FontSize',20);
set(ax(1),'YLim',[-0.5 0.5]);
set(ax(2),'YLim',[0 2]);
set(ax(2),'ytick',[0 roundoff(min(plotRi),2) 2]);
set(ax(1),'YTick',[roundoff(min(Ro),2) 0 roundoff(max(Ro),2)]);

set(get(ax(1),'ylabel'),'String','Ro');
set(get(ax(2),'ylabel'),'String','Ri');
set(get(ax(1),'ylabel'),'FontSize',20);
set(get(ax(2),'ylabel'),'FontSize',20);

xlabel('X (km)', 'FontSize', 20);
set(handles.sub1,'FontSize', 20);
set(handles.sub2,'FontSize', 20);

%export_fig('-painters',[imdir 'u-ri-ro.png']);

%% u-snapshot + sigma
 dir = ['E:\Work\instability\ROMS\si_part\edge\3D\run01\'];
file = [dir 'ocean_pv.nc'];
file2 = [dir 'ocean_his.nc'];
clear handles

tindex = 9; % time of max. growth rate

temp = ncread(file2,'u',[1 240 1 tindex], [Inf 1 Inf 1]);
roms_grid = roms_get_grid(file2,file2,0,1);
pv = ncread(file,'pv',[1 1 1 1],[Inf Inf Inf 1]);
xpv = ncread(file,'x_pv');
zpv = ncread(file,'z_pv');
xmid = ceil(size(pv,1)/2);

pvmid = squeeze(pv(:,size(pv,2)/2,ceil(size(pv,3)/2),1));
pv0l = xpv(find_approx(pvmid(1:xmid),0))/1000;
pv0r = xpv(find_approx(pvmid(xmid:end),0) + xmid-1)/1000;

pvnl = 20.88;%xpv(find_approx(pvmid(1:xmid),pvmid(xmid),1))/1000;
pvnr = 49.3;%xpv(find(pvmid(xmid:end) == pvmid(xmid),1,'last') + xmid-1)/1000;
pvpos = pvmid; pvpos(pvmid<0)=0;
pvneg = pvmid; pvneg(pvmid>0) = 0;

% plot u
hfig2 = figure;
set(hfig2,'Units','pixels');
set(hfig2,'Position',[0 0 1600 900]);
handles.sub1 = subplot(211);
[c2,handles.contourf] = contourf(roms_grid.x_u(1,:)/1000,roms_grid.z_u(:,1,1),squeeze(temp(:,:))',30); shading flat;
fancy_map = flipud(cbrewer('div', 'RdYlGn', 32));
colormap(fancy_map);
handles.h_cbar = colorbar;

title('u','FontSize',24);

set(get(handles.h_cbar,'title'),'String', 'm/s');
set(get(handles.h_cbar,'Title'),'FontSize',20)
set(handles.h_cbar,'YTick',[-0.05:0.01:0.05]);

set(gca,'YTick',[-50:10:0]);
ylabel('Z (m)', 'FontSize',20);
xtick = get(handles.sub1,'XTick');
xtick(3) = str2num(num2str(pv0l,'%.2f'));
xtick(6) = str2num(num2str(pv0r,'%.2f')); 
xtick(end+1) = str2num(num2str(pvnl,'%.2f'));
xtick(end+1) = str2num(num2str(pvnr,'%.2f'));
xtick = sort(xtick);
set(handles.sub1,'XTick',xtick);

% Now do Ri, Ro

figure(hfig2)
handles.sub2 = subplot(212)
hold on

[x,y,z,~,~,~] = roms_var_grid(file2,'v');

zindex = 10;

dx = x(2)-x(1); 
f0 = 1e-4;
N2 = 1e-5; %%%%%%%%%%%%%%%%%%%%%%%%% HARD CODED VALUE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = double(squeeze(ncread(file2,'v',[1 ceil(length(y)/2) 1 1],[Inf 1 Inf 1])));
Ro = (diff(v(:,zindex),1,1)./dx./f0);
vz = bsxfun(@rdivide,diff(v,1,2),diff(z)');
Ri = N2./vz.^2;

plotRi = addnan(Ri(:,zindex),2);
sigma = real(sqrt(f0^2 * (avg1(1./plotRi) - 1 -Ro)))*86400;

hplot = plot(avg1(x)/1000,sigma,'LineWidth',1.5);
ax(1) = gca;
% Move subplots into position
set(handles.sub1,'Units','pixels')
set(handles.sub2,'Units','pixels')

pos1 = get(handles.sub1,'Position');
set(handles.sub1,'Position',pos1 + [0 -100 0 100])
set(handles.sub2,'Position', pos1 + [0 -400 0 -100]);

% mimic the axis above
xtick = get(handles.sub1,'XTick');
xlim = get(handles.sub1,'XLim');
set(gca,'XTick',xtick);
set(gca,'XLim',xlim);

set(ax(1),'FontSize',20);
set(ax(1),'YLim',[0 10]);
set(ax(1), 'ytick',sort([get(ax(1),'ytick') roundoff(max(sigma),2)]));

ylabel('\sigma_{max} (per day)');

xlabel('X (km)', 'FontSize', 20);
set(handles.sub1,'FontSize', 20);
set(handles.sub2,'FontSize', 20);

export_fig('-painters',[imdir 'u-sigma.png']);

%% no-shear v profile

load v-fullshear

 dir = ['E:\Work\instability\ROMS\si_part\edge\3D\run01\'];
file2 = [dir 'ocean_his.nc'];

vel = v(:,240,end);

plot(xvmat(:,1,end)/1000,vel,'LineWidth',2);
set(gca,'FontSize',20);
xlabel('X (km)');
ylabel('v (m/s)');
xlim([0 70]);
pause;
export_fig([imdir 'v-fullshear.png']);
 
 %% baroclinic instability - before & after temp
 
 mod_movie('ocean_his.nc','temp',[1 1],{'x' '8000' '62000'},'y','mid','nocaxis;fancy_cmap');
 set(gcf,'Position',[0 0 1600 900])
 export_fig('-nocrop',[imdir 'temp-before.png']);
 mod_movie('ocean_his.nc','temp',[Inf Inf],{'x' '8000' '62000'},'y','mid','nocaxis;fancy_cmap');
 set(gcf,'Position',[0 0 1600 900])
 export_fig('-nocrop',[imdir 'temp-after.png']);
 
 %% energy diagnostics
 
 cd('E:\Work\instability\ROMS\si_part\edge\3D\run01');
 
 load energy-avg-y-full.mat
 
 fontSize = 20;
 
 figure;
 
 subplot(311)
 plot(time_A./86400,A*86400,'LineWidth',2);
 title('Growth Rate (per day)','FontSize',20);
 set(gca,'FontSize',fontSize);
 subplot(312)
 plot(t_en./86400,EKE,'b','LineWidth',2); 
 title('EKE per unit mass','FontSize',20);
 set(gca,'FontSize',fontSize);
 subplot(313)
 plot(t_en./86400,PE,'b','LineWidth',2);
 title('PE per unit mass','FontSize',20);
 xlabel('Time (days)','FontSize',fontSize);
 set(gca,'FontSize',fontSize);
 
 set(gcf,'Position',[0 0 1600 900])
 export_fig('-painters',[imdir 'energy-diag.png']);
 
 %% Compare 2D vs 3D temp
 
 for ii=1:31
    mod_movie('E:\Work\instability\ROMS\si_part\edge\3D\run01\ocean_his.nc','temp',[ii ii],{},'y','mid','pcolor;shading interp;fancy_cmap;nocaxis');
    mod_movie('E:\Work\instability\ROMS\si_part\edge\3D\run01-2D\ocean_his.nc','temp',[ii ii],{},'y','mid','pcolor;shading interp;fancy_cmap;nocaxis');
    pause;
    close all
 end
 
 %% bfric snapshots - look in analyze.m
 
 %% 3D slices

varname = 'u';
fname = 'ocean_his.nc';

xloc = 15; yloc = 40; % km
zloc = -10; %m
tlim = 20;

volume = {'x' num2str(xloc*1000) '55000';
          'y' num2str(yloc*1000) Inf;
          'z' '-50' num2str(zloc)};

[x,y,z,vol] = roms_extract(fname,varname,volume);
[read_start,read_count] = roms_ncread_params(4,0,1,100,[1 tlim],1,vol);
var = ncread(fname,varname,read_start,read_count);
[xmat,ymat,zmat] = meshgrid(x/1000,y/1000,z);

time = ncread(fname,'ocean_time');

make_video = 1;

if make_video
    mm_instance = mm_setup;
    mm_instance.pixelSize = [1600 900];
    mm_instance.outputFile = 'mm_output.avi';
    mm_instance.ffmpegArgs = '-q:v 1 -g 1';
    mm_instance.InputFrameRate = 3;
    mm_instance.frameRate = 3;
end

hfig = figure; set(gcf,'Renderer','opengl');set(hfig,'Position',[0 0 1600 900]);
fancy_map = flipud(cbrewer('div', 'RdYlGn', 32));
for i=2:tlim
    var1=permute(var(:,:,:,i),[2 1 3]);
    h = slice(xmat,ymat,zmat,double(var1),xloc,yloc,z(end));
    set(h,'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
    set(gca,'xlim',[xloc max(x)]); 
    set(gca,'zlim',[-50 zloc]);
    set(gca,'ztick',[-45:15:-15]);
    axis tight;
    xlabel('X (km)','FontSize', 20);
    ylabel('Y (km)','FontSize', 20);
    zlabel('Z (m)','FontSize', 20);
    title(sprintf('u | t = %.2f / %.2f days',time(i)/86400,time(end)/86400),'FontSize',20);
    set(gca,'FontSize',20);
    view(-37.5,50);
    
    if make_video, caxis([-0.05 0.05]); end
    
    colormap(fancy_map);
    h_cbar = colorbar;
    set(get(h_cbar,'title'),'String', 'm/s');
    set(get(h_cbar,'Title'),'FontSize',20)
    
    set(h_cbar,'YTick',[-0.05:0.01:0.05]);
    set(h_cbar,'FontSize',20);
    
    
    if make_video
        mm_addFrame(mm_instance,gcf); 
    else 
        pause(0.1);
    end
end

if make_video
    mm_render(mm_instance);
    system(['move mm_output.avi ' imdir '\slice-u.mpg']);
end

%% random transparency (alpha) code

%export_fig([imdir 'tide-temps.png']);
% set(gcf,'Renderer','opengl');
% alphable = findobj(h_plot, '-property', 'FaceAlpha');
% for i=1:length(alphable)
%  set(alphable, 'FaceAlpha', 0.6);
% end
% cb = findobj(gcf,'Type','axes','Tag','Colorbar');
% cbIm = findobj(cb,'Type','image');
% alpha(cbIm,0.8)

%set(h_plot,'facecolor', 'interp', 'FaceAlpha',0.8) ;
%alpha(h_plot,0.8);

%% random isosurface / patch code

%s
% 
%     p= patch(isosurface(xmat,ymat,zmat,u1,20));
%     isonormals(xmat,ymat,zmat,u1, p)
%     set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
%     daspect([1 1 1]); axis tight; 
%     view([-190.5 -10])
%     colormap(prism(28))
%     camlight; lighting gouraud
%     pause(0.01)
% end


 
 