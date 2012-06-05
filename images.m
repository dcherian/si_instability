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

v = ncread(file2,'v',[1 1 1 1], [Inf Inf Inf 1]);
roms_grid = roms_get_grid(file2,file2,0,1);
hold on
[c2,h2] = contour(roms_grid.x_v(1,:)/1000,roms_grid.z_v(:,1,1),squeeze(v(:,ceil(size(v,2)/2),:))','b','LineWidth',1.5);
clabel(c2,h2,[-0.35:0.1:0],'FontSize',16,'Color','b','LabelSpacing',72*6,'Rotation',0);

pv = ncread(file,'pv',[1 1 1 1],[Inf Inf Inf 1]);
xpv = ncread(file,'x_pv');
zpv = ncread(file,'z_pv');
pvmid = squeeze(pv(:,size(pv,2)/2,:,1));
pvpos = pvmid; pvpos(pvmid<0)=0;
pvneg = pvmid; pvneg(pvmid>0) = 0;
% figure
% [c,h] = contour(xpv,zpv,pvpos',20,'k');
% %clabel(c,h);
hold on
[c,h] = contour(xpv/1000,zpv,pvneg',10,'k','LineWidth',1.5);
set(h,'LineWidth',1.5);
%clabel(c,h);

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
%handles = mod_movie(file2,'temp',[1 1],{},'y','mid','shading flat; fancy_cmap; topresent');

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
xmid = ceil(size(pvmid,1)/2);

pvmid = squeeze(pv(:,size(pv,2)/2,ceil(size(pv,3)/2),1));
pv0l = xpv(find_approx(pvmid(1:xmid),0))/1000;
pv0r = xpv(find_approx(pvmid(xmid:end),0) + xmid-1)/1000;

pvnl = xpv(find_approx(pvmid(1:xmid),pvmid(xmid),1))/1000;
pvnr = xpv(find(pvmid(xmid:end) == pvmid(xmid),1,'last') + xmid-1)/1000;
pvpos = pvmid; pvpos(pvmid<0)=0;
pvneg = pvmid; pvneg(pvmid>0) = 0;

figure(hfig2)
handles.sub2 = subplot(212)
hold on
plot(xpv/1000,fillnan(pvpos,0),'b','LineWidth',1.5);
plot(xpv/1000,fillnan(pvneg,0),'k--','LineWidth',1.5);
%linex([pv0l pv0r],'');
set(handles.sub1,'Units','pixels')
set(handles.sub2,'Units','pixels')

pos1 = get(handles.sub1,'Position');
set(handles.sub1,'Position',pos1 + [0 -100 0 100])
set(handles.sub2,'Position', pos1 + [0 -400 0 -100]);
title('PV','FontSize',24);
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

%% temperature animation map - 2D

cd('E:\Work\instability\ROMS\si_part\edge\3D\run01-2D\');
 mod_movie('ocean_his.nc','temp',[],{'x' '30000' '40000'},'y','mid','fancy_cmap;movieman;');
 
 %%  tilted rolls movie - temp
 
 cd('E:\Work\instability\ROMS\si_part\edge\3D\run01\');
 
 mod_movie('ocean_his.nc','temp',[],{'x' '8000' '62000'},'z','mid','fancy_cmap;movieman;');
 system(['move mm_output.avi ' imdir  '\3D-temp.mpg']);
 
 %% tilted rolls movie - u
 cd('E:\Work\instability\ROMS\si_part\edge\3D\run01\');
 
 mod_movie('ocean_his.nc','u',[2 20],{'x' '8000' '62000'},'z','mid','fancy_cmap;movieman;caxis([-0.05 0.05])');
 system(['move mm_output.avi ' imdir  '\3D-u.mpg']);
 
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
 
 fontSize = 16;
 
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
 export_fig([imdir 'energy-diag.png']);
 
 %% Compare 2D vs 3D temp
 
 for ii=1:31
    mod_movie('E:\Work\instability\ROMS\si_part\edge\3D\run01\ocean_his.nc','temp',[ii ii],{},'y','mid','pcolor;shading interp;fancy_cmap;nocaxis');
    mod_movie('E:\Work\instability\ROMS\si_part\edge\3D\run01-2D\ocean_his.nc','temp',[ii ii],{},'y','mid','pcolor;shading interp;fancy_cmap;nocaxis');
    pause;
    close all
 end
 