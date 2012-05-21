%% set variables

imdir = 'E:\Work\instability\graphs\talk\';

%% Ken's run cross section - u
dir = 'E:\Work\instability\ken102\';
file = [dir 'STABB102_his_0001.nc'];
mod_movie(file,'u',[10 10],{},'y','mid','pcolor;shading interp;nocaxis;fancy_cmap');

%% Ken's run cross section - pv
file = [dir 'ocean_pv.nc'];
file2 = [dir 'STABB102_his_0001.nc'];
h_plot = mod_movie(file2,'temp',[1 1],{},'y','mid','shading flat; fancy_cmap');
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
[c2,h2] = contour(roms_grid.x_v(1,:)/1000,roms_grid.z_v(:,1,1),squeeze(v(:,ceil(size(v,2)/2),:))','b');

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
[c,h] = contour(xpv/1000,zpv,pvneg',10,'k');
set(h,'LineWidth',1.5);
%clabel(c,h);

export_fig([imdir 'ken-cross.png']);