% figures out aspect ratio / slope of SI cells

gcm = 1;

gcm_dir  = 'E:\Work\instability\mitgcm\2d\results\';
roms_dir = 'E:\Work\instability\ROMS\si_full\2D cases\';

roms_index = [10:21];
gcm_index  = [62:65 68:76];
gcm_visc = [4 4 3 2 10 20 0.9 2.4 2 2 2 2 4];
roms_visc = [5 5 5 5 5 5 5 2 1 2 3 4];
lgcm = length(gcm_index);
lroms = length(roms_index);

clear Ar

%%

warning off

dx  = 80;
dz  = 2;
xjump = 10;
zjump = 5;

if gcm, limit = lgcm; index = gcm_index; visc = gcm_visc;
else    limit = lroms; index = roms_index; visc = roms_visc; end

figure

cpb = progressbar();
for kk=1:limit
    
    if gcm
        fname = [gcm_dir 'mnc_00' num2str(gcm_index(kk)) '\state*'];
        vars = rdmnc(fname);
        txz = squeeze(vars.Temp(:,1,:,end));
        dx = mean(diff(vars.Xp1));
        dz = mean(diff(vars.Z));
    else
        fname = [roms_dir 'run' num2str(roms_index(kk)) '\ocean_his.nc'];
        %xr = double(squeeze(ncread(fname,'x_rho');
        grid = roms_get_grid(fname,fname,0,1);
        txz = double(squeeze(ncread(fname,'temp',[1 3 1 1],[Inf 1 Inf Inf])));
        dx = mean(mean(diff(grid.x_rho,1,2)));
        dz = 2;
        txz = txz(:,:,end);
        [x,z] = meshgrid(grid.z_r(:,1,1),grid.x_rho(1,:));
        [xi,zi] = meshgrid([grid.z_r(1):dz:grid.z_r(end)],grid.x_rho(1,:));
        txz=interp2(x,z,txz,xi,zi);
    end
    s   = size(txz);

    clear covx covz lx lz

    for i=1:zjump:s(2)
        ind = ceil(i/zjump);
        [covx(ind,:),lagx] = xcov(txz(:,i),'biased');
        lx(ind) = mean(abs(lagx(find_approx(covx(ind,:),0,3))*dx));
    end

    for i=1:xjump:s(1)
        ind = ceil(i/xjump);
        [covz(ind,:),lagz] = xcov(txz(i,:),'biased');
        lz(ind) = mean(abs(lagz(find_approx(covz(ind,:),0,3))*dz));
    end

    Ar(kk) = median(lx)./median(lz);
    txt = sprintf(' Progress: kk=%d',kk);
    progressbarupdate(cpb,kk/limit*100,txt);
end
cpb.stop();

%%
plot(visc,Ar,'*')
text(visc,Ar,num2str(index'));
if gcm, title('MITgcm'); else title('ROMS'); end
xlabel('visc');
ylabel('Ar = L_x/L_z');
xlim([0 10]);

warning on