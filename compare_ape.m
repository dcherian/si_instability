% compares APE for 2D runs
dir = 'E:\Work\instability\ROMS\si_part\edge\2D\';
dirs = {'run01','run02','run03','run04','run05','run06','run07','run08','run09','run10_2','run11','run12_2','run13_2','run14_2','run15', ...
        'run16','run17','run18','run19','run21','run23','run24'}; % 4 and 7 are outliers (PVmin / PVmid > 1.1 - greater wavelength too)
runx = [01 02 03 04 05 06 07 08 09 10.2 11 12.2 13.2 14.2 15 16 17 18 19 21 23 24];
fname = 'ocean_his.nc';
volume = {};

enname = 'energy-avg-x.mat';
pvname = 'ocean_pv.nc';

redo_en = 1;
flag2D = 1;

for ii=1:length(dirs)
    cd([dir dirs{ii}]);
    
    % find APE
    all_ape(ii) = roms_ape(fname);
    
    if flag2D % then compare final PE - Initial PE
        if ~exist(enname,'file') || redo_en == 1, roms_energy(fname,[],volume,4,1); close all; end
                
        load(enname);
        
        dPE(ii) = PE(end)-PE(1);
    else % fr 3D runs
        % figure out SI and BI peaks in growth rate
        % choose average of two as endpoint for SI
    end
    
    % figure out PVmid 
    if ~exist(pvname,'file')
        pv = roms_pv(fname,[1 1]); 
    else
        pv = ncread(pvname,'pv',[1 1 1 1],[Inf Inf Inf 1]);
    end
    
    s = size(pv);
    pvmid(ii) = pv(ceil(s(1)/2),ceil(s(2)/2),ceil(s(3)/2));
end

%% plot

figure(1);
subplot(211)
plot(abs(all_ape),abs(dPE),'*');
xlabel('APE');
ylabel('\Delta PE');
beautify; box off;

subplot(212)
plot(abs(all_ape),abs(dPE./all_ape),'*');
xlabel('APE (J/kg)');
ylabel('\Delta PE / APE');
beautify; box off;