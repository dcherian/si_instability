% plots for bunch of edge effects runs

folder = 'E:\Work\instability\ROMS\si_part\edge\';
%runs = {'run02','run06','run07','run08'};
%aa = [1.1 1.5 1.3 1.4];
%runs={'run11','run13','run12','run14','run18','run19','run120'}; %'run06','run10' - 30,40
%aa = [50 60 70 100 101 150 120];
runs = {'run50_1','run50_2','run50_3','run50_4','run50_5','run50_6','run50_7''run50_8'};
aa = [32.445040 40.231040 46.231040 32.2310040 40.345040 40.345060 36.455040];
%runs={'run120'}
colors = distinguishable_colors(length(aa));

redo_l = 0;
cutoff = 1/100; % % of SI EKE below which diagnostics are ignored

volume = {'x' '12000' '35000'};

clear Aruns t_Aruns t_Eruns Lxruns Lyruns Lxmax Lymax tmax Amax

figure(1); clf;
figure(2); clf;
figure(3); clf;

tmax = 7 * 86400; % time range to study (in sec). Because some runs go to 20 days or so
fname = 'ocean_his.nc';

for i=1:length(runs)
    dir = [folder runs{i} '\'];
    cd(dir);
    
    if ~exist('length_scales_u.mat','file') || redo_l == 1
        roms_length_scales(fname,'u',[],volume,0,'ignore3;ignore2;'); 
        title(runs{i});
    end
    if ~exist('energy-avg-x.mat','file'), roms_energy(fname,[],4,1,0); end
    
    % read data
    load length_scales_u
    load energy-avg-x
    
    warning off
    roms_grid = roms_get_grid(fname,fname,0,1);
    warning on
    misc = roms_load_misc(fname);
    time = ncread(fname,'ocean_time');
    f0 = mean(misc.f(:));
    g = 9.81;
    
    % growth rates
    tind = find_approx(time_A,tmax,1);
    Aruns{i} = A(1:tind);
    pk = findpeaks(A(1:tind),'SORTSTR','descend');
    try 
        mind = find(A(1:tind) == pk(1));
    catch ME
        mind = find(time_A(1:tind)./86400 == 1.5);
    end
    Amax(i) = A(mind);
    tAmax(i) = time_A(find(A(1:tind) == Amax(i)));
    t_Aruns{i} = time_A(1:tind);
    
    figure(1); subplot(411); hold on;
    plot(time_A(1:tind)./86400,A(1:tind)*86400,'Color',colors(i,:));
    plot(tAmax(i)./86400,Amax(i)*86400,'*','LineWidth',2,'Color',colors(i,:));
    xlim([0 tmax/86400+1]);
    %xlabel('Time (days)');
    ylabel('Growth Rate (d^{-1})');
    
    % EKE
    tind = find_approx(t_en,tmax,1);
    Eruns{i} = EKE(1:tind);
    t_Eruns{i} = t_en(1:tind);
    pk = findpeaks(Eruns{1},'sortstr','descend');
    cutind = find_approx(Eruns{i},EKE(1)+cutoff*pk(1),1);
    t_cut = t_en(cutind);
    
    figure(1); subplot(412); hold on;
    plot(t_en(1:tind)./86400,EKE(1:tind),'Color',colors(i,:));
    xlim([0 tmax/86400+1]);
    %xlabel('Time (days)');
    ylabel('EKE');
    
    % Length Scales
    tind = find_approx(time_L,tmax,1);
    cutind2 = find_approx(time_L,t_cut,1);
    Lxruns{i} = conv(L(1,cutind2:tind),[1 1 1]/3,'same');
    Lyruns{i} = conv(L(2,cutind2:tind),[1 1 1]/3,'same');
    ind = find_approx(time_L,tAmax(i),1);
    Lxmax(i) = L(1,ind);
    Lymax(i) = L(2,ind);
    
    figure(1); subplot(413); hold on;
    plot(time_L(cutind2:tind)./86400,L(1,cutind2:tind)/1000,'Color',colors(i,:));
    plot(time_L(ind)./86400,Lxmax(i)/1000,'*','LineWidth',2,'Color',colors(i,:));
    xlim([0 tmax/86400+1]);
    %xlabel('Time (days)');
    ylabel('L_x(km)');
    
    figure(1); subplot(414); hold on;
    plot(time_L(cutind2:tind)./86400,L(2,cutind2:tind)/1000,'Color',colors(i,:));
    plot(time_L(ind)./86400,Lymax(i)/1000,'*','LineWidth',2,'Color',colors(i,:));
    xlim([0 tmax/86400+1]);
    xlabel('Time (days)');
    ylabel('L_y(km)');
    
    figure(3); hold on
    plot(time_L(cutind2:tind)./86400,L(2,cutind2:tind)/1000,'Color',colors(i,:));
    %plot(time_L(ind)./86400,Lymax(i)/1000,'*','LineWidth',2,'Color',colors(i,:));
    xlim([0 tmax/86400+1]);
    xlabel('Time (days)');
    ylabel('L_y(km)');
    % Tx, vz, vx - Estimate M^2, N^2, 

    dx = mean(mean(diff(roms_grid.x_rho,1,2)));
    dy = mean(mean(diff(roms_grid.y_rho,1,1)));
    dz = squeeze(diff(roms_grid.z_r(:,1,1)));

    xmid = ceil(size(roms_grid.x_rho,2)/2);
    ymid = ceil(size(roms_grid.y_rho,1)/2);
    zmid = ceil(size(roms_grid.z_r  ,1)/2);
    
    temp = double(ncread('ocean_his.nc','temp',[1 1 1 1],[Inf Inf Inf 1]));
    M2 = g*misc.Tcoef*avg1(squeeze(diff(temp(:,ymid,:,:),1,1)),2)./dx;
    N2 = g*misc.Tcoef*avg1(bsxfun(@rdivide,diff(squeeze(temp(:,ymid,:,:)),1,2),permute(dz,[3 1 2])),1);
    
    M2mid = M2(xmid,zmid);
    N2mid = N2(xmid,zmid);
    
    % Estimate theoretical most unstable length scale
    m = 2*pi/100;
    k = (-M2mid)/N2mid * m;
    figure(1); subplot(413);
    liney(2*pi/k/4/1000);
    
    % Calculate slope of isotherms and richardson number using both uz and vz
    
    % Estimate initial APE
    
    % plots
    
end

figure(1); 
subplot(411); title('Summary');
subplot(412); legend(num2str(aa'),'Location','best');

figure(3); legend(num2str(aa'),'Location','best');


%%
figure(2);
subplot(211);
plot(aa,Lxmax/1000,'r*');
hold on
plot(aa,Lymax/1000,'b*');
legend('L_x','L_y','Location','Best');
ylabel('Length scale at max. growth rate');
%xlim([1 2]);
subplot(212);
plot(aa,Amax*86400,'k*');
ylabel('Growth Rate (d^{-1})');
xlabel('Run');
%xlim([1 2]);
%legend('A');

