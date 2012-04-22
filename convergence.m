file  = 'ocean_his.nc';
%dir   = '/home/poison/deepak/ROMS/runs/si/si_full/3d/';
dir = 'E:\Work\instability\ROMS\si_full\3d\';
%dir ='/media/Data/Work/instability/ROMS/si_full/3d/';
list  = {'run41','run42','run43','run44','run46','run47','run48','run50-3','run61','run62'};
%list = {'run50','run50-2','run50-3'};
%list = {'run48','run54'};
vars  = {'energy'; 'length_scales'; 'growth_rate'};
color = {'r','g','b','k','c','m','y','k'};
kencode = 0;
energy_x = 0;
volume = {'x' '14000' '22000'};

clear label exmax px eymax Aymax Axmax LX LY LZ EKEx EKEy T

redo_ex = 0;
redo_ey = 0;
redo_l  = 0;
figure(1); clf; hold on; title('EKE/Area');
figure(2); clf; hold on; title('Growth Rates');

if kencode, engyfname = 'energy-ken.mat'; else engyfname = 'energy-avg-y.mat'; end

for i = 1:size(list,2)
    cd([dir list{i}]);
    
    if ~exist('energy-avg-x.mat','file') || redo_ex ==  1 && ~kencode && energy_x
        roms_energy('ocean_his.nc',[],4,1);
    end
        
    if ~kencode
        if ~exist('energy-avg-y.mat','file') || redo_ey == 1
            roms_energy('ocean_his.nc',[],4,2);
        end
    else
        if ~exist('energy-ken.mat','file') || redo_ey == 1
            [tim,eke] = engtsflc('ocean_his.nc',1,0.25,19,0,29.5,5);
        end
    end
    
    if ~exist('length_scales.mat','file') || redo_l == 1
        roms_length_scales('ocean_his.nc','u',[],volume);
    end
    
    warning off;
    grid = roms_get_grid('ocean_his.nc','ocean_his.nc',0,1);
    warning on
    
    xl = size(grid.x_rho,2);
    yl = size(grid.x_rho,1);
    zl = size(grid.z_r,1);
    
    dv(i) = max(grid.x_rho(:)-grid.x_rho(1))*max(grid.y_rho(:)-grid.y_rho(1))*max(abs(grid.z_r(:)))/(xl*yl*zl);
    run = list{i};
    label{i} = [ run(4:5) ' : ' num2str(xl-2) ' x ' num2str(yl-2) ' x ' num2str(zl)];
    
    clear EKE MKE OKE PE A time Lx Ly Lz L time_A time_E time_L
    
    load('length_scales.mat');
    LX(i) = Lx;
    LY(i) = Ly;
    LZ(i) = Lz;
   
    T{i} = time_L/86400;
    
    if ~kencode & energy_x
        load('energy-avg-x.mat');
        dPEx(i) = PE(end) - PE(1);
        intEKEx(i) = trapz(t_en,EKE)./t_en(end);
        %EKEx{i} = EKE;
        Ax{i} = A;
        Axmax(i) = max(A(2:end));
        exmax(i) = max(EKE(:));
        lxmax(i,:) =  mean(L(:,find_approx(time_L,time_A(find(A == max(A))),2))');

        figure(1); 
        plot(t_en/86400,EKE,'b');

        figure(2);
        plot(time_A/86400,A,'b');
    end

    load(engyfname);
    figure(1); 
    plot(t_en/86400,EKE,'r');
    if zl == 40, plot(t_en/86400,EKE,'g'); end
  
    figure(2);
    plot(time_A/86400,A*86400,'r');
    
    tlim = find_approx(t_en,18*86400,1); % days    
    intEKEy(i) = trapz(t_en(1:tlim),EKE(1:tlim));%./t_en(end);
    dPEy(i) = PE(tlim)-PE(1);
    %EKEy{i} = EKE;
    Ay{i} = A;
    Aymax(i) = max(A(2:end));
    eymax(i) = max(EKE(:));
    [peaks,ind] = findpeaks(A,'THRESHOLD',1e-8,'NPEAKS',2);
    tmaxgr(i,:) = time_A(ind(time_A(ind)/86400 < 7));
    lymax(i,:,1) =  mean(L(:,find_approx(time_L, tmaxgr(i,1),2))'); %  first peak - SI
    lymax(i,:,2) =  mean(L(:,find_approx(time_L, tmaxgr(i,2),2))'); % second peak - BI
   
end

figure(1)
legend('EKE-x','EKE-y');
figure(2)
legend('A-x','A-y');

%%

figure(3)
clf;
px = dv;
subplot(162);
plot(lymax(:,1,1)./max(lymax(:,1,1))*100,px,'k.','MarkerSize',14);
hold on
plot(lymax(:,2,1)./max(lymax(:,2,1))*100,px,'g.','MarkerSize',14);
l1 = xlim; l1(1) = 0; xlim(l1);
lim = ylim;
xlabel('SI Length Scales (%)');
ylabel('Grid Cell Volume (m^3)');
set(gca,'YTick',sort(unique(px)))
set(gca,'LineWidth',1.5)
box off
legend('L_x','L_y','Location','Best');
linex(100)
set(gca,'XMinorTick','off')
set(gca,'YMinorTick','off')

subplot(163)
plot(lymax(:,1,2)./max(lymax(:,1,2))*100,px,'k.','MarkerSize',14);
hold on
plot(lymax(:,2,2)./max(lymax(:,2,2))*100,px,'g.','MarkerSize',14);
xlabel('BI Length Scales (%)');
l1 = xlim; l1(1) = 0; xlim(l1);
ylim(lim);
set(gca,'YTick',sort(unique(px)))
set(gca,'LineWidth',1.5)
box off
legend('L_x','L_y','Location','Best');
linex(100)
set(gca,'XMinorTick','off')
set(gca,'YMinorTick','off')

subplot(164) 
if ~kencode & energy_x, plot(exmax./max(exmax)*100,px,'b.','MarkerSize',14); end
hold on
plot(eymax./max(eymax)*100,px,'r.','MarkerSize',14);
ylim(lim);
xlabel('Max. EKE/Horizontal Area (%)');
set(gca,'YTick',sort(unique(px)))
set(gca,'LineWidth',1.5)
legend('x avg','y avg','Location','Best');
l1 = xlim; l1(1) = 0; xlim(l1);
box off
linex(100)
set(gca,'XMinorTick','off')
set(gca,'YMinorTick','off')

subplot(165)
if ~kencode & energy_x, plot(abs(intEKEx)./max(intEKEx)*100,px,'b.','MarkerSize',14); end
hold on
plot(abs(intEKEy)./max(intEKEy)*100,px,'r.','MarkerSize',14);
xlabel('\int EKE dt (%)');
ylim(lim);
set(gca,'YTick',sort(unique(px)))
set(gca,'LineWidth',1.5)
legend('x avg','y avg','Location','Best');
l1 = xlim; l1(1) = 0; xlim(l1);
box off
linex(100)
set(gca,'XMinorTick','off')
set(gca,'YMinorTick','off')

subplot(166) 
if ~kencode & energy_x, plot(Axmax./max(Axmax)*100,px,'b.','MarkerSize',14); end
hold on
plot(Aymax./max(Aymax)*100,px,'r.','MarkerSize',14);
xlabel('Max. Growth Rate (%)');
l1 = xlim; l1(1) = 0; xlim(l1);
set(gca,'YTick',sort(unique(px)))
set(gca,'LineWidth',1.5)
box off
linex(100)
set(gca,'XMinorTick','off')
set(gca,'YMinorTick','off')

subplot(161)
ylim(lim);
for i=1:size(list,2)
   text(0.5,px(i),label{i},'VerticalAlignment','middle', ...
        'HorizontalAlignment','center','FontSize',10);
end
box off
set(gca,'Visible','off')

%cd /home/poison/deepak/ROMS/runs/si/analysis/

% plotyy = abs(bsxfun(@rdivide,ploty,max(abs(ploty))));
% figure
% hold on
% for j = 1:size(vars,1)
%     plot(plotx(:,j),plotyy(:,j),'*','Color',color{j})
%     text(plotx(:,j),plotyy(:,j), label, 'VerticalAlignment','bottom','HorizontalAlignment','right');
% end
% legend(cellstr(vars));
% xlabel('Grid Cell Volume');
% ylabel('Fraction of maximum');
% title(['Convergence : t = ' num2str(tstep) ' days']);
% info'

%% old check fields convergence code

% file  = 'ocean_his.nc';
% %dir   = '/home/poison/deepak/ROMS/runs/si/si_full/3d/';
% dir = 'E:\Work\instability\ROMS\si_full\3d\'
% list  = {'run44'}; %{'run41','run42','run43','run44','run45'};
% vars  = {'energy'; 'length_scales'; 'growth_rate'};
% color = {'r','g','b','k','c','m','y','k'};
% offset = 20;
% tstep = 120; % in days. What makes most sense?
% xloc  = 5000;
% yloc  = 30000;
% zloc  = 50;
% clear label;
% 
% k=0;
% cpb = progressbar();
% 
% for i = 1:size(list,2)
%         fname    = strcat(dir,strcat(list{i},'/',file));
%     for j = 1:size(vars,1)
%         % get grid for particular variable        
%         [x,y,z]  = roms_var_grid(fname,vars{j});
%         dx       = max(max(diff(x(:,1))));
%         dy       = max(max(diff(y(1,:))));
%         dz       = max(diff(z(:,1,1)));
%         
%         nx = size(x,1);
%         ny = size(y,2);
%         nz = size(z,1);
%         
%         % figure out timestep
%         dt = double(ncread(fname,'dt')); % in seconds
%         if isempty(findstr(file,'_avg'))
%             nfile = ncread(fname,'nHIS');
%         else
%             nfile = ncread(fname,'nAVG');
%         end
%         
%         xloc = mean(x(:));
%         yloc = mean(y(:));
%         zloc = mean(z(:));
%         
%         read_tstep = (double(tstep*86400/(dt*(nfile))));
%         read_xloc  = find_approx(x(:,1),xloc,1);
%         read_yloc  = find_approx(y(1,:),yloc,1);
%         read_zloc  = find_approx(z(:,1,1),zloc,1);
% 
%         time = ncread(fname,'ocean_time');
%         if length(time) < read_tstep, continue; end
%         
%         ploty(i,j) = double(ncread(fname,vars{j},[read_xloc read_yloc read_zloc read_tstep],[1 1 1 1]));
%         plotx(i,j) = double(dx*dy*dz/1000);
%         label(i,:) = num2str(i+offset);
%         
%         info(i) = cellstr([label(i,:) ': ' list{i} '  :  ' num2str(nx) ...
%                             ' x ' num2str(ny) ' x ' num2str(nz) ':' num2str(dt)]);
%         
%         k=k+1;
%         txt = sprintf(' Progress: i=%d, j=%d',i,j);
%         progressbarupdate(cpb,k/(size(list,2)*size(vars,1))*100,txt);
%     end
% end
% cpb.stop();
% 
% %%
% 
% plotyy = abs(bsxfun(@rdivide,ploty,max(abs(ploty))));
% figure
% hold on
% for j = 1:size(vars,1)
%     plot(plotx(:,j),plotyy(:,j),'*','Color',color{j})
%     text(plotx(:,j),plotyy(:,j), label, 'VerticalAlignment','bottom','HorizontalAlignment','right');
% end
% legend(cellstr(vars));
% xlabel('Grid Cell Volume');
% ylabel('Fraction of maximum');
% title(['Convergence : t = ' num2str(tstep) ' days']);
% info'

