%%
t_end = 6*86400;

load energy-avg-y.mat
figure(1)
tind = find_approx(time_A,t_end,1);
plot(time_A(1:tind)./86400,A(1:tind)*86400,'b');hold on
figure(2)
tind = find_approx(t_en,t_end,1);
plot(t_en(1:tind)./86400,EKE(1:tind)./max(EKE(1:tind)),'b'); hold on

load energy-ken.mat
figure(1)
tind = find_approx(time_A,t_end,1);
plot(time_A(1:tind)./86400,A(1:tind),'r');
legend('me','ken')
figure(2)
tind = find_approx(t_en,t_end,1);
plot(t_en(1:tind)./86400,EKE(1:tind)./max(EKE(1:tind)),'r');
linex([t_en(12) t_en(24) t_en(end)]./86400,' ');
legend('me','ken')
