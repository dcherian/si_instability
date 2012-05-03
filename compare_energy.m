%%
load energy-avg-y.mat
figure(1)
plot(time_A./86400,A./max(A),'b');hold on
figure(2)
plot(t_en./86400,EKE./max(EKE),'b'); hold on

load energy-ken.mat
figure(1)
plot(time_A./86400,A./max(A),'r');
legend('me','ken')
figure(2)
plot(t_en./86400,EKE./max(EKE),'r');
linex([t_en(12) t_en(24) t_en(end)]./86400,' ');
legend('me','ken')
