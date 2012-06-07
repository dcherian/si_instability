function eig4bvp

sigma = 0.05;
z = linspace(0,1,100);
zsolve = [z(1) z(2) z(2:end)];
solinit = bvpinit(zsolve,@mat4init,sigma);
sol = bvp4c(@stoneode,@stonebc,solinit);

fprintf('Eigenvalue: %7.3f.\n',...
        sol.parameters)

Sxint = deval(sol,zsolve);

phi = [0 Sxint(1,:)]
plot(z,phi)
title(['\phi | \sigma = ' num2str(sol.parameters)]) 
xlabel('z')
% ------------------------------------------------------------
function dydz = stoneode(z,y,sigma,k,l,Ri)
        k=0.6;
    l=8;
    Ri=1;
h = sigma + k*z;
dydz = [  y(2)
         -((k^2*(Ri*h.^4 + (3-Ri)*h.^2 -2) + l^2 * (Ri*h.^4 + (1-Ri)*h.^2))./(h.^2 .* (1-h.^2).^2)).*y(1) ];
% ------------------------------------------------------------
function res = stonebc(yleft,yright,lambda)
res = [  yleft(1,1) 
         yright(1,1)-1
         yleft(1,2)-1
         yright(1,2)
         yleft(2,1)-1 ];
% ------------------------------------------------------------
function yinit = mat4init(z,k,l,Ri)
yinit = [  cos(z)
          -sin(z) ];