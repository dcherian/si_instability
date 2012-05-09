% Compare scalings for Taylor & Ferrari (2009) setup

f0 = 1e-4;
M2 = -4.2e-7;
N2 = 9e-6;
nu = 1e-4;
H = 50; % Mixed layer depth!!!
v0 = abs(M2/f0)*H;
Ri0 = 0.5;
m = 0.27;2*pi/H;

% wavelengths are problematic
lraf = 2*pi/sqrt(f0/nu * sqrt(1/Ri0 - 1)*(1+ N2^2/M2^2)^(-1)) % raf
lst = 2*pi/(pi/sqrt(1-Ri0)*f0/v0) % stone
lst15 = 2*pi/(15*f0/v0) % stone
ldc = abs(2*pi/(2*pi/H * M2/N2))

% growth rates agree
sigmast = sqrt(1/Ri0-1)*f0
sigmaraf = sqrt(M2^2/N2 - f0^2)

% onset of KH instability
u0 = 9.4e-6;
tkh = 1/sigmaraf * (log(abs(v0/u0)) + log(Ri0/m/H * (sqrt(3) - sqrt(1/Ri0 - 1)))) / 86400*24