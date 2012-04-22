% Generate parameters for oceanic case similar to Jones and Thorpe (1992)

f  = 1e-4;
N2 = 1E-5;
N  = sqrt(N2);

% From paper
Vzjt = 5E-3;
Ljt  = 132e3; % wavelength 132 km
Xjt  = 2*Ljt;
Yjt  = 960E3;
Hjt  = 5040;
Ujt  = 0.5*Vzjt*Hjt; % max. velocity (0 at center)

Ldjt = N*Hjt/f; % deformation radius
Rojt = Ujt/(f*Ljt);
Rijt = N2/Vzjt^2;
Bjt  = Ldjt/Ljt; % Burger Number

% for ocean
Ho  = 100; 
Ldo = N*Ho/f
Lo  = Ldo/Bjt
Uo  = 2*Rojt*f*Lo
Vzo = Uo/Ho

Xo = 2*Lo;
Yo = Yjt/Ldjt * Lo;

delta1 = 410/Yjt * Yo;
delta2 = 550/Yjt * Yo;

fprintf('\n Domain: X = %f , Y = %f, Z = %f, pertub: %f to %f \n\n', Xo,Yo,Ho, delta1, delta2);

