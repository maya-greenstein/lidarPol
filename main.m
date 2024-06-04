%% ASEN 6365 Final Project 
% Maya Greenstein 05/01/23
clear; clc; close all;

%% Define Constants 

c = 299792458; %[m/s]
kB = 1.380649e-23; %[(m^2kg)/(s^2K)]
h = 6.62607015e-34; %[m^2kg/s]
av = 6.0221408e23; %[1/mol]
ep0 = 8.854187817e-12; %[F/m]
me = 9.1093837e-31; %[kg]
e = 1.60217663e-19; %[coulumbs]
Ru = 8.314; %[kgm^2/(s^2Kmol)]

delta = .00365; % air depolarization ratio on central Cabannes line at 532nm
D = 2*delta/(1+delta); 

%% Solve Output Signal 

lp = lidar_params; % load in potassium lidar parameters 
Matrix = FindMuellerMatrices;

% Term 1: Number of laser shots out
NL = lp.pE/((h*c)/lp.wl); 

% Term 3: Telescope params 
A = (pi*(lp.D/2)^2); % primary mirror area 

%% Rayleigh Counts from Alts 

targetHeight = 10; %[m]
lidarHeight = 40; %[m]
% targetHeight = 2010; %[m]
% lidarHeight = 2040; %[m]

% Boulder lat/lon
lat = 40.016869;
lon = -105.279617;

% Today 
year = 2023;
time = 12*3600;
doy = day(datetime('now'),'dayofyear');

dR = lp.BinSize;
alts = lidarHeight:-dR:(targetHeight-3); %High to low 
nBins = length(alts);

% MSIS datas
molCols = [1, 2, 3, 4, 5, 7, 8, 9];

%Rayleigh 
[TR, rhoR] = atmosnrlmsise00(alts,lat,lon,year,doy,time, 'None');
TR = TR(:,2); %Just want atmospheric temp
nR = sum(rhoR(:,molCols), 2); 
n_VR = nR./av; %n/V in ideal gas law [mols/m^3]
PR = (n_VR.*Ru.*TR)/100; %pressure [Pa->mbar]

betaRL = (2.938e-32)*(PR./TR).*(1./((lp.wl).^(4.0117))); %Rayleigh backscatter

P = 0.7629*(1+0.9324*cosd(180)^2); %Rayleigh scattering phase function
betaTotal = (betaRL.*4*pi)./P;   
Trans = exp(-cumtrapz(-alts,betaTotal));


% NRl = NL*betaRL'*dR.*(A./((alts-targetHeight).^2)).*(lp.t^2)*lp.eta;
% NRl=flip(NRl); %The Rayleigh counts will happen right in front of the detector (why we need overlap function)
% This will look diff that for Na bc NL (from pulse energy), 
% A (reciever area) and eta are ALL much lower 

% Overlap 
[G] = Overlap(alts, lp);

%% Stokes Stuff 



% ROUGH SURFACE 

% Define ground
roughGroundPulseWidth = .1;
reflectionGround1 =  normpdf(alts, targetHeight, roughGroundPulseWidth);

refAng = pi; % assuming lambertian scattering
wGround = 0.2; %Ground albedo (reflectivity)
BRDF =  wGround/refAng;

% Add dimension to surface using normal distrobution
reflectionGround1 = BRDF*reflectionGround1./max(reflectionGround1);

% Depolarization
groundD1 = .9;

% Solve SVLE
[Sco1, Scx1] = SVLE(alts, groundD1, D, NL, dR, reflectionGround1, Trans, betaRL, G, Matrix, lp, []); 

% figure()
% title('Rough Surface')
% hold on
% % plot(Sco1(1,1000:end), alts(1000:end), 'b');
% % plot(Scx1(1,1000:end), alts(1000:end), 'r');
% plot(Sco1(1,:), alts, 'b');
% plot(Scx1(1,:), alts, 'r');
% legend('Co-Pol', 'Cx-Pol');
% xlabel('S_{Rx} Intensity');
% ylabel('Altitude [m]');



% SMOOTH SURFACE (LAMBERTIAN) 

% Define ground
smoothGroundPulseWidth = .03;
reflectionGround2 =  normpdf(alts, targetHeight, smoothGroundPulseWidth);

reflectionGround2 = wGround*reflectionGround2./max(reflectionGround2);


% Depolarization
groundD2 = .2;

[Sco2, Scx2] = SVLE(alts, groundD2, D, NL, dR, reflectionGround2, Trans, betaRL, G, Matrix, lp, []); 
% 
figure()
sgtitle('Surface Scattering Polarization')

subplot(1,2,1)
title('Smooth Surface')
hold on
plot(Sco2(1,:), alts, 'b');
plot(Scx2(1,:), alts, 'r');
xlabel('S_{Rx} Intensity');
ylabel('Altitude [m]');
ylim([8 40])

subplot(1,2,2)
title('Rough Surface')
hold on
% plot(Sco1(1,1000:end), alts(1000:end), 'b');
% plot(Scx1(1,1000:end), alts(1000:end), 'r');
plot(Sco1(1,:), alts, 'b');
plot(Scx1(1,:), alts, 'r');
legend('Co-Pol', 'Cx-Pol');
xlabel('S_{Rx} Intensity');
ylabel('Altitude [m]');
ylim([8 40])



% WATER SURFACE (FRESNEL) 

% Compute elements for water surface fresnel reflec mat
WAlpha = 0.5*(tan(lp.Pointing - lp.WTransAngle)./ ...
                   tan(lp.Pointing + lp.WTransAngle)).^2;
               
WEta = 0.5*(sin(lp.Pointing - lp.WTransAngle)./ ...
                   sin(lp.Pointing + lp.WTransAngle)).^2;

WGamma = - (tan( lp.Pointing - lp.WTransAngle).* ...
                           sin(lp.Pointing - lp.WTransAngle))./ ...
                  (tan( lp.Pointing + lp.WTransAngle).* ...
                             sin(lp.Pointing + lp.WTransAngle));

oceanPulseSize = .08;                         
reflectionWater = normpdf(alts, targetHeight, oceanPulseSize);
reflectionWater = reflectionWater./ max(reflectionWater);
                         
                         
groundD3 = .0392; %delta =.02

[Sco3, Scx3] = SVLE(alts, groundD3, D, NL, dR, reflectionWater, Trans, betaRL, G, Matrix, lp, [WAlpha, WEta, WGamma]); 

figure()
title('Water Surface')
hold on
% plot(Sco3(1,1000:end), alts(1000:end), 'b');
% plot(Scx3(1,1000:end), alts(1000:end), 'r');
plot(Sco3(1,:), alts, 'b');
plot(Scx3(1,:), alts, 'r');
legend('Co-Pol', 'Cx-Pol');
xlabel('S_{Rx} Intensity');
ylabel('Altitude [m]');


