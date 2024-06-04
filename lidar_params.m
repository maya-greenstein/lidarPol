% function to call all of the _____ lidar parameters
function [lp] = lidar_params
    
    lp.pE = 5e-6;           % [J] pulse energy
    lp.rr = 19000;          % [Hz] repitition rate
    lp.wl = 532e-9;         % [m] laser wavelength in a vacuum
    lp.D = .05;             % [m] telescope diameter
    lp.A = (pi*(lp.D/2)^2); % [m] telescope area 
    lp.BinSize = 0.0005;    % [m] range bin size
    
%     lp.alt = 30;            % [m] altitude of lidar
%     lp.eta = 2.06e-4;       % full system efficiency
%     lp.eta = .05;       % full system efficiency
    lp.rxEff = .0094;
    lp.txEff = .022;
    
    % transmitter
    lp.Lambda                   = 532.00e-9;    % laser wavelength               [m]
    lp.LaserW0                  = 0.005;        % laser output beam radius       [m]
    lp.Divergence               = 12e-3;        % beam divergence              [rad]
    
    % receiver 
    lp.FNumTele                 = 13.9;         % F number of the telescope      [ ]
    lp.FieldStop                = 20e-3;         % Diameter of the field stop     [m]
    lp.Separation               = 75e-3;        % laser-telescope separation     [m]
    lp.BandPassFilt             = 0.0234;       % pass band of the filter used  [nm]
    lp.TransPol                 = 2.5;          % Rx/Tx offset angle           [deg]
    lp.FOV                      = (lp.FieldStop)/(lp.FNumTele*lp.D); %8e-3;
%     
    lp.Pointing = .001;         
    lp.WTransAngle = asin(sin(lp.Pointing)/1.33);
    
%     lp.t = .9999;              % lower atm transmission at 589nm
    lp.t = 1;  
    
    lp.CoOffset = 0;        % [deg] Rx/Tx offset angle for co-pol
    lp.CxOffset = 90;       % [deg] Rx/Tx offset angle for cx-pol
end