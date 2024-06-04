% Written By:  Robert Stillwell
% Written For: ARSENL
% This function calculates the overlap function of the laser beam and the
% telescope assuming they are separated by some small distance. It assumes
% that the beam and telescope field of view expand as cones.

function [G] = Overlap(Distances, SP)
   
    % Inputs: Altitudes:  The altitudes of each lidar measurement
    %         SP:         The lidar system's operational parameters
 
    % Outputs: Overlap:   An array containing expected overlap percentage
    %                     given the system configuration of Supr
    
    Distances = -Distances + Distances(1,1);
    
    %% Temp variables 
    phi1   = SP.FOV/2;        % half of the telescope's FOV
    phi0   = SP.Divergence/2; % half of the beam divergence angle
    rsep   = SP.Separation;   % space between telescope and laser
    
    %% Determining overlap proportion
    % This is from Halldorsson and Langerholc, App. Opt. 17 2 (1978)
    G=zeros(length(Distances),1);
    
    for m=1:length(Distances)
        rfov   = SP.D/2 + phi1*Distances(m);
        rlaser = sqrt(SP.LaserW0^2 + (phi0*Distances(m))^2);

        if rsep >= rfov + rlaser % before any overlap
            Ainter = 0;

        elseif rsep <= abs(rfov-rlaser) % complete overlap with laser < fov
            Ainter = (pi*rlaser^2);

        elseif rsep <= rfov + rlaser && rsep >= abs(rfov-rlaser) %middle ground
            Ainter = rfov^2*acos((rsep^2+rfov^2-rlaser^2)/(2*rsep*rfov)) + ...
                rlaser^2*acos((rsep^2+rlaser^2-rfov^2)/(2*rsep*rlaser)) - ...
                1/2*sqrt(((rfov+rlaser)^2-rsep^2)*(rsep^2-(rfov-rlaser)^2));
        end

        G(m)=Ainter/(pi*rlaser^2);
    end
end