%% Function to Compute the Stokes Vectors for a Dual-Polarized System

function [Srxco, Srxcx] = SVLE(alts, DSurf, Dair, NL, dR, surfReflec, Trans, betaRL, G, Matrix, lp, wps) 
    % Inputs: 
        % alts:        array of altitudes the signal is propogating through
        % DSurf:       depolarization value of the surface of reflection
        % DSurf:       depolarization value of the air
        % NL:          number of transmitted photons
        % dR:          integration range bin 
        % surfReflec:  array of amount of surface reflectance at altitude*
        %              *0s except for at surface
        % Trans:       amount of transmission at altitudes 
        % betaRL:      rayleigh backscatter coefficient at altitudes 
        % Matrix:      mueller matrix descriptions
        % lp:          parameters of the lidar system being used 
        % wps:         water surface parameters: empty if not water surface
        %              [alpha, eta, gamma] 

    nBins = length(alts);
    
    % Transmitter Path                    (Mtx)
    Mtx   = Matrix.Trans(lp.txEff);
    
    % Receiver path                        (Mrx)
    MrxCo = Matrix.Trans(lp.rxEff)*Matrix.Pol(1,lp.CoOffset);
    MrxCx = Matrix.Trans(lp.rxEff)*Matrix.Pol(1,lp.CxOffset);
    
    % preallocate 
    Srxco = zeros(4, nBins);
    Srxcx = zeros(4, nBins);

    for m = 1:nBins

        % Atmospheric Scattering Matrix: Depolarization/Rayleigh beta
        F = Matrix.Depol(Dair).*betaRL(m);

        % First atmospheric transmission       (Tatm1)
        Tatm1 = Matrix.Trans(Trans(m));

        % Second atmospheric transmission      (Tamt2)
        Tatm2 = Matrix.Trans(Trans(m));

        % Transmitted Vector                   (Stx)
        Stx   = Matrix.Vec(lp.CoOffset);

        

        ScaleAtmos = lp.A./((alts(m)-alts(1))^2)*dR*NL*G(m); % can add more shots & overlap function here 
        ScaleGround = lp.A./((alts(m)-alts(1))^2)*dR*NL*G(m);
%         ScaleAtmos = lp.A./((alts(m)-alts(1))^2)*dR; % can add more shots & overlap function here 
%         ScaleGround = lp.A./((alts(m)-alts(1))^2)*dR;

        if isempty(wps)
            % Surface Depolarizing Lambertian Scattering
            Rg = Matrix.Depol(DSurf).*surfReflec(m);
            
            Srxco(:,m) = (ScaleAtmos*MrxCo*Tatm2*F*Tatm1*Mtx*Stx)+...
                (ScaleGround*MrxCo*Tatm2*Rg*Tatm1*Mtx*Stx);
        %     Srxco(:,m) = MrxCo*(ScaleAtmos*Tatm2*F*Tatm1*Mtx*Stx)+...
        %         MrxCo*(ScaleGround*Tatm2*Rg*Tatm1*Mtx*Stx);
            Srxcx(:,m) = (ScaleAtmos*MrxCx*Tatm2*F*Tatm1*Mtx*Stx)+...
                (ScaleGround*MrxCx*Tatm2*Rg*Tatm1*Mtx*Stx);
        else 
            
            Rw = -Matrix.KatRefl(wps(1),wps(3),wps(3)).* ...
                Matrix.Depol(DSurf).*surfReflec(m);
            
            Srxco(:,m) = (ScaleAtmos*MrxCo*Tatm2*F*Tatm1*Mtx*Stx)+...
                (ScaleGround*MrxCo*Tatm2*Rw*Tatm1*Mtx*Stx);
        %     Srxco(:,m) = MrxCo*(ScaleAtmos*Tatm2*F*Tatm1*Mtx*Stx)+...
        %         MrxCo*(ScaleGround*Tatm2*Rg*Tatm1*Mtx*Stx);
            Srxcx(:,m) = (ScaleAtmos*MrxCx*Tatm2*F*Tatm1*Mtx*Stx)+...
                (ScaleGround*MrxCx*Tatm2*Rw*Tatm1*Mtx*Stx);
            
        end

    end

end