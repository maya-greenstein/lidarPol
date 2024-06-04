% Written By:  Robert Stillwell & Rory Barton-Grimley
% Written For: ARSENL
% Updated: 2021-10-15 Kevin Sacca
%
% This functions defines a structure of function handles which are used to
% define the canonical Mueller matrices which are used in the SVLE.m
% function to calculate the complete polarization resolved measurement.
% Note here that hte definition of depolairzation d is different than
% depolarization delta per Matt Hayman's thesis.  
 
function [Matrix] = FindMuellerMatrices
%
% Inputs: none
%
% Outputs: Matrix: A structure containing function handles defining the
%                  necessary Mueller matrix types
%

%% Depolarizer 
Matrix.Depol = @(d) [1, 0 , 0 ,  0   ;
                     0,1-d, 0 ,  0   ;
                     0, 0 ,d-1,  0   ;
                     0, 0 , 0 ,2*(d-1)];
%% Uniform Depolarizer 
Matrix.UniformDepol = @(d) [1, 0 , 0 ,  0   ;
                     0, d, 0 ,  0   ;
                     0, 0 , -d,  0   ;
                     0, 0 , 0 , -d];
%% Depolarizer 
Matrix.NonUniformDepol = @(d1,d2,d3) [1, 0 , 0 ,  0   ;
                                     0,d1, 0 ,  0   ;
                                    0, 0 , -d2,  0   ;
                                    0, 0 , 0 , -d3];
%% Kattwar Reflection  
Matrix.KatRefl = @(a,e,g) [a + e, a - e, 0,  0   ;
                        a - e, a + e, 0,  0   ;
                        0,       0,   g,  0   ;
                        0,       0,   0,  g];
%% Kattwar Transmission  a -> w
Matrix.KatTransAW = @(a,e,g,k) k.*[a + e, a - e, 0,  0   ;
                                 a - e, a + e, 0,  0   ;
                                 0,       0,   g,  0   ;
                                 0,       0,   0,  g]  ;
%% Kattwar Transmission  w -> a 
Matrix.KatTransWA = @(a,e,g,k) k.*[a + e, a - e, 0,  0   ;
                                   a - e, a + e, 0,  0   ;
                                     0,       0,   g,  0   ;
                                     0,       0,   0,  g]  ;
%% Half Wave Plate (vertical)
Matrix.HalfWave = [1,0, 0, 0;
                   0,1, 0, 0;
                   0,0,-1, 0;
                   0,0, 0,-1];
%% Linear Polarizer
Matrix.LinPol = .5*[1,1,0,0;
                    1,1,0,0;
                    0,0,0,0;
                    0,0,0,0];
%% Rotation matrix in radians
Matrix.Rot   = @(T) [1,   0    ,    0    ,0;
                     0,cos(2*T),-sin(2*T),0;
                     0,sin(2*T), cos(2*T),0;
                     0,   0    ,    0    ,1];
%% Rotation matrix in degrees
Matrix.RotD  = @(T) [1,    0    ,    0     ,0;
                     0,cosd(2*T),-sind(2*T),0;
                     0,sind(2*T), cosd(2*T),0;
                     0,    0    ,    0     ,1];

%% Transmission matrix                 
Matrix.Trans = @(T) [T,0,0,0;
                     0,T,0,0;
                     0,0,T,0;
                     0,0,0,T];
%% Rotating a matrix in radians                 
Matrix.MatRot  = @(M,T) Matrix.Rot(T)*M*Matrix.Rot(-T);
%% Rotating a matrix in degrees                 
Matrix.MatRotD = @(M,T) Matrix.RotD(T)*M*Matrix.RotD(-T);
%% Functional half wave plate
Matrix.HWP = @(Trans,Theta) Matrix.Trans(Trans)*Matrix.MatRotD(Matrix.HalfWave,Theta);
%% Functional polarizer
Matrix.Pol = @(Trans,Theta) Matrix.Trans(Trans)*Matrix.MatRotD(Matrix.LinPol,Theta);
%% Transmitted Vector
Matrix.Vec   = @(T) Matrix.RotD(T)*[1;1;0;0];
%% Background vector
Matrix.Back  = [1;0;0;0];
end

