% This function is part of the optimization method used in the following 
% paper:
%
% G. Di Domenico, D. Weisman, A. Panichella, D. Rroitman and A. Arie,
% "Planar on-Chip Mode Sorter" (2021)
%
% ------------------------------- Copyright -------------------------------
% Copyright (C) 2021 Annibale Panichella and Giuseppe Di Domenico. You are 
% free to use this code for research purposes. All publications which use
% this code should acknowledge the use and add the referece to the paper
% listed above.
%
% This project is free software: you can redistribute it and/or modify it
% under the terms of the GNU Lesser General Public License as published by
% the Free Software Fundation, either version 3.0 of the Licence, or (at 
% your opinion) any later version.
% -------------------------------------------------------------------------

function Fits = BlackBoxfunctionThreeBeams(chromosome)
% This function encapsulate the Beam propagation method simulation
% takes in imput the variables to be sotimized in the form of a vector 
% (chromosome) that encodes one device to be simulated and output three 
% number rapresenting how well the device encode the desired function for 
% the different input fields.

global k0 beta deltaz devStepz devStepx Z0 InputFieldG InputFieldH InputFieldT
global TargetFieldG TargetFieldH TargetFieldT Propagator mult loss

device = kron(reshape(chromosome,devStepz/mult,devStepx/mult),ones(mult));  % the chromosome is reshaped and the resolution is increased for the simulation
dn=[Z0, beta*device, Z0];                                                   % sets the refractive index for the full simulation domain 

Field=InputFieldG;                                                          % simulation for the gaussian field
for hh=1:devStepz
    Field=ifft(fft(Field).*Propagator).*loss;                             %%% computation of the first half step diffraction of the absorption boundary
    Field=Field.*exp(-1i*k0*deltaz*dn(hh,:));                             %%% computation of the phase due to an inomogeneus refractive index
    Field=ifft(fft(Field).*Propagator);                                   %%% computation of the second half step diffraction
end
FitsG=norm(TargetFieldG-abs(Field).^2);                                     % measure of the difference between the output field and a target field defined by the user

Field=InputFieldH;                                                          % simulation for the Hermite Gauss of the first and second orfer fields folowing the same steps
for hh=1:devStepz
    Field=ifft(fft(Field).*Propagator).*loss;
    Field=Field.*exp(-1i*k0*deltaz*dn(hh,:));
    Field=ifft(fft(Field).*Propagator);
end
FitsH=norm(TargetFieldH-abs(Field).^2);

Field=InputFieldT;
for hh=1:devStepz
    Field=ifft(fft(Field).*Propagator).*loss;
    Field=Field.*exp(-1i*k0*deltaz*dn(hh,:));
    Field=ifft(fft(Field).*Propagator);
end
FitsT=norm(TargetFieldT-abs(Field).^2);

Fits=[FitsG,FitsH,FitsT];                                                   % results values to be minimize
end
