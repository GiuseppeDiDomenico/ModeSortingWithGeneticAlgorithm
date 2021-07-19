% This is the main implementation of the optimization method used
% in the following paper:
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

clear all; close all;  clc
global k0 beta deltaz devStepz devStepx Z0 InputFieldG InputFieldH InputFieldT
global TargetFieldG TargetFieldH TargetFieldT Propagator mult loss

%%% PARAMETERS
n0=1.008856;                                    % Plasmonic refractide index for air silver interface
lambda=1.064e-6/n0;                             % wavelength of the plasmonic field
k0=2*pi/lambda;                                 % vacuum wavevector
w0=8e-6;                                        % input beam width
beta=1.049321-n0;                               % difference of index of refraction (for the case of 50nm of PMMA)

w0C=1.6e-6;                                     % Target Field width
FieldDisplacement=(12e-06);                     % Target Field Displacement

%%% COMPUATIONAL PARAMETERS AND GRID SET-UP
deltax=200e-9;                                  % resolution of the computational domain in the x direction (transvers)
deltaz=200e-9;                                  % resolution of the computational domain in the z direction (propagation)
mult=round((400e-9)/deltax);                    % number of grid point for every device pixel 
devSizex=36e-6;                                 % device size (transvers)
devSizez=72e-6;                                 % device size (propagation)
Lx=64e-6;                                       % computational width of system
Lz=devSizez;                                    % computational length of system

nstepx=round(Lx/deltax);                        % number of transverse grid points
nstepz=round(Lz/deltaz);                        % number of longitudinal propagation steps
devStepx=round(devSizex/deltax);                % device grid size
devStepz=round(devSizez/deltaz);                % device grid size
z=0:deltaz:Lz;                                  % z array
x=(-Lx/2+deltax/2+((1:nstepx)-1)*deltax);       % x array

%%% FIELDS DEFINITION
InputFieldG=hermiteH(0,sqrt(2)*x./w0).*exp(-(x/w0).^2);                    % inizial Hermite Gauss Field zero order (gaussian)
InputFieldH=hermiteH(1,sqrt(2)*x./w0).*exp(-(x/w0).^2);                    % inizial Hermite Gauss Field first order
InputFieldT=hermiteH(2,sqrt(2)*x./w0).*exp(-(x/w0).^2);                    % inizial Hermite Gauss Field second order
InputFieldG=InputFieldG*10/norm(InputFieldG);                              % fields normalization
InputFieldH=InputFieldH*10/norm(InputFieldH);                              % 
InputFieldT=InputFieldT*10/norm(InputFieldT);                              % 

loss=[erf(.05*(0:59)) ones(1,length(x)-120) erf(.05*(59:-1:0))];           % this adds an absorbtion boundary condition

TargetFieldG=exp(-(x/w0C).^2);                                             % target field definitions
TargetFieldG=abs(TargetFieldG*10/norm(TargetFieldG)).^2;                   % target field normalization
TargetFieldH=circshift(TargetFieldG,-round(FieldDisplacement/deltax));     % target field displacement for Gaussian field
TargetFieldT=circshift(TargetFieldG,round(FieldDisplacement/deltax));      % target field displacement for Hermite Gauss field

%%% DISPERSIVE STEP SETUP
deltaf=1/Lx; kx=zeros(1,nstepx); ntx=0;                                    % blackboxfunction needs the propagator vector to work. It compute the
for nn=1:nstepx                                                            % diffraction for a given spacial frequency step (deltaf)
    ikx=nn-ntx-1; kx(nn)=2*pi*deltaf*ikx;                                  % 
    if(nn-nstepx/2==0), ntx=nstepx; end                                    % 
end                                                                        % 
Propagator=exp(1i*(((k0*n0)^(-1)*kx.^2)*(deltaz./4.)));                    % 

Z0=zeros(devStepz,round((nstepx-devStepx)/2));                             % used in the BlackBoxFunction

%% Start L2-NSGA2
NObjective=3;                                                              % number of objectives, it rapresents the number of modes to sort
popsize=1000;                                                              % inizial population size
nVariables=devStepx*devStepz/(mult^2);                                     % number of variables to be optimized i.e. number of PMMA pixels in the device 
fun=@BlackBoxfunctionThreeBeams;                                           % function that perform the beam propagation method simulation
iterations=2005;                                                           % number of generation of the genetic algorithm
mutationProb=0.005;                                                        % probability of random mutations
nFOS=120;                                                                 % Number of Linkage structures (FOSs) to consider for the crossover
nFrequency=20;                                                             % Frequency (number of generations) with which the linkage model should be re-computed
tic

%%%% Main optimization algorithm
[xBest, fval] = L2NSGA(NObjective, popsize, nVariables, fun, iterations, mutationProb, nFrequency, nFOS, []);

evaluationTime=toc;
disp(['Evaluation time is ',num2str(evaluationTime),'  seconds'])
plot3(fval(:,1), fval(:,2), fval(:, 3), 'ro');  title('Fitness values'), xlabel('Hermite-Gauss0'), ylabel('Hermite-Gauss1'), zlabel('Hermite-Gauss2'), grid   %Pareto front plot

distance=zeros(length(fval),1);                                            % amongh the final devices we choice the one with the more ballanced performance between the beams
for gg=1:length(fval)  
    distance(gg)=norm(cross([100 100 100]-[0 0 0],fval(gg,:)-[0 0 0]))/norm([100 100 100]-[0 0 0]);  
end
[M, ind]=min(distance);
device=kron(reshape(xBest(ind,:),devStepz/mult,devStepx/mult),ones(mult)); % device reshaped and increased in resolution 

ExactTime=datestr(now,'mm-dd-yy HH-MM-SS');
if ~exist('Results', 'dir'),   mkdir('Results'), end                       % create the 'Results' folder if it does not exists

save(sprintf('Results/Workspace_%s.mat', ExactTime))                       % save all the variables
save(sprintf('Results/device_%s.mat', ExactTime),'device')                 % save the best device

%% Display the results

dn=[Z0, beta*device, Z0];
[nstepzc,nstepxc]=size(dn);

% Simulations to produce the figures
FinalFieldG=zeros(nstepzc,nstepx); Field=InputFieldG;
for hh=1:nstepzc
    Field=ifft(fft(Field).*Propagator);
    Field=Field.*exp(-1i*k0*deltaz*dn(hh,:));
    Field=ifft(fft(Field).*Propagator);
    FinalFieldG(hh,:)=Field;
end

FinalFieldH=zeros(nstepzc,nstepx); Field=InputFieldH;
for hh=1:nstepzc
    Field=ifft(fft(Field).*Propagator);
    Field=Field.*exp(-1i*k0*deltaz*dn(hh,:));
    Field=ifft(fft(Field).*Propagator);
    FinalFieldH(hh,:)=Field;
end

FinalFieldT=zeros(nstepzc,nstepx); Field=InputFieldT;
for hh=1:nstepzc
    Field=ifft(fft(Field).*Propagator);
    Field=Field.*exp(-1i*k0*deltaz*dn(hh,:));
    Field=ifft(fft(Field).*Propagator);
    FinalFieldT(hh,:)=Field;
end

figure(2)
subplot(1,2,1), imagesc([abs(FinalFieldG.').^2; abs(FinalFieldH.').^2; abs(FinalFieldT.').^2]), colorbar, axis equal
subplot(1,2,2), imagesc(device'), colorbar, axis equal
F1 = getframe(gcf);

figure(3)
subplot(3,1,1), plot(abs(FinalFieldG(end,:)).^2), hold on
plot(TargetFieldG), hold off, title('Hermite-Gauss 0')
subplot(3,1,2), plot(abs(FinalFieldH(end,:)).^2), hold on
plot(TargetFieldH), hold off, title('Hermite-Gauss 1')
subplot(3,1,3), plot(abs(FinalFieldT(end,:)).^2), hold on
plot(TargetFieldT), hold off, title('Hermite-Gauss 2')
F2 = getframe(gcf);

imwrite(F1.cdata, sprintf('Results/deviceImm1_%s.jpg', ExactTime))
imwrite(F2.cdata, sprintf('Results/deviceImm2_%s.jpg', ExactTime))
