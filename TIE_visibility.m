% Original paper: https://doi.org/10.1364/OE.16.003223

% Author: Kiarash Tajbakhsh, 
% Swiss Federal Laboratories for Materials Science, Empa, Überlandstrasse 120, 8600 Dübendorf, Switzerland.
% University of Fribourg, Faculty of Science and Medicine, Chemin du Musée 8, 1700 Fribourg, Switzerland.

% Please also cite this paper upon using the code:
% A survey of micro-CT for 3D virtual histology of FFPE blocks,
% Kiarash Tajbakhsh et. al. 2024
%% Input
close all; clear all;
% All lenght units are in um unless stated otherwise
pixel_size = 127; % Physical pixel size of the detector
min_SOD = 30; % in mm due to geometrical contrains.
obj_thickness = 2.5e3; %Thickness of the object of interest.
I0 = 1; % Incident wave amplitude
lambda = 8.2e-5; k=2*pi/lambda; % Mean energy wavelenght (lambda), and wavenumber (k)

SDD = linspace(60,1100,1000); %in um; 
SOD = SDD'*logspace(-1.57,0,1000);
SDD_2D = SDD'*ones(1,size(SDD,2));
M = SDD_2D./SOD;

R2 = (SDD_2D-SOD)*1e3; %in um
x = (-100:100); %detector span in um

sig_src = 5; sig_det = 2*pixel_size;    % change these according to your scanner
sig_sys = sqrt((M-1).^2.*(sig_src./M).^2+(sig_det./M).^2); % imaging system FWHM of PSF
sig_obj = 1;
sigma = sqrt(sig_sys.^2+sig_obj.^2);
delta = 1e-6; %real part of the refractive index.
phi = -k*(delta)*obj_thickness; % imposed phase-shift by the specimen.

opt_SOD = SDD.*sqrt(sig_src^2/(sig_det^2+sig_src^2)); % optimal SOD for...
% your scanner and given SDD in (mm)

%% Body
constant = (R2)./M./(k*sigma.^3*sqrt(2*pi))*abs(phi);
s = zeros(size(SOD,1),size(SOD,2),size(x,2));
V = zeros(size(SOD,1),size(SOD,2));
SNR = zeros(size(SOD,1),size(SOD,2));
Nf = zeros(size(SOD,1),size(SOD,2));
for i = 1:size(SOD,1) % i goes as SDD
    for j = 1:size(SOD,2) % j goes as SOD (fraction of SDD)
        s(i,j,:) = I0*(1-constant(i,j)*x.*exp(-(0.5/sigma(i,j)^2).*x.^2))/M(i,j)^2; % The modulated intensity distribution
        Nf(i,j) = sigma(i,j)^2*k*M(i,j)/R2(i,j); % Fresnel number of the setup
        V(i,j) = 0.2420*((lambda*R2(i,j))./(sigma(i,j)^2*2*pi*M(i,j)))*abs(phi); % visibility
    end
end

%% Visualization

t = rot90(V,-1);
t = flip(t,2);

[C,h]=contour(SDD,pixel_size./M(500,:),t);
set(gca,'YScale', 'log')
ylabel('Voxel size (\mum)','FontSize',18); xlabel('SDD (mm)','FontSize',18); 
set(gca, 'FontSize', 18);

clabel(C,h,"FontSize",18,"fontname","times")
% highlighting valid geometries for FFPE blocks.
xx = [60:1100 1100 60]; yy = [(pixel_size*min_SOD)./(60:1100) pixel_size pixel_size]; 
patch(xx,yy,'k', 'FaceAlpha', 0.05,'EdgeColor','none')

grid on