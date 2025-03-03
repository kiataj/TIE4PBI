%% Micro-CT Imaging Simulation for 3D Virtual Histology
% This code simulates a micro-CT setup used for 3D virtual histology of FFPE blocks.
% Please cite the following paper if you use or extend this code:
% "A survey of micro-CT for 3D virtual histology of FFPE blocks,"
% Kiarash Tajbakhsh et al., 2024. 
% DOI: https://doi.org/10.1364/OE.16.003223

%% Initialization
close all; clear; clc;

%% Define Physical and Geometric Parameters
% Units: micrometers (µm) unless stated otherwise
pixel_size   = 18;           % Detector pixel size [µm]
min_SOD      = 30;            % Minimum source-to-object distance [mm] (geometrical constraint)
obj_thickness = 2.5e3;         % Thickness of the object [µm]
I0           = 1;             % Incident wave amplitude

lambda = 8.2e-5;             % Wavelength [µm]
k = 2*pi / lambda;           % Wavenumber [1/µm]

%% Define Imaging Geometry
% SDD: Source-to-Detector Distance range (µm)
SDD = linspace(60, 1100, 1000);

% SOD: Source-to-Object Distance grid (each row corresponds to an SDD, values logarithmically spaced)
SOD = SDD' * logspace(-1.57, 0, 1000);

% Create a 2D grid for SDD (used in element-wise operations)
SDD_2D = SDD' * ones(1, numel(SDD));

% Magnification factor: M = SDD / SOD
M = SDD_2D ./ SOD;

% R2: Distance from object to detector in µm (note conversion from mm)
R2 = (SDD_2D - SOD) * 1e3;

% Detector coordinate range (µm)
x = -100:100;

%% System Resolution Parameters
% Define source and detector resolution characteristics (FWHM values)
sig_src = 5;                     % Source width [µm]
sig_det = 2 * pixel_size;        % Detector PSF width [µm]

% Compute the system PSF combining source and detector contributions, scaled by magnification
sig_sys = sqrt((M - 1).^2 .* (sig_src ./ M).^2 + (sig_det ./ M).^2);

% Intrinsic object PSF
sig_obj = 1;                     
sigma = sqrt(sig_sys.^2 + sig_obj.^2);  % Combined PSF (system and object)

% Phase shift introduced by the specimen (depends on refractive index decrement)
delta = 1e-6;                    % Refractive index decrement
phi = -k * delta * obj_thickness;  % Imposed phase-shift

% Compute optimal SOD for a given SDD and scanner parameters
opt_SOD = SDD .* sqrt(sig_src^2 / (sig_det^2 + sig_src^2));

%% Simulation: Intensity, Fresnel Number, and Visibility
% Preallocate arrays for the modulated intensity distribution (s), visibility (V),
% and Fresnel number (Nf).
s = zeros(size(SOD,1), size(SOD,2), numel(x));
V = zeros(size(SOD,1), size(SOD,2));
Nf = zeros(size(SOD,1), size(SOD,2));

% Constant term used in the intensity modulation
constant = R2 ./ M ./ (k * sigma.^3 * sqrt(2*pi)) * abs(phi);

% Loop over all SDD and SOD combinations to compute:
% - The modulated intensity distribution on the detector (s)
% - The Fresnel number (Nf)
% - The visibility metric (V)
for i = 1:size(SOD,1)        % Loop over SDD values
    for j = 1:size(SOD,2)    % Loop over SOD values
        % Calculate the modulated intensity distribution across the detector
        s(i,j,:) = I0 * (1 - constant(i,j) * x .* exp(-0.5 * (x.^2) / sigma(i,j)^2)) / M(i,j)^2;
        
        % Compute the Fresnel number of the imaging setup
        Nf(i,j) = sigma(i,j)^2 * k * M(i,j) / R2(i,j);
        
        % Calculate the visibility metric (contrast measure)
        V(i,j) = 0.2420 * ((lambda * R2(i,j)) / (sigma(i,j)^2 * 2*pi * M(i,j))) * abs(phi);
    end
end

%% Visualization: Contour Plot of Visibility
% Adjust the visibility matrix orientation for plotting
V_vis = rot90(V, -1);
V_vis = flip(V_vis, 2);

% Plot contour of visibility vs. SDD (mm) and voxel size (pixel_size/M at mid SDD)
figure;
[contourMatrix, hContour] = contour(SDD, pixel_size ./ M(500,:), V_vis);
set(gca, 'YScale', 'log');
ylabel('Voxel size (\mum)', 'FontSize', 18);
xlabel('SDD (mm)', 'FontSize', 18);
set(gca, 'FontSize', 18);
clabel(contourMatrix, hContour, "FontSize", 18, "fontname", "times");

% Highlight valid geometries for FFPE blocks with a transparent patch
xx = [60:1100, 1100, 60];
yy = [(pixel_size * min_SOD) ./ (60:1100), pixel_size, pixel_size];
patch(xx, yy, 'k', 'FaceAlpha', 0.05, 'EdgeColor', 'none');

grid on;
