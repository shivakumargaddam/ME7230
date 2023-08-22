% input filenames ---------------------------------------------------------
filename = 'results\shearbox_h0o05_20d';
gmshfile = 'mesh_files\shearbox_h0o05_20d.m';
% input constants ---------------------------------------------------------
YoungsMod = 208.2e9;
RatPoisson = 0.3;
EnergyGc = 2700;
Charlen = 0.022e-3;
% input simulation parameters ---------------------------------------------
duincMAX = 1.1e-7;
duincMIN = 1e-8;
duinc = 1.1e-7;
ntsteps = 400;
tolerance = 0.0001;
plotmesh = true;
mode = 2;                                                                   % Tension:1 , Shear:2
 