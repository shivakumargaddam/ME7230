clear all
close all
clc
newsim = true;
diary temp.log

%% PROGRAM TITLE
disp('==================================================================');
disp('       Phase field method for quasi-static brittle fracture       ');
disp('         Applied Finite Element Methods(ME7230) - Project         ');
disp('                  Shiva Kumar Gaddam - MM22D014                   ');
disp('==================================================================');

%% Input
filename = 'results\tripatch';
resultfile = strcat(filename,'_results.mat');
logfile = strcat(filename,'.log');
gmshfile = 'mesh_files\patchtri.m';
eType = "tri3";                                                            % tir3 or quad4
% constants ---------------------------------------------------------------
YoungsMod = 210e9;
RatPoisson = 0.3;
EnergyGc = 5000;
Charlen = 0.1e-3;
% simulation parameters ---------------------------------------------------
duincMAX = 1e-6;
duincMIN = 1e-6;
duinc = 1e-6;
ntsteps = 100;
tolerance = 0.0000001;
plotmesh = false;

%% Main program

% ---------------------------------------------------------- Importing mesh
disp("GMESH IMPORT:")
[nodes,elements,noN,noE,uvfixed,uvmoving,phinode] = PreProc.gmsh_import(eType,gmshfile,filename,plotmesh);
disp(' ');
% ----------------------------------------------------- boundary conditions
uebcnodes = [uvfixed*2-1;uvfixed*2;uvmoving*2-1;uvmoving*2];
uebcvalues = [0; 0; 0; 0; 0; 0; 1; 1];
reacnodes = find(~uebcvalues);
phiebcnodes = [];
phiebcvalues = [];

% ---------------------------------------------------------- Save Variables
if newsim
    save("tempresults.mat","nodes","elements","noN","noE","uebcvalues","uebcnodes","reacnodes","uvmoving","uvfixed","phiebcnodes","phiebcvalues","phinode",...
        "eType","Charlen","EnergyGc","RatPoisson","YoungsMod","gmshfile","ntsteps","tolerance");
end
% ------------------------------------------------------ Staggered Approach
disp("Simulation begins! ")
Solver.stagPFM(duincMAX,duincMIN,duinc);
disp('Simulation is completed! ');

%% Post processing
load('tempresults.mat')

% phi vs axial strain plot ------------------------------------------------
phivstrain = figure('Name','phivstrain','NumberTitle','off');
epsilony = Disp/1e-3;
phfield = phievol(1,:)';
plot(epsilony(1:3:end),phfield(1:3:end),'bs','LineWidth',1.2,'MarkerSize',8,'DisplayName','MATLAB')
hold on
E = YoungsMod;
nu = RatPoisson;
stressana = (E*(1-nu)/(1+nu)/(1-2*nu).*epsilony.^2)./(EnergyGc/Charlen + E*(1-nu)/(1+nu)/(1-2*nu).*epsilony.^2);
plot(epsilony,stressana,'k','LineWidth',1.2,'DisplayName','analytical')

set(gca,'FontSize',12)
box on
pbaspect([1 1 1]);
xlabel('Axial strain [\epsilon_y]')
ylabel('Phase-field parameter [\phi]')
legend('Location','southeast')
% set(gca,'fontname','calibri')
set(gcf,'units','pixels','position',[400 300 500 500]);
title(" Benchmark problem: quadrilateral element")
subtitle('phase-field parameter as a function of axial strain')

% axial stress vs axial strain plot ---------------------------------------
stressvstrain = figure('Name','stressvstrain','NumberTitle','off');
epsilony = Disp/1e-3;
stress = zeros(size(epsilony));
for i = 1:length(stress)
    dummy = PostProc.getStress(strain(:,1,i),phfield(i),E,nu);
    stress(i) = dummy(2);
end
plot(epsilony(1:end),stress(1:end),'bs','LineWidth',1.2,'MarkerSize',8,'DisplayName','MATLAB')
hold on
E = YoungsMod;
nu = RatPoisson;
stressana = (1-phfield).^2*E*(1-nu)/(1+nu)/(1-2*nu).*epsilony;
plot(epsilony,stressana,'k','LineWidth',1.2,'DisplayName','analytical')

set(gca,'FontSize',12)
box on
pbaspect([1 1 1]);
xlabel('Axial strain [\epsilon_y]')
ylabel('Axial stress [\sigma_y]')
legend('Location','southeast')
% set(gca,'fontname','calibri')
set(gcf,'units','pixels','position',[1000 300 500 500]);
title(" Benchmark problem: quadrilateral element")
subtitle('axial stress as a function of applied axial strain')

%% End of the program
diary off
movefile('temp.log', logfile)
movefile('tempresults.mat', resultfile)

