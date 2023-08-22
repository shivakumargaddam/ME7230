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
input_template

%% Main program
% ---------------------------------------------------------- Importing mesh
disp("GMESH IMPORT:")
[eType,nodes,elements,noN,noE,uvfixed,uvmoving,phinode,crack] = PreProc.gmsh_import(gmshfile,filename,plotmesh);
disp(' ');

% ---------------------------------------------------------------- Precrack
[cracknodes,crackvalues] = PreProc.crackgen(nodes,crack,5e-5);

% ----------------------------------------------------- boundary conditions
if mode == 1
    disp("Mode I : Tension")
    uebcnodes = [uvfixed*2-1;uvfixed*2;uvmoving*2-1;uvmoving*2];
    uebcvalues = [zeros(length(uvfixed),1);zeros(length(uvfixed),1);zeros(length(uvmoving),1);ones(length(uvmoving),1)];
elseif mode == 2
    disp("Mode II : Shear")
    uebcnodes = [uvfixed*2-1;uvfixed*2;uvmoving*2-1;uvmoving*2];
    uebcvalues = [zeros(length(uvfixed),1);zeros(length(uvfixed),1);ones(length(uvmoving),1);zeros(length(uvmoving),1)];
end

% ---------------------------------------------------------- reaction nodes
reacnodes = find(~uebcvalues);
phiebcnodes = [];
phiebcvalues = [];
disp("Press any key to start simulation!")
pause

% ---------------------------------------------------------- Save Variables
if newsim
    save("tempresults.mat","nodes","elements","noN","noE","uebcvalues","uebcnodes","reacnodes","uvmoving","uvfixed","phiebcnodes","phiebcvalues","phinode",...
        "eType","Charlen","EnergyGc","RatPoisson","YoungsMod","gmshfile","ntsteps","tolerance","filename","crackvalues","cracknodes","duinc","duincMIN","duincMAX");
end

% ------------------------------------------------------ Staggered Approach
disp("Simulation begins! ")
Solver.stagPFM(duincMAX,duincMIN,duinc);
disp('Simulation is completed! ');

%% Post processing
PostProc.visdisp
PostProc.visphi

%% End of the program
diary off
resultfile = strcat(filename,'_results.mat');
logfile = strcat(filename,'.log');
movefile('temp.log', logfile)
movefile('tempresults.mat', resultfile)

