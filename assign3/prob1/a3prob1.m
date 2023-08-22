clear all
close all
clc

meshtype = [1 2]; % the (integer)values indicate the amount of distortion; 1-square and >1-trapezoid
h = 1./2.^[1:6]; % different mesh sizes
L2 = zeros(length(h),length(meshtype));
H1 = zeros(length(h),length(meshtype));

for imeshtype = 1:length(meshtype)
    disp(['mesh type: ',num2str(meshtype(imeshtype))])
    for imesh = 1:length(h)
        [nodes,elements,ubcnode] = mesh_generation(meshtype(imeshtype),h(imesh));
        [L2(imesh,imeshtype),H1(imesh,imeshtype)] = poisson(nodes,elements,ubcnode,true);
    end
end

%% L2 norm
L2norm = figure('Name','L2norm','NumberTitle','off');
loglog(h,L2(:,1),'*-b','LineWidth',1.5,'DisplayName','structured')
dummy = polyfit(log(h),log(L2(:,1)),1);
L2_square = dummy(1,1)
hold on
loglog(h,L2(:,2),'*-r','LineWidth',1.5,'DisplayName','trapezoidal')
dummy = polyfit(log(h),log(L2(:,2)),1);
L2_trpzdl = dummy(1,1)
set(gca,'FontSize',12)
box on
pbaspect([1 1 1]);
xlabel('log(h)')
ylabel('log(L_2)')
legend('Location','northwest')
% set(gca,'fontname','calibri')
set(gcf,'units','pixels','position',[300 300 500 500]);
title("Structured vs. Trapezoidal meshes")
subtitle('convergence study : L_2 norm')
%% H1 norm
H1norm = figure('Name','H1norm','NumberTitle','off');
loglog(h,H1(:,1),'*-b','LineWidth',1.5,'DisplayName','structured')
dummy = polyfit(log(h),log(H1(:,1)),1);
H1_square = dummy(1,1)
hold on
loglog(h,H1(:,2),'*-r','LineWidth',1.5,'DisplayName','trapezoidal')
dummy = polyfit(log(h),log(H1(:,2)),1);
H1_trpzdl = dummy(1,1)
set(gca,'FontSize',12)
box on
pbaspect([1 1 1]);
xlabel('log(h)')
ylabel('log(H_1)')
legend('Location','northwest')
% set(gca,'fontname','calibri')
set(gcf,'units','pixels','position',[1000 300 500 500]);
title("Structured vs. Trapezoidal meshes")
subtitle('convergence study : H_1 semi-norm')


function [nodes,elements,ubcnode] = mesh_generation(imeshtype,h)
noE = (1/h)^2;
% Defining required matrices
nodes3D = zeros(8,2,noE);    % nodal coordinates in 3D array in which each 2D array represents one element
nodes = zeros(noE*8,2);      % nodal coordinates of entire mesh
elements = zeros(noE,8);     % nodes of each element in the mesh
% the template
nodes3D(1,:,1) = [0 0];
nodes3D(2,:,1) = [0.5 0];
nodes3D(3,:,1) = [0.5 imeshtype/(imeshtype+1)];
nodes3D(4,:,1) = [0 1/(imeshtype+1)];
nodes3D(1,:,2) = nodes3D(2,:,1);
nodes3D(2,:,2) = [1 0];
nodes3D(3,:,2) = [1 1/(imeshtype+1)];
nodes3D(4,:,2) = nodes3D(3,:,1);
nodes3D(1,:,3) = nodes3D(3,:,1);
nodes3D(2,:,3) = nodes3D(3,:,2);
nodes3D(3,:,3) = [1 1];
nodes3D(4,:,3) = [0.5 1];
nodes3D(1,:,4) = nodes3D(4,:,1);
nodes3D(2,:,4) = nodes3D(3,:,1);
nodes3D(3,:,4) = nodes3D(4,:,3);
nodes3D(4,:,4) = [0 1];
if 2*h<=1
    nodes3D(:,:,1:4) = nodes3D(:,:,1:4)*2*h;
    for ih = 2:1/2/h
        nodes3D(:,:,(1:4)+4*(ih-1)) = nodes3D(:,:,1:4);
        nodes3D(:,1,(1:4)+4*(ih-1)) = nodes3D(:,1,(1:4)+4*(ih-1)) + 2*h*(ih-1);
    end
    for jh = 2:1/2/h
        nodes3D(:,:,(1:4/(2*h))+4/(2*h)*(jh-1)) = nodes3D(:,:,1:4/(2*h));
        nodes3D(:,2,(1:4/(2*h))+4/(2*h)*(jh-1)) = nodes3D(:,2,(1:4/(2*h))+4/(2*h)*(jh-1)) + 2*h*(jh-1);
    end
end
for i = 1:noE
    nodes3D(5,:,i) = (nodes3D(1,:,i)+nodes3D(2,:,i))/2;
    nodes3D(6,:,i) = (nodes3D(2,:,i)+nodes3D(3,:,i))/2;
    nodes3D(7,:,i) = (nodes3D(3,:,i)+nodes3D(4,:,i))/2;
    nodes3D(8,:,i) = (nodes3D(4,:,i)+nodes3D(1,:,i))/2;
end
% Combining nodal coordinates of each element into a 2D array
for i= 1:noE
    nodes(((i-1)*8+1):i*8,:) = nodes3D(:,:,i);
end
% Sorting and removing duplicates of the nodal coordinates
nodes = unique(nodes,'rows');
% Extracting the node numbers of each element -- ElNodes
for i = 1:noE
    [~, elements(i,:)] = ismember(nodes3D(:,:,i),nodes,'rows');
end
ubcnodeN = unique(find(ismember(nodes(:,2),1)));
ubcnodeE = unique(find(ismember(nodes(:,1),1)));
ubcnodeW = unique(find(ismember(nodes(:,1),0)));
ubcnodeS = unique(find(ismember(nodes(:,2),0)));
ubcnode = unique([ubcnodeN;ubcnodeE;ubcnodeW;ubcnodeS]);
end