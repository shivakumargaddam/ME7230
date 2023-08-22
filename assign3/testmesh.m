clear
clc

a = 0.002;
b = 0.001;
P = -1e6;
E = 1e6;
nu = 0.3;

mesh = importdata("mesh1m2d_n88_e142.xlsx");
nodes = mesh.data.nodes;
elements = mesh.data.elements;
length(nodes)
length(elements)


if true
    meshfig = figure('Name','Mesh','NumberTitle','off');
    grid off;
    axis on;
    xlabel ('Length [m]');
    ylabel ('Width [m]');
    % xlim([0 5]);
    % ylim([0 3]);
    daspect([1 1 1]);
    set(gca, 'color', 'none');                                              % To remove background
    set(gcf,'units','pixels','position',[500 300 800 600]);                 % To change the size of the figure
    patch('Faces',elements,'Vertices',nodes,'FaceColor','cyan','FaceAlpha',1);
    set(gca,'FontSize',12)
    box on
end