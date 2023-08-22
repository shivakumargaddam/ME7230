clear all
close all
clc

c = 340;
meshes = ["mesh1m2d_n88_e142.xlsx" "mesh1m4d_n304_e542.xlsx" "mesh1m8d_n1145_e2160.xlsx" "mesh1m16d_n4508_e8758.xlsx"];
h = [1/2 1/4 1/8 1/16];
Length = 5;
Width = 3;
whwa50 = [];

%% analytical frequencies
natfana = [];
m = 0:100;
n = 0:100;
mncombo = table2array(combinations(m,n));
for icombo = 1:length(mncombo)
    im = mncombo(icombo,1);
    in = mncombo(icombo,2);
    natfana = [natfana;c/2*sqrt((im/Length)^2+(in/Width)^2)];
end
[natfana,I] = sort(natfana);
mn50 = mncombo(I(1:50),:);

%%
for imesh = 1:length(h)
    mesh = importdata(meshes(imesh));
    nodes = mesh.data.nodes;
    elements = mesh.data.elements;
    [natf_fem] = acoustic (nodes,elements,c,true);
    natf_fem = natf_fem/2/pi;           % in hertz
    whwa50 = [whwa50;natf_fem(50)/natfana(50)];
end
whwa = natf_fem(1:50)./natfana(1:50);
%% plot
figure;
plot(h,whwa50,'.b-','MarkerSize',19,'LineWidth',2,'DisplayName','\omega^h/\omega^a')
yline(1.05,'-','LineWidth',2,'DisplayName','\omega^h/\omega^a=1.05');
set(gca,'FontSize',12)
box on
pbaspect([1 1 1]);
xlabel('h [m]')
ylabel('\omega^h/\omega^a')
% set(gca, 'fontname','calibri');                       % To remove background and change font
set(gcf,'units','pixels','position',[400 300 500 500]);
legend('Location','northwest')
title("Acoustic Problem")
subtitle('mesh convergence: \omega^h/\omega^a<1.05')

figure;
plot(1:50,whwa,'.b-','MarkerSize',19,'LineWidth',2,'DisplayName','\omega^h/\omega^a')
set(gca,'FontSize',12)
xlim([1,50])
box on
pbaspect([1 1 1]);
xlabel('mode number')
ylabel('\omega^h/\omega^a')
% set(gca, 'fontname','calibri');                       % To remove background and change font
set(gcf,'units','pixels','position',[950 300 500 500]);
title("Acoustic Problem")
subtitle('\omega^h/\omega^a vs. mode number')



