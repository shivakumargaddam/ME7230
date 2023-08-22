clear all
close all
clc

data = readmatrix('data.csv');
% Problem1 - strain energy vs dof
ql1 = data(:,[1 2]);
qq1 = data(:,[3 4]);
tl1 = data(:,[5 6]);
tq1 = data(:,[7 8]);
% Problem 2 - max. stress vs hole radius
ql2 = data(:,[9 10]);
tl2 = data(:,[11 12]);
% Problem 3 - strain energy vs nu - plane stress
ql3 = data(:,[13 14]);
qq3 = data(:,[15 16]);
tl3 = data(:,[17 18]);
tq3 = data(:,[19 20]);
% Problem 3 - strain energy vs nu - plane strain
ql31 = data(:,[21 22]);
qq31 = data(:,[23 24]);
tl31 = data(:,[25 26]);
tq31 = data(:,[27 28]);

figure
set(gca,'FontSize',12)
hold on
plot(2*ql1(:,1),ql1(:,2),'.-','MarkerSize',17,'LineWidth',1.5,'DisplayName','quadrilateral-linear')
plot(2*qq1(:,1),qq1(:,2),'.-','MarkerSize',17,'LineWidth',1.5,'DisplayName','quadrilateral-quadratic')
plot(2*tl1(:,1),tl1(:,2),'.-','MarkerSize',17,'LineWidth',1.5,'DisplayName','triangle-linear')
plot(2*tq1(:,1),tq1(:,2),'.-','MarkerSize',17,'LineWidth',1.5,'DisplayName','triangle-quadratic')
plot([0,1e7]',[1987.2,1987.2]','r--','DisplayName','convergence limit')
% ylim([1975 1988.2])
xlim([0 250000])
xlabel('degrees of freedom')
ylabel('total strain energy [J]')
legend('Location','southeast')
set(gcf,'units','pixels','position',[100 100 600 600]);
box on
title("Strain Energy vs 'dof' ")
% pbaspect([1 1 1])

figure
set(gca,'FontSize',12)
hold on
plot(0.5*ql2(:,2),ql2(:,1)/1e6,'.-','MarkerSize',17,'LineWidth',1.5,'DisplayName','quadrilateral-linear')
plot(0.5*tl2(:,2),tl2(:,1)/1e6,'.-','MarkerSize',17,'LineWidth',1.5,'DisplayName','triangle-linear')
xlabel('hole radius [m]')
ylabel('max stress [MPa]')
legend('Location','southeast')
set(gcf,'units','pixels','position',[200 100 600 600]);
box on
title("Max. stress vs hole radius")
% pbaspect([1 1 1])

figure
set(gca,'FontSize',12)
hold on
plot(ql3(:,1),ql3(:,2),'.-','MarkerSize',17,'LineWidth',1.5,'DisplayName','quadrilateral-linear')
plot(qq3(:,1),qq3(:,2),'.-','MarkerSize',17,'LineWidth',1.5,'DisplayName','quadrilateral-quadratic')
plot(tl3(:,1),tl3(:,2),'.-','MarkerSize',17,'LineWidth',1.5,'DisplayName','triangle-linear')
plot(tq3(:,1),tq3(:,2),'.-','MarkerSize',17,'LineWidth',1.5,'DisplayName','triangle-quadratic')
xlabel("poisson's ratio")
ylabel('total strain energy [J]')
legend('Location','northeast')
set(gcf,'units','pixels','position',[300 100 600 600]);
box on
title("Strain Energy vs \nu - Plane stress")
% pbaspect([1 1 1])

figure
set(gca,'FontSize',12)
hold on
plot(ql31(:,1),ql31(:,2),'.-','MarkerSize',17,'LineWidth',1.5,'DisplayName','quadrilateral-linear')
plot(qq31(:,1),qq31(:,2),'.-','MarkerSize',17,'LineWidth',1.5,'DisplayName','quadrilateral-quadratic')
plot(tl31(:,1),tl31(:,2),'.-','MarkerSize',17,'LineWidth',1.5,'DisplayName','triangle-linear')
plot(tq31(:,1),tq31(:,2),'.-','MarkerSize',17,'LineWidth',1.5,'DisplayName','triangle-quadratic')
xlabel("poisson's ratio")
ylabel('total strain energy [J]')
legend('Location','northeast')
set(gcf,'units','pixels','position',[400 100 600 600]);
box on
title("Strain Energy \nu - Plane strain")
% pbaspect([1 1 1])