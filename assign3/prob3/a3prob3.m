clear all
close all
clc

E = 1e6;
P = -1e6;
nu = [0.3];
meshes = ["mesh1mm10d_n1045_e320.xlsx"];
h = [0.1*1e-3];

for inu = nu
    L2norm = figure('Name','L2norm','NumberTitle','off');
    hold on
    H1norm = figure('Name','H1norm','NumberTitle','off');
    hold on
    for imethod = 1:1
        L2 = zeros(size(meshes));
        H1 = zeros(size(meshes));
        for imesh = 1
            mesh = importdata(meshes(imesh));
            nodes = mesh.data.nodes;
            elements = mesh.data.elements;
            ubcnode = find(ismember(nodes(:,1),0));
            vbcnode = find(ismember(nodes(:,2),0));
            R = nodes(:,1).^2+nodes(:,2).^2;
            fbcnode = find(round(sqrt(R),9) == 0.001);
            [~,bnodes] = sortrows(nodes(fbcnode,:));
            belements = [bnodes(1:2:end-2) bnodes(2:2:end-1) bnodes(3:2:end)];
            if true
                % Plot mesh
                meshfig = figure('Name','Mesh','NumberTitle','off');
                grid off;
                axis on;
                xlabel ('x-axis [mm]');
                ylabel ('y-axis [mm]');
                xlim([-0.00005 0.00205]*1e3);
                ylim([-0.00005 0.00205]*1e3);
                daspect([1 1 1]);
                set(gca, 'color', 'none');                                              % To remove background
                set(gcf,'units','pixels','position',[500 300 800 600]);                 % To change the size of the figure
                patch('Faces',elements(:,[1 5 2 6 3 7 4 8]),'Vertices',nodes*1e3,'FaceColor','cyan','FaceAlpha',1);
                hold on
                plot(nodes(vbcnode,1)*1e3,nodes(vbcnode,2)*1e3,'*b');
                plot(nodes(ubcnode,1)*1e3,nodes(ubcnode,2)*1e3,'*r');
                patch('Faces',belements,'Vertices',nodes*1e3,'EdgeColor','yellow','FaceColor','none','LineWidth',2);
                set(gca,'FontSize',12)
                box on
            end
            [L2(imesh),H1(imesh)] = cylipipe (nodes,elements,ubcnode,vbcnode,belements,inu,E,P,imethod);
        end
        % L2 norm
        figure(L2norm)
        loglog(h,L2,'*-')
        polyfit(log(h),log(L2),1)
        % H1 norm
        figure(H1norm)
        loglog(h,H1,'*-')
        polyfit(log(h),log(H1),1)
    end
    %% L2 norm
    figure(L2norm)
    set(gca,'FontSize',12)
    box on
    pbaspect([1 1 1]);
    xlabel('log(h)')
    ylabel('log(L_2)')
    legend({'normal','selective integration technique','u-p formulation'},'Location','northwest')
    set(gcf,'units','pixels','position',[500 300 800 600]);
    title("Cylindrical Pipe: convergence study - L_2 norm")

    %% H1 norm
    figure(H1norm)
    set(gca,'FontSize',12)
    box on
    pbaspect([1 1 1]);
    xlabel('log(h)')
    ylabel('log(H_1)')
    legend({'normal','selective integration technique','u-p formulation'},'Location','northwest')
    set(gcf,'units','pixels','position',[500 300 800 600]);
    title("Cylindrical Pipe: convergence study - H_1 semi-norm")
end