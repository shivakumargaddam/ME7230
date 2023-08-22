classdef PostProc
    % HELP - PreProc (Mesh Import)
    %   -- Methods:
    %               getStress   [calculates 2d stress at a point]
    %               visdisp     [visualize the dispalcement at final step and generate animated gif]
    %               visphi      [visualize the damage at final stepa and generate animated gif]

    properties
    end

    methods (Static)
        function [stress] = getStress(strain,phi,E,nu)
            arguments
                strain (3,1) {mustBeNumeric, mustBeFinite}
                phi    (1,1) {mustBeNumeric, mustBeFinite}
                E      (1,1) {mustBeNumeric, mustBeFinite}
                nu     (1,1) {mustBeNumeric, mustBeFinite}
            end
            krs = 1e-9;                                                     % Small Residual Stiffness
            Cpsn = (E/(1+nu)/(1-2*nu))*[1-nu nu 0;
                nu 1-nu 0;
                0 0 (1-2*nu)/2];                                            % Plane strain stiffness tensor
            stress = ((1-phi)^2+krs)*Cpsn*strain;
        end

        function visdisp(finalstate,animation,simdatafile)
            arguments
                finalstate = true
                animation = true
                simdatafile {mustBeFile} = 'tempresults.mat'
            end
            % Extraction of Variables -------------------------------------
            simdata         = matfile(simdatafile,'Writable',true);
            nodes           = simdata.nodes;
            elements        = simdata.elements;
            uvdispevol      = simdata.uvdispevol;
            ntsteps         = simdata.ntsteps;
            filename        = simdata.filename;
            % Declaration of Variables ------------------------------------
            noN = size(nodes,1);                                            % No of nodes
            noE = size(elements,1);                                         % number of elements
            dispmagx = 1;                                                   % Displacement multiplier - x
            dispmagy = 1;                                                   % Displacement multiplier - y

            if finalstate
                % ---------------------------------------------------------- Final Image
                fuvdisp = figure('Name','Displacement','NumberTitle','off');
                displacednodes(:,1) = nodes(:,1) + dispmagx*uvdispevol(1:2:2*noN,ntsteps+1);  % displaced nodal coordinates: x-coordinate
                displacednodes(:,2) = nodes(:,2) + dispmagy*uvdispevol(2:2:2*noN,ntsteps+1);  % displaced nodal coordinates: y-coordinate
                patch('Faces',elements,'Vertices',displacednodes,'FaceColor','cyan','FaceAlpha',1);

                grid off;
                box on;
                axis off;
                xlabel x-axis/Length;
                ylabel y-axis/Height;
                xlim([0-0.05 1+0.05]*1e-3)
                ylim([0-0.05 1+0.05]*1e-3)
                daspect([1 1 1]);
                set(gca, 'color', 'none');                                  % To remove background
                set(gcf,'units','pixels','position',[300 300 800 600]);     % To change the size of the figure
                dispfile = strcat(filename,'_disp');
                savefig(dispfile)
            end

            if animation
                % ------------------------------------------------------------------ GIF
                fuvdisp_gif = figure('Name','Displacement_gif','NumberTitle','off');
                grid off;
                box on;
                axis off;
                xlabel x-axis/Length;
                ylabel y-axis/Height;
                xlim([0-0.05 1+0.05]*1e-3)
                ylim([0-0.05 1+0.05]*1e-3)
                axis off
                daspect([1 1 1]);
                set(gca, 'color', 'none');                                  % To remove background
                set(gcf,'units','pixels','position',[300 300 800 600]);     % To change the size of the figure

                animfile = strcat(filename,'_disp.gif');
                for tstep = 1:ntsteps+1
                    cla;
                    animatednodes(:,1) = nodes(:,1) + dispmagx*uvdispevol(1:2:2*noN,tstep);    % displaced nodal coordinates: x-coordinate
                    animatednodes(:,2) = nodes(:,2) + dispmagy*uvdispevol(2:2:2*noN,tstep);    % displaced nodal coordinates: y-coordinate
                    patch('Faces',elements,'Vertices',animatednodes,'FaceColor','cyan','FaceAlpha',1);
                    % Capture the plot as an image
                    frame = getframe(fuvdisp_gif);
                    im = frame2im(frame);
                    [imind,cm] = rgb2ind(im,256);
                    % Write to the GIF File
                    if tstep == 1
                        imwrite(imind,cm,animfile,'gif', 'Loopcount',inf,'DelayTime',0.08);
                    else
                        imwrite(imind,cm,animfile,'gif','WriteMode','append','DelayTime',0.08);
                    end
                end
            end
        end

        function visphi(finalstate,animation,simdatafile)
            arguments
                finalstate = true
                animation = true
                simdatafile {mustBeFile} = 'tempresults.mat'
            end
            % Extraction of Variables -------------------------------------
            simdata         = matfile(simdatafile,'Writable',true);
            nodes           = simdata.nodes;
            elements        = simdata.elements;
            phievol         = simdata.phievol;
            ntsteps         = simdata.ntsteps;
            filename        = simdata.filename;
            % Declaration of Variables ------------------------------------
            noN = size(nodes,1);                                            % No of nodes
            noE = size(elements,1);                                         % number of elements

            if finalstate
                % ---------------------------------------------------------- Final Image
                fphi = figure('Name','phase-field','NumberTitle','off');
                phifinal = phievol(:,ntsteps+1);
                patch('Faces',elements,'Vertices',nodes,'FaceVertexCData',phifinal,'FaceColor','interp','EdgeColor','none');                     % Rigid
                colorbar;
                clim([0 1]);

                grid off;
                box on;
                axis off;
                xlabel x-axis/Length;
                ylabel y-axis/Height;
                daspect([1 1 1]);
                set(gca, 'color', 'none');                                  % To remove background
                set(gcf,'units','pixels','position',[300 300 800 600]);     % To change the size of the figure
                colormap turbo;
                dispfile = strcat(filename,'_phi');
                savefig(dispfile)
            end

            if animation
                % ------------------------------------------------------------------ GIF
                fphi_gif = figure('Name','phase-field_gif','NumberTitle','off');
                grid off;
                box on;
                axis off;
                xlabel x-axis/Length;
                ylabel y-axis/Height;
                daspect([1 1 1]);
                set(gca, 'color', 'none');                                  % To remove background
                set(gcf,'units','pixels','position',[300 300 800 600]);     % To change the size of the figure
                colormap turbo;

                animfile = strcat(filename,'_phi.gif');
                for tstep = 1:ntsteps+1
                    cla;
                    phifinal = phievol(:,tstep);
                    patch('Faces',elements,'Vertices',nodes,'FaceVertexCData',phifinal,'FaceColor','interp','EdgeColor','none');
                    colorbar;
                    clim([0 1]);
                    % Capture the plot as an image
                    frame = getframe(fphi_gif);
                    im = frame2im(frame);
                    [imind,cm] = rgb2ind(im,256);
                    % Write to the GIF File
                    if tstep == 1
                        imwrite(imind,cm,animfile,'gif', 'Loopcount',inf,'DelayTime',0.08);
                    else
                        imwrite(imind,cm,animfile,'gif','WriteMode','append','DelayTime',0.08);
                    end
                end
            end
        end
% -------------------------------------------------------------------------
    end
end

