classdef PreProc
    % HELP - PreProc (Mesh Import)
    %   -- Methods:
    %               gmsh_import [obtain elements and nodes from gmsh]
    
    properties
    end
    
    methods (Static)
        function [eType,nodes,elements,noN,noE,uvfixed,uvmoving,phinode,crack] = gmsh_import(gmshfile,filename,plotmesh)
            % HELP - gmsh_import
            %   -- Nodes, Elements of the domain and physical groups are imported from 'Gmsh'
            %   -- input file is xxxx.m
            %   -- [nodes,elements,uvfixed,uvmoving,phinode] = gmsh_import(eType,gmshfile,filename,plotmesh);
            
            arguments
                gmshfile {mustBeFile}
                filename {mustBeText}
                plotmesh = true
            end

            % Mesh Import -------------------------------------------------
            run(gmshfile);
            nodes3D = msh.POS;
            if isfield(msh,'TRIANGLES')
                eType = "tri3";
                eldetails = msh.TRIANGLES ;
            elseif isfield(msh,'QUADS')
                eType = "quad4";
                eldetails = msh.QUADS;
            end
            bnodes = msh.LINES;

            uvfixed = bnodes(bnodes(:,3)==1,1:2);
            uvmoving = bnodes(bnodes(:,3)==2,1:2);
            phinode = bnodes(bnodes(:,3)==3,1:2);
            crack = bnodes(bnodes(:,3)==4,1:2);
            uvfixed = unique(uvfixed);
            uvmoving = unique(uvmoving);
            phinode = unique(phinode);
            crack = unique(crack);

            nodes = nodes3D(:,1:end-1);                                     % since it is 2D
            elements = eldetails(:,1:end-1);                                % the last column is physical group ID -- Grain ID

            noN = size(nodes,1);                                            % No of nodes
            noE = size(elements,1);                                         % No of elements

            % Mesh Visualization ------------------------------------------
            fprintf("\tNo of 'NODES' in the generated mesh\t\t: %d\n", noN);
            fprintf("\tNo of 'ELEMENTS' in the generated mesh\t: %d\n", noE);
            fprintf("\tElement type: %s\n", eType);
            

            if plotmesh
                fmesh = figure('Name','Mesh Visualization','NumberTitle','off');
                patch('Faces',elements,'Vertices',nodes,'FaceColor','cyan','FaceAlpha',1);
                hold on
                plot(nodes(uvfixed,1),nodes(uvfixed,2),'.r','MarkerSize',15);
                plot(nodes(uvmoving,1),nodes(uvmoving,2),'.b','MarkerSize',15);
                plot(nodes(phinode,1),nodes(phinode,2),'.k','MarkerSize',15);
                grid off;
                axis on;
                box on
                xlabel x-axis;
                ylabel y-axis;
                daspect([1 1 1]);
                set(gca, 'color', 'none');                                              % To remove background
                set(gcf,'units','pixels','position',[100 100 800 600]);                 % To change the size of the figure
                colormap lines
                meshfile = strcat(filename,'_mesh');
                savefig(meshfile)
                % meshfile = strcat(filename,'_mesh.png');
                % exportgraphics(fmesh,meshfile,'Resolution',500);
            end
        end

        function [cracknodes,crackvalues] = crackgen(nodes,crack,crackwidth)
            noN = size(nodes,1);
            nocracknodes = length(crack);
            distance = crackwidth/2;
            cracknodes = [];
            crackvalues = [];
            for inode = 1:noN
                for icracknode = 1:nocracknodes
                    if norm(nodes(inode,:)-nodes(crack(icracknode),:)) <= distance
                        cracknodes = [cracknodes; inode];
                        crackvalues = [crackvalues; 1-norm(nodes(inode,:)-nodes(crack(icracknode),:))/distance];
                    end
                end
            end
        end
        
    end
end

