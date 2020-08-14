% M = load_mesh_freesurfer(datapath_freesurfer,id,hemi)
%   Loading meshes from a freesurfer file directory
%     datapath_freesurfer - path where data files are
%     id - subject id where to get the mesh
%     hemi - 'lh' or 'rh' for left/right hemisphere (default: 'lh')
%
% Herve Lombaert, May 2015
%
function M = load_mesh_freesurfer(datapath_freesurfer,id,hemi)

    if nargin<3;  hemi = 'lh';  end;

    %% Read directly in freesurfer hierarchy
    polar           = '.polar.mat';
    eccen           = '.eccen.mat';

    addpath([get_functiondir('load_mesh_freesurfer') '/freesurfer']);

    % Mesh Basic
    infile = [datapath_freesurfer '/' id '/surf/' hemi '.white'];
    if exist(infile,'file');
        [M.X,M.F]   = read_surf(infile);  M.F = M.F+1; % FS are 0-indexed
        M.n = size(M.X,1);
    end
    
    % Mesh Sulcal Depth
    infile = [datapath_freesurfer '/' id '/surf/' hemi '.sulc'];
    if exist(infile,'file');
        M.C = read_curv(infile);
    end

    % Mesh Parcellation
    infile = [datapath_freesurfer '/' id '/label/' hemi '.aparc.annot'];
%     infile = [datapath_freesurfer '/' id '/label/' hemi '.aparc.a2005s.annot'];
    %infile = [datapath_freesurfer '/' id '/label/' hemi '.aparc.DKTatlas40.annot'];
    if exist(infile,'file');
        M.P  = read_parcellation(infile);
    end
    % Mesh Cortical Thickness
    infile = [datapath_freesurfer '/' id '/surf/' hemi '.thickness'];
    if exist(infile,'file');
        M.T = read_curv(infile);
    end
    
    
    % Mesh Retinotopy    
    infile = [datapath_freesurfer '/' id '/surf/' hemi '.polar.mat'];
    if exist(infile,'file')
        tmp = load(infile);
        if isfield(tmp,'S1')
            M.R                   = tmp.S1;
        elseif isfield(tmp,'topo_new2')
            M.R(:,1)              = tmp.topo_new2(:,1);
            M.R_significance(:,1) = tmp.topo_new2(:,2);
        end
    end
    infile = [datapath_freesurfer '/' id '/surf/' hemi '.eccen.mat'];
    if exist(infile,'file')
        tmp = load(infile);
        if isfield(tmp,'S1')
            M.R(:,2)              = tmp.S1;
        elseif isfield(tmp,'topo_new2')
            M.R(:,2)              = tmp.topo_new2(:,1);
            M.R_significance(:,2) = tmp.topo_new2(:,2);
        end
    end

    % Mesh Retinotopy Mask
    infile = [datapath_freesurfer '/' id '/surf/' id '_' hemi '_mask.mat'];
    if exist(infile,'file')
        tmp = load(infile);
        %M.Mask = tmp.full_data;
        th = 7;
        if max(M.R_significance(:)) < 10;  th = .27;  end;  % such as subject marc
        M.Mask = 1*(M.R_significance>th);
        M.R              = bsxfun(@times, M.R,              M.Mask);
        M.R_significance = bsxfun(@times, M.R_significance, M.Mask);
        M.R             (isnan(M.R))              = 0;
        M.R_significance(isnan(M.R_significance)) = 0;
    end

    % Check if there are visual areas
    infile = [datapath_freesurfer '/' id '/surf/' 'allvisual-' hemi '.1D.dset'];
    if exist(infile,'file')
        fid = fopen(infile,'r');
        if 1
            C = textscan(fid,' %d %d', 'CommentStyle','#');
            C = [C{1}+1 C{2}];
        else
            i = 1;
            C = [];
            while 1
                line = fgets(fid);
                if feof(fid);       break;     end;
                if line(1) == '#';  continue;  end;
                c = sscanf(line,'%f %f');
                C(i,:) = [c(1)+1 c(2)];
                i = i+1;
            end
        end
        fclose(fid);

        M.VA = zeros(M.n,1);
        M.VA(C(:,1)) = C(:,2);

        M.VA(M.VA==50) = 7;
        M.VA(M.VA==51) = 8;
        M.VA(M.VA==52) = 9;
        M.VA(M.VA==53) = 9;
        M.VA(M.VA==100) = 10;
        M.VA(M.VA==101) = 11;
        M.VA(M.VA==102) = 12;
        M.VA(M.VA==103) = 13;
        M.VA(M.VA==150) = 14;
        M.VA(M.VA==151) = 15;
        M.VA(M.VA==152) = 16;
        M.VA(M.VA==153) = 17;
        M.VA(M.VA==154) = 18;
    end

end
