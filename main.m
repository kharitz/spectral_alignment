clc
close all;
clear all;

addpath(genpath('./utils/'))

%% Initializationsof paths and parameters
% Folder paths to load and save data
data_name = 'mindboggle'; % Change the path for adni
datapath = ['./dataset/' data_name];
path_to_manlabels = 'path_to_vtk_file_for_mindboggle'; % Manaual parcels
path_to_csvfile = 'path_to_data.csv'; % Phenotype data
save_path_embed        = ['./output/' data_name '_embedding/'];
save_path_trans        = ['./output/' data_name '_tranformation/'];
save_path_forpy        = ['./output/' data_name '_forpy/'];
if ~exist(save_path_embed,'dir');  mkdir(save_path_embed);  end;
if ~exist(save_path_trans,'dir');  mkdir(save_path_trans);  end;
if ~exist(save_path_forpy,'dir');  mkdir(save_path_forpy);  end;


% Parameters 
hemi             = { 'lh', 'rh'}; % Left or Right Hemisphere
ne               = 5;     % total number of eigenvectors to decompose  (to warp basis)
hemi_type        = 1;     % Choose 1 for left hemishpere
krot             = 5;     % eigenvectors used to transformation
matching_samples = 10000;   % number of points used to find basis transformation
w_sulcal         = 1;      % sulcal weights
ki               = 5;      % eigenvectors used for matching
niter            = 2*krot; % nb iteration to refine basis transformation matrix (sign flip, reordering)
use_basis_transf = 1;      % refine basis transformation
matching_mode    = 'partial';  % 'partial', 'complete' - Matching using partial/complete set of mesh points
do_save          = 1;      % To save files and results
%% Embeddings and Transformation

% Load the reference mesh
       
id0 = 'HLN-12-6'; % Can pick any subject from the dataset
disp(['*** Loading reference' id0]);
filename = [datapath '/' id0 ];

if exist(filename,'dir')
    M0  = load_mesh_freesurfer(datapath, id0, hemi{hemi_type});
    % Embeddings
    M0.E = embedding(M0,ne,'randomwalk');
else
    disp(['*** Reference_' id0 'does not exist']);
end
if do_save;  save([save_path_embed id0 '.mat'], 'M0', '-v7.3');  end;


disp('*** Saving aligned embeddings ***');
files = dir(datapath);
all_files = files(3:end);
for iter_subject = 1:numel(all_files)
               
    % Load Mesh
    id = all_files(iter_subject).name;
    disp(['*** Loading ' id]);
    filename = [datapath '/' id ];

    if exist(filename,'dir')
        M  = load_mesh_freesurfer(datapath, id, hemi{hemi_type});
        % Embeddings
        M.E = embedding(M,ne,'randomwalk');
    else
        disp(['***' id 'does not exist']);
    end


    M = flipeigen(M0,M,3);  % Flip eigenvector sign -> should be common with M0
    if strcmp(id,id0);  continue;  end;  % skip self-registration

    % Transform Embedding - Find transformation
    Mw = M;
    disp(['*** Warping ' id ' to ' id0]);

    min_ssd  = 1e10;
    best_R12 = [];
    best_R21 = [];

    R12  = eye(krot+1);  % Initial Transformation (identity)
    R21  = eye(krot+1);  % Initial Transformation (identity)

    % Get basis
    U1 = [ones(M0.n,1)/sqrt(M0.n) M0.E.U(:,1:krot)];
    U2 = [ones(M.n, 1)/sqrt(M.n)  M. E.U(:,1:krot)];

    % Refine Transformation
    for iter_match = 1:niter

        % Eigenvectors to use in this iteration
        trunc_id = min( 5+1 + ceil((iter_match-1)*.65), krot+1 );  % which eigenvectors to keep
        trunc    = zeros(krot+1,krot+1);
        trunc(1:trunc_id,1:trunc_id) = eye(trunc_id);

        rng(1);
        ki_match = ki;  if iter_match == 1;  ki_match = 3;  end;  % start with 3 (less ambiguous);

        switch matching_mode
            case 'complete'  % Complete mesh

                % Rough Correspondence Between M0 and M
                [c1,c2] = match([w_sulcal*M0.C M0.E.X(:,1:ki_match)],  [w_sulcal*M.C Mw.E.X(:,1:ki_match)]);
                c = struct('corr12',c1, 'corr21',c2);
                if do_display;  subplot(121);  displaymesh(Ms.X,Ms.F,Mw.E.X(:,2)); view(-90,0); link_light;  end;

                % Change base
                U12 = U1(c.corr21,:);  %U12 = bsxfun(@rdivide, U12, sqrt(sum(U12.^2)));  % U1 on mesh 2
                R12 = U2' * U12;  % product U1 (on mesh2) with U2 (on mesh2)

                U21 = U2(c.corr12,:);  %U21 = bsxfun(@rdivide, U21, sqrt(sum(U21.^2)));  % U2 on mesh 1
                R21 = U1' * U21;  % product U2 (on mesh1) with U1 (on mesh1)

            case 'partial' % Partial mesh

                % Rough Correspondence Between M0 and M
                n    = min([matching_samples  M0.n  M.n]);  % nb of points to use in matching
                idx1 = randperm(M0.n);  idx1 = idx1(1:n);
                idx2 = randperm(M .n);  idx2 = idx2(1:n);

                [c1,c2] = match([w_sulcal*M0.C(idx1) M0.E.X(idx1,1:ki_match)],  [w_sulcal*Mw.C(idx2) Mw.E.X(idx2,1:ki_match)]);
                c = struct('corr12',c1, 'corr21',c2);

                % Change base
                U12 = U1(idx1(c.corr21),:);  %U12 = bsxfun(@rdivide, U12, sqrt(sum(U12.^2)));  % U1 on mesh 2
                U22 = U2(idx2,:);            %U22 = bsxfun(@rdivide, U22, sqrt(sum(U22.^2)));  % U2 on mesh 2
                R12 = U22' * U12;  % product U1 (on mesh2) with U2 (on mesh2)

                U21 = U2(idx2(c.corr12),:);  %U21 = bsxfun(@rdivide, U21, sqrt(sum(U21.^2)));  % U2 on mesh 1
                U11 = U1(idx1,:);            %U11 = bsxfun(@rdivide, U11, sqrt(sum(U11.^2)));  % U1 on mesh 1
                R21 = U11' * U21;  % product U2 (on mesh1) with U1 (on mesh1)

                % Rescale weights, based on downsampling
                R12 = R12 * sqrt(M0.n*M.n)/n;
                R21 = R21 * sqrt(M0.n*M.n)/n;

        end

        % Use first eigenvectors only - identity in higher dimension
        R12 = trunc*R12*trunc + (eye(krot+1)-trunc) * eye(krot+1) * M.n /sqrt(M.n*M0.n);
        R21 = trunc*R21*trunc + (eye(krot+1)-trunc) * eye(krot+1) * M0.n/sqrt(M.n*M0.n);

        % Make it symmetric        
        nR12 = (R12 + R21')/2;
        nR21 = (R21 + R12')/2;
        R12  = nR12;
        R21  = nR21;
        

        % Transformed basis
        X2 = U2 * R12 * diag([1; M.E.V(1:krot)].^-.5);
        Mw.E.X(:,1:krot) = X2(:,2:end);

        % Check Energy
        ki_ssd = 5;
        switch matching_mode
            case 'complete'
                ssd_fwd(iter_match) = sum(sum(([w_sulcal*M0.C                 M0.E.X(:,1:ki_ssd)             ] -  [w_sulcal*Mw.C(c.corr12)       Mw.E.X(c.corr12,1:ki_ssd)      ]).^2,2));
                ssd_bwd(iter_match) = sum(sum(([w_sulcal*M0.C(c.corr21)       M0.E.X(c.corr21,1:ki_ssd)      ] -  [w_sulcal*Mw.C                 Mw.E.X(:,1:ki_ssd)             ]).^2,2));

            case 'partial'
                ssd_fwd(iter_match) = sum(sum(([w_sulcal*M0.C(idx1)           M0.E.X(idx1,1:ki_ssd)          ] -  [w_sulcal*Mw.C(idx2(c.corr12)) Mw.E.X(idx2(c.corr12),1:ki_ssd)]).^2,2));
                ssd_bwd(iter_match) = sum(sum(([w_sulcal*M0.C(idx1(c.corr21)) M0.E.X(idx1(c.corr21),1:ki_ssd)] -  [w_sulcal*Mw.C(idx2)           Mw.E.X(idx2,1:ki_ssd)          ]).^2,2));
        end
        ssd(iter_match) = ssd_fwd(iter_match) + ssd_bwd(iter_match);

        % Retain best transformation
        cur_ssd = ssd_fwd(iter_match) + ssd_bwd(iter_match);
        if cur_ssd < min_ssd
            best_R12 = R12;
            best_R21 = R21;
            min_ssd = cur_ssd;
        end

        % Check stopping criterium
        if iter_match>5  &&  ssd(iter_match) > ssd(iter_match-2)
            disp(['Stopping iterations']);
            disp(['Before : SSD ' num2str(ssd_fwd(1))          ' ' num2str(ssd_bwd(1))          ' ******']);
            disp(['After  : SSD ' num2str(ssd_fwd(iter_match)) ' ' num2str(ssd_bwd(iter_match)) ' ******']);
            break;  % stop if energy is increasing
        end

    end % for iter_match

    % Retain best transformation
    R12 = best_R12;
    R21 = best_R21;
    X2  = U2 * R12 * diag([1; M.E.V(1:krot)].^-.5);
    Mw.E.X(:,1:krot) = X2(:,2:end);
    if do_save;  disp('Saving embedding and transformation data'); end;
    if do_save;  save([save_path_trans '/transformation-' id '-to-' id0 '.mat'],'R12','R21','-v7.3');  end;
    if do_save;  save([save_path_embed id '.mat'], 'Mw', '-v7.3');  end;
    
    
    % Allocate the data for saving
    X = Mw.E.X(:,1:3); % aligned spectral coordinates
    U = Mw.E.U(:,1:3); % un-aligned spectral coordinates
    A = Mw.E.W;        % weighted adjacency
    EUC = Mw.X;        % Eucledian coordinates
    C = Mw.C;          % sulcul depth
    T = Mw.T;          % cortical thickness
    F = Mw.F;          % mesh face
    [Y, GT] = processparcels(Mw, path_to_manlabels);
    [SX, AG] = processphenotype(data_name, path_to_csvfile); 
       
    % Saving
    if do_save;  disp('Saving for py'); end;
    if do_save;  save([save_path_forpy id '.mat'], 'X', 'U', 'A', 'EUC', 'C', 'T', 'F','Y','GT','SX','AG'); end;
    
                
end  % iter_subject ids








    

