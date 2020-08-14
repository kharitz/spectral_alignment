% [corr12,corr21] = match(E1,E2)
%   Find correspondences between embeddings E1 and E2
%     E1 - embedding 1
%     E2 - embedding 2
%     corr12     - correspondence from E1 onto E2 (n1 x 1 column vector)
%     corr21     - correspondence from E2 onto E1 (n1 x 1 column vector)
%
function [corr12,corr21] = match(E1,E2)

%    if ~do_verbose;  disp(['*** Matching ' num2str(size(E1,1)) ' x ' num2str(size(E2,1)) ' points in ' num2str(size(E1,2)) 'D']);  end;

%     mytic;
    
    % Find correspondences betwen 1 onto 2
    if nargout>0
        corr12 = zeros(size(E1,1),1);
        tree = kdtree(E2);
        if do_multithread
            n = matlabpool('size'); % split the transformation into many smaller computations
            parfor i=1:n
                M{i} = kdtree_closestpoint(tree, E1(i:n:end,:));
            end
            for i=1:n
                corr12(i:n:end,:) = M{i};
            end
            clear 'M';
        else
            corr12 = kdtree_closestpoint(tree, E1); % all corresponding points from embedding 1 onto embedding 2
        end
    end

    % Find correspondences betwen 2 onto 1
    if nargout>1
        corr21 = zeros(size(E2,1),1);
        tree = kdtree(E1);
        if do_multithread
            n = matlabpool('size'); % split the transformation into many smaller computations
            parfor i=1:n
                M{i} = kdtree_closestpoint(tree, E2(i:n:end,:));
            end
            for i=1:n
                corr21(i:n:end,:) = M{i};
            end
            clear 'M';
        else
            corr21 = kdtree_closestpoint(tree, E2); % all corresponding points from embedding 2 onto embedding 1
        end
    end

%     mytoc;
    
end


