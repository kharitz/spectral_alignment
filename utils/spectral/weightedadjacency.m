% W = weightedadjacency(X,F)
%   Compute weighted adjacency matrix
%   X - vertex position (n x ndim)
%   F - triangulation (n x 3, vertex indices for each triangle)
%   opts - structures with set different types of weighting
%          'kernel' - 'l2', affinity is inverse distance between nodes (default)
%                   - 'l2inv', affinity is distance between nodes
%                   - 'gaussian', affinity is distance in Gaussian kernel
%                        with options:
%                          'sigma', kernel width (default: 'max')
%                             scalar - scalar for every edge weight
%                             'max'  - max (weights) / 3
%                             'mean' - mean(weights) / 3
%          'eps'    - add a small epsilon to all weights (default: 0)
%                     penalizes small distances (avoid division by zero)
%
%   W - weighted adjacency matrix (n x n)
%       if nodes i and j are connected, W(i,j) gives their affinity
%       affinity is the inverse distance W(i,j) = dist(X(i,:) - X(j,:)))^-1
%
function W = weightedadjacency(X,F,opts)

    if nargin<3;  opts = struct();  end;
    if ~isfield(opts,'kernel');  opts.kernel = 'l2';  end;
    if ~isfield(opts,'eps');     opts.eps    = 0;     end;

    % compute weights for all links (euclidean distance)
    weights = [sum((X(F(:,1),:)-X(F(:,2),:)).^2,2); ...
               sum((X(F(:,1),:)-X(F(:,3),:)).^2,2); ...
               sum((X(F(:,2),:)-X(F(:,1),:)).^2,2); ...
               sum((X(F(:,2),:)-X(F(:,3),:)).^2,2); ...
               sum((X(F(:,3),:)-X(F(:,1),:)).^2,2); ...
               sum((X(F(:,3),:)-X(F(:,2),:)).^2,2)].^0.5;
           
    switch(opts.kernel)
        case 'l2'
            %weights = weights.^-1;  % (BUG which was enabled in earlier time, unknown when this was introduced)
            weights = (weights + opts.eps).^-1;  % penalize small distances (avoid division by zero)
            
        case 'l2inv'
            weights = (weights + opts.eps).^+1;  % penalize small distances (avoid division by zero)

        case 'gaussian'
            if ~isfield(opts,'sigma');  opts.sigma = 'max';  end;
            if isfield(opts,'sigma') && isequal(opts.sigma,'max');   opts.sigma = max (weights) / 3;  end;
            if isfield(opts,'sigma') && isequal(opts.sigma,'mean');  opts.sigma = mean(weights) / 3;  end;
            weights = exp(-(weights + opts.eps).^2 / (2*opts.sigma^2));
            
        otherwise;  error(['Unknown kernel type: ' opts.kernel]);  
    end
    
    % remove duplicated edges
    rows = [F(:,1); F(:,1); F(:,2); F(:,2); F(:,3); F(:,3)];
    cols = [F(:,2); F(:,3); F(:,1); F(:,3); F(:,1); F(:,2)];
    [rc,idx] = unique([rows,cols], 'rows','first');
    weights = weights(idx);

    % convert to double if necessary
    if ~isa(rc,     'double');  rc      = double(rc);       end;
    if ~isa(weights,'double');  weights = double(weights);  end;

    % fill adjacency matrix
    W = sparse(rc(:,1), rc(:,2), weights);
    
end

