% [E] = embedding(M,ne,opts)
%    Compute the graph embedding of L,
%      which are eigenmodes weighted by \lambda^-1
%    M       - Mesh (with M.X, M.F)
%              If M.G exist, this is used as node weighting in the decomposition
%    k       - number of eigenvectors to compute, including null vector (default: 6)
%    mode    - type of graph laplacian, either:
%                'randomwalk': D^-1 (D-W)
%                'symmetric' : D^-.5 (D-W) D^-.5
%                'bare'      : (D-W)
%    opts    - options to pass to eigs solver, it is also passed to weightedadjacency
%              (default: 'disp=0','isreal=1','issym=0','maxit=50')
%    E       - Embedding with E.X,E.F, and sorted eigenvectors in columns (n x k)
%              and sorted eigenvalues (column vector)
%
function E = embedding(M,k,mode,opts)

    if nargin<2;                 k = 6;  end;
    if nargin<3;                 mode = 'randomwalk';  end;
    if nargin<4;                 opts = struct(); end;
    
    if ~isfield(opts,'disp');    opts.disp   = 0;     end;
    if ~isfield(opts,'isreal');  opts.isreal = 1;     end;
    if ~isfield(opts,'issym');   opts.issym  = 0;     end;
    if ~isfield(opts,'maxit');   opts.maxit  = 4000;  end;

    if ~isfield(M,'n');  M.n = size(M.X,1);  end;
    
    % Graph Laplacian
    switch mode
        case 'randomwalk'
            % Laplacian
            E.W = weightedadjacency(M.X,M.F,opts);
            [E.D,E.Dinv] = degree(E.W);
            E.L = E.Dinv * (E.D - E.W);
            
            % Node-weighted Laplacian
            if isfield(M,'G')
                E.Ginv = spdiags(M.G.^-1,      0,speye(M.n));
                E.L = E.Ginv * E.L;
            end
            
            % Spectrum
            [E.U,E.V] = spectrum(E.L,k,[],opts);
            
            % Embedding
            E.X = E.U;
            E.X(:,2:end) = E.U(:,2:end) * diag(abs(E.V(2:end)).^-.5);
            
        case 'symmetric'
            % Laplacian
            E.W = weightedadjacency(M.X,M.F,opts);
            [E.D,E.Dinv,E.Dsqrt,E.Dinvsqrt] = degree(E.W);
            E.L = E.Dinvsqrt * (E.D - E.W) * E.Dinvsqrt;

            % Node-weighted Laplacian
            if isfield(M,'G')
                E.Ginvsqrt = spdiags(M.G.^-.5,     0,speye(M.n));
                E.Gsqrt    = spdiags(M.G.^+.5,     0,speye(M.n));
                E.L = E.Ginvsqrt * E.L * E.Ginvsqrt;
            end
            
            % Spectrum
            opts.issym = 1;  % Faster
            [E.U,E.V] = spectrum(E.L,k,[],opts);
            
            % Embedding
            E.X = E.Dinvsqrt * E.U;
            if isfield(M,'G');  E.X = E.Ginvsqrt * E.X;  end;  % If Node-weighted
            mag = sqrt(sum(E.X.^2,1));
            E.X = E.X * diag(mag.^-1);
            E.X(:,2:end) = E.X(:,2:end) * diag(abs(E.V(2:end)).^-.5);

        case 'bare'
            % Laplacian
            E.W = weightedadjacency(M.X,M.F,opts);
            [E.D,E.Dinv] = degree(E.W);
            E.L = (E.D - E.W);

            % Node-weighted Laplacian
            if isfield(M,'G')
                E.Ginv = spdiags(M.G.^-1,      0,speye(M.n));
                E.L = E.Ginv * E.L;
            end
            
            % Spectrum
            [E.U,E.V] = spectrum(E.L,k,[],opts);
            
            % Embedding
            E.X = E.U;
            E.X(:,2:end) = E.U(:,2:end) * diag(abs(E.V(2:end)).^-.5);
    end
    
    
end
