% [U,Lambda] = spectrum(L,ne,opts)
%    Compute the graph spectrum of L
%    L       - general graph Laplacian
%    k       - number of eigenvectors to compute, including null vector (default: 6)
%    v0      - initial eigenvector for eigs (default: [])
%    opts    - options to pass to eigs solver
%              (default: 'disp=0','isreal=1','issym=0','maxit=50')
%    E       - sorted eigenvectors in columns (n x k)
%    lambda  - sorted eigenvalues (column vector)
%
function [U,Lambda] = spectrum(L,k,v0,opts)

    if nargin<2;                 k = 6;  end;
    if nargin<4;                 opts = struct('disp',0,'isreal',1,'issym',0,'maxit',4000); end;
    if nargin>2 && ~isempty(v0); v0 = sum(v0,2); opts.v0 = v0 / norm(v0);  end;

    disp(['*** Spectral decomposition of a ' num2str(size(L,1)) '^2 matrix']);
    
    if  isfield(opts,'srnd');        srnd(opts.srnd);         end;
    if ~isfield(opts,'eigs_sigma');  opts.eigs_sigma = 'sm';  end;
    
    
    % Eigenvectors
    [U,Lambda,flag] = eigs(L,k,opts.eigs_sigma,opts);
    Lambda = diag(Lambda);
    
    % Failsafe (not interested in imaginary part)
    U = real(U);
    Lambda = real(Lambda);
    
    % Sort ascending order
    [~,idx] = sort(Lambda,1,'ascend');
    Lambda = Lambda(idx);
    U = U(:,idx);
    
    if flag;  warning('Spectral decomposition did not converge');  end;
    if sum(Lambda<1e-10)>1;  warning('Multiple null vectors');  end;
    
    % Arbitrary signing
    signf = 1 - 2*(U(1,:)<0);  % first point in X always has positive eigenmode values
    U = bsxfun(@times,U,signf);

    
end
