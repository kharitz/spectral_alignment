% Dinv = invdegree(W)
%   Compute the inverse degree matrix from an adjacency
%   W - a matrix
%   Dinv - diagonal inverse degree matrix (each elt is the column sum of W inversed)
%
function Dinv = invdegree(W)

    n = size(W,1);
    d = sum(W)';
    Dinv = spdiags(d.^-1,0,speye(n)); %D = diag(sum(W));
    
end
