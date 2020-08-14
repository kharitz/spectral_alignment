% D = degree(W)
%   Compute the degree matrix from an adjacency
%   W    - a matrix
%   D    - diagonal degree matrix (each elt is the column sum of W)
%   Dinv - inverse diagonal degree matrix
%   Dsqrt    - squared diagonal degree matrix
%   Dinvsqrt - inverse squared diagonal degree matrix
%
function [D,Dinv,Dsqrt,Dinvsqrt] = degree(W)

    n    = size(W,1);
    D    = spdiags(sum(W)',    0,speye(n)); %D = diag(sum(W));
    
    if nargout>1
        Dinv = spdiags(sum(W)'.^-1,0,speye(n)); %D = diag(sum(W));
    end
    
    if nargout>2
        Dsqrt = spdiags(sum(W)'.^+.5,0,speye(n)); %D = diag(sum(W).^.5);
    end
    
    if nargout>3
        Dinvsqrt = spdiags(sum(W)'.^-.5,0,speye(n)); %D = diag(sum(W).^-.5);
    end
    
end
