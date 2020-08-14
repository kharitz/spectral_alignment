% L = generallaplacian(W,G)
%   Compute the general Laplacian
%   W - weighted adjacency matri
%   G - degree matrix (default: G=D, degree matrix)
%   L = G^-1 * (D-W);
function L = generallaplacian(W,G)

    D = degree(W);
    
    if nargin<2; G = D; end;
    
    L = G^-1 * (D-W);

end
