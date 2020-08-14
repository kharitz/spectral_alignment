% A = adjacency(F)
%   Build adjacency matrix from triangle list
%      F - triangle list
%      A - adjacency matrix
%
function A = adjacency(F)

    n = max(F(:));
    
    % remove duplicated edges
    rows = [F(:,1); F(:,1); F(:,2); F(:,2); F(:,3); F(:,3)];
    cols = [F(:,2); F(:,3); F(:,1); F(:,3); F(:,1); F(:,2)];
    [rc,idx] = unique([rows,cols], 'rows','first');

    % fill adjacency matrix
    A = sparse(rc(:,1),rc(:,2),1,n,n);    

end
