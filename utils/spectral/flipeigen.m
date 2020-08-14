% M2 = flipeigen(M1,M2)
%   Flip eigenvector signs of E2 such that they match those of E1
%   Signing over meshes should be coherent spatially (matching with barycenter
%   of eigenvector poles)
%
%     M1 - Reference embedding
%     M2 - Embedding to flip sign
%     ne - Number of eigenvectors to check sign
%
% herve.lombaert@gmail.com - Jan. 2015
%
function M2 = flipeigen(M1,M2,ne)

    if nargin < 3;    ne = min(M1.E.ne, M2.E.ne);  end;
    
    % Find centers of posivite/negative modes - check spatial consistency
    for iter_eigen = 1:ne
        %
        offset1 = mean(M1.X(:,1:3));
        offset2 = mean(M2.X(:,1:3));
        X1 = bsxfun(@minus, M1.X(:,1:3), offset1);  % centered mesh
        X2 = bsxfun(@minus, M2.X(:,1:3), offset2);  % centered mesh

        % Weighted barycenters of poles
        w1 = sign(M1.E.X(:,iter_eigen)).*abs(M1.E.X(:,iter_eigen).^3);
        w = w1;  w(w<0) = 0;  w = w/sum(w);  
        avgX1p = sum(bsxfun(@times,X1,w));  % positive pole barycenter
        w = w1;  w(w>0) = 0;  w = w/sum(w);  
        avgX1m = sum(bsxfun(@times,X1,w));  % negative pole barycenter

        w2 = sign(M2.E.X(:,iter_eigen)).*abs(M2.E.X(:,iter_eigen).^3);
        w = w2;  w(w<0) = 0;  w = w/sum(w);  
        avgX2p = sum(bsxfun(@times,X2,w));  % positive pole barycenter
        w = w2;  w(w>0) = 0;  w = w/sum(w);  
        avgX2m = sum(bsxfun(@times,X2,w));  % negative pole barycenter

        % distances between matched poles
        distp = sum((avgX1p - avgX2p).^2) + sum((avgX1m - avgX2m).^2);  % check where centers are spatially
        distm = sum((avgX1p - avgX2m).^2) + sum((avgX1p - avgX2p).^2);  % choose the one combination that makes more sense (red spot should be closer to each others

        if distm < distp
            if do_verbose;  disp(['Flip ' num2str(iter_eigen)]);  end;
            M2.E.U(:,iter_eigen) = -M2.E.U(:,iter_eigen);
            M2.E.X(:,iter_eigen) = -M2.E.X(:,iter_eigen);
        end
    end

end