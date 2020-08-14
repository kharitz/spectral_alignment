% displayeigenmaps(X, F, E, pos, total_rows, nbeigenvectors, ishorizontal)
%
% display eigenmaps on a row
%
%  X - vertices
%  F - faces
%  E - eigenmaps
%  pos - on which row (or colum) to display the eigenmaps
%  total - total number of columns (if we will be showing other eigenmaps)
%  nbeigenvectors - number of eigenvectors to show (default, show all vectors)
%  isvertical - show horizontally the eigenmodes (default, 1)
%
% Herve Lombaert, April 28th, 2010

function displayeigenmaps(X, F, E, pos, total_rows, nbeigenvectors, ishorizontal)

if nargin<4; pos            = 1;  end;
if nargin<5; total_rows     = 1;  end;
if nargin<6; nbeigenvectors = size(E,2); end;
if nargin<7; ishorizontal   = 1;  end;

if ~ishorizontal

    % Vertically
    for i = 1:nbeigenvectors
        subplot(nbeigenvectors,total_rows,pos + (i-1)*total_rows);
        if isempty(F)
            plot4(X(:,1), X(:,2), X(:,3), E(:,i), '.');
        else
            trisurf(F, X(:,1), X(:,2), X(:,3), E(:,i));
            shading interp;
        end
        axis image; axis off;
        title(['Eigen ' num2str(i)]);
    end
    
else
    
    % Horizontally
    for i = 1:nbeigenvectors
        subplot(total_rows,nbeigenvectors,(pos-1)*nbeigenvectors + i);
        if isempty(F)
            plot4(X(:,1), X(:,2), X(:,3), E(:,i), '.');
        else
            trisurf(F, X(:,1), X(:,2), X(:,3), E(:,i));
            shading interp;
        end
        axis image; axis off;
        title(['Eigen ' num2str(i)]);
    end
    
end

end