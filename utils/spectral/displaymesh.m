function h = displaymesh(X,F,C)

    if nargin<3;  C = 1:size(X,1);  end;
    
    h = trisurf(F, X(:,1), X(:,2), X(:,3), C);  
    shading interp;  daspect([1 1 1]);  axis on;  material dull;  
    xlabel('x');  ylabel('y');  zlabel('z');
    view(0,0);  
    
    if numel(findobj(gca,'type','light')) == 0
        camlight;
        link_light;
    end
    
    set(gca,'cameraviewanglemode','manual');
    
    if nargout < 1;  clear h;  end;
end
