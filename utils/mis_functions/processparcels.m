function [y_one_hot,man_one_hot] = processparcels(M, path)
%Process the FreeSurfer or manual labels
%  Removes certain parcels and standardizes
if exist(path,'file')
    man_lab= vtk_polydata_read(path);
    mangt = man_lab.point_data.data;
    mangt = mangt-1000;
    mangt(mangt==-1001)=0;
    mangt(mangt==33) = 0;
    mangt(mangt==32) = 0;
    mangt(mangt==1) = 0;
    mangt(mangt==4) = 0;
    if length(unique(mangt))==31
        mangt(Mw.P==0)=0;
    end
    [~, loc] = ismember(mangt, unique(mangt));
    man_one_hot = full(ind2vec(loc')');
    Y = M.P;
    Y = Y-1; Y(Y==-1)=0;
    Y(Y==1)=0;Y(Y==32)=0;Y(Y==33)=0;
    [~, loc] = ismember(Y, unique(Y));
    y_one_hot = full(ind2vec(loc')');
    
else
    man_one_hot = [];
    Y = M.P;
    Y = Y-1; Y(Y==-1)=0;
    Y(Y==1)=0;Y(Y==32)=0;Y(Y==33)=0;
    [~, loc] = ismember(Y, unique(Y));
    y_one_hot = full(ind2vec(loc')');
    
end


end

