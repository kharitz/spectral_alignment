% do_multithread(0|1)
%    Global constant variable
%       Set with:    do_multithread(1);
%       Query with:  if do_multithread;  disp('multithreading');  end;
%
function o = do_multithread(i)

    if ~exist('matlabpool','file');
        if nargin>0;  warning('Matlabpool not present');  end;
        o = 0;  
        return;
    end

    persistent varval;
    %if isempty(varval);  varval = 1;  activate(varval);  end;
    if isempty(varval);  varval = 0;  end;
    if nargin>0;         varval = i;  activate(varval);  end;
    o = varval;
    
end

function activate(i)
    if i  % Activate
        if matlabpool('size')==0;  matlabpool;  end;
    else  % Deactivate
        if matlabpool('size')>0;   matlabpool close;  end;
    end
end
