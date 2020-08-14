function [SX, AG] = processphenotype(data_name, path) 
% Extract age and sex of the subject
    if (strcmp(data_name, 'mindboggle'))
        mindsourceinfo = readtable(path);
        [idx,~] = find(strcmp(no_file(id).name(1:end-7),mindsourceinfo{:,1})==1);
        if (strcmp(mindsourceinfo{idx,2},'M'))
           SX = 0;
        else
            SX = 1;
        end
        AG = mindsourceinfo{idx,4};
        
    end
    if (strcmp(data_name, 'mindboggle'))
        % Load the phenotype data for classification
        label_table = readtable(path);
        [idx,val]= find(strcmp(files_pat(i).name,label_table{:,:})); 
        idx = idx(1);  
        %sex
        if(strcmp(label_table{idx,4}{:},'F'));SX=1; elseif(strcmp(label_table{idx,4}{:},'M'));SX=0; end;       
        %age      
        AG = str2num(label_table{idx,5}{:});        
    end
end

