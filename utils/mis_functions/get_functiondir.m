%  dir = get_functiondir(functioname)
%     Get directory of functioname.m
%
%  Herve Lombaert, Jan. 4th, 2013
%
function dir = get_functiondir(filename)
    if ~isempty(filename)
        [dir,name,ext] = fileparts(which(filename));
    else
        dir = pwd;
    end
end
