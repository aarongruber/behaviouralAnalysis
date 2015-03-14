function ca = dataset2cell_convertNominal(ds)
%
%  cellArray = dataset2cell_convertNominal(ds)
%   Convert a dataset into a cell array, but converting nominals to string
%

ca = dataset2cell(ds);
ca = cellfun(@convertNominal, ca, 'UniformOutput',false);


function y = convertNominal(x)
    if(isnumeric(x))
        y = x;
    else
        y = char(x);
end


