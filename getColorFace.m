function cmat = getColorFace(cindx, info)
%
%  function cmat = getColors(cindx, inf)
%
%  Gets colors from table. Periodic over 20 colours. The purpose is to
%  generate text tags corresponding to plot colours that can be used in
%  titles, text,... .
%
%  This has a different color order so as to allow independent control over
%  the marker face color seperate from the marker/line color.
%
% INPUTS:
%       cindx:  index numbers from table (e.g. 2:6)
%        info:  indicates what info to return'num' returns matric whereas
%               'name' returns text names (for plotting)
% OUTPUTS:
%       cmat:   Either matrix of numbers (if 'info' = 'num') or cell array
%               of names (if 'info'='name')
%
% 
% examples:
%       getColor([1,4,6666],'num')
%
%           ans =
%               0         0    1.0000
%               0    0.7500    0.7500
%               0.7500    0.7500         0
%
%       getColor([1,4,6666],'name')
%
%           ans = 
%               'blue'    'cayan'    'gold'
%
%   Aaron Gruber   2011_9_19

%% define colors and their names
colororder = [
	0.00  0.00  0.00    
	1.00  1.00  1.00      
	0.00  0.50  0.00
	0.00  0.75  0.75
	0.75  0.00  0.75
	0.75  0.75  0.00
	0.25  0.25  0.25
	0.75  0.25  0.25
	0.95  0.95  0.00
	0.25  0.25  0.75
	0.75  0.75  0.75
	0.00  1.00  0.00
	0.76  0.57  0.17
	0.54  0.63  0.22
	0.34  0.57  0.92
	1.00  0.10  0.60
	0.88  0.75  0.73
	0.10  0.49  0.47
	0.66  0.34  0.65
	0.99  0.41  0.23
];

cnames = {'black', 'white', 'red', 'cayan',  'purple', 'gold', 'black', ...
         'maroon','yellow','dkblue','gray','ltgreen','tan','pea','ltblue'...
         'pink','peach','dkgreen','dkmaroon','ltred'};
     
%%
[numColors, ~] = size(colororder);
switch info
    case 'num'
        for k=1:length(cindx)
            ind = mod(cindx(k)-1,numColors)+1; % periodic over 1 to num colors
            cmat(k,:) = colororder(ind,:);
        end   
        
    case 'name'
        for k=1:length(cindx)
            ind = mod(cindx(k)-1,numColors)+1; % periodic over 1 to num colors
            cmat(k) = cnames(ind);
        end
        
    otherwise
        error('second argument must be either ''num'' or ''name'' ')
end

