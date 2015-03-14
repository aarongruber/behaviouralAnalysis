function E=CompEntropyOnBinSeq(x,len)
%
% E=CompEntropyOnBinSeq(x,len)
%
% Compute entropy (corrected for finite sample) from a binary sequense, parsing the sequense into chunks
% of lentgh 'len'. 
%
% INPUTS:
%       x:  sequense of binary digits or chars (should be column if eash indicator is more than one
%           character)
%     len:  number of consecutive elements of x to group together to form
%           lexicon (e.g. len=3 would reduce x into a new sequense of
%           lentgh numel(x)/len, whith each element coding for 3
%           consecutive elements in x)



x =reshape(x,numel(x),1);
x = grp2idx(x)-1;
if(numel(unique(x))>2)  % this check may not be needed - might work fine with more levels
    error('Input string has too many levels to compute entropy with this function -- make a different parser'), 
end

numTruncEl = length(x)-mod(length(x),len); % find correct length (divisible by 'len')
for k=1:len
    xm(k,:)=x(k:len:numTruncEl);           % deal into a m x len vec
end
%xm=xm';
xmc=num2cell(xm');               
xmcs = cellfun(@(x) num2str(x),xmc);        % convert to array of strings
x_dec =  bin2dec(xmcs);                     % convert to decimal

% %% this is for a newer entropy estimator - ESB: Nemenman, Bialex, deRuyter van Steveninck, 'Entropy and information in 
% %                                       neural spike trains: Progress on the sampling problem' PHYSICAL REVIEW E 69, 056111(2004)
%  
% [h_name,h_freq] = hist(nominal(x_dec)); % cast data as categorical, then histogram counts frequency.
% K = numel(x_dec);
% qfun=1;
% precision=.1;
% try
%     [S_nsb, dS_nsb, S_cl, dS_cl, xi_cl, S_ml,errcode]=find_nsb_entropy (double(h_freq), h_name, K, precision,qfun);
%     E = S_nsb;
% catch
%     E = NaN
%    disp('-----------------Cant find entropy----------------')
% end


%% this code is for the biased estimator with the classic Miller-Madow estimator: outdated
E = Entropy(x_dec);
k=2^len;
E = E + (k-1)/(1.3863*length(x));   % this is a correction due to finite samples; per Lee et al 2004, 