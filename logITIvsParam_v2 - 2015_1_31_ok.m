function d = logITIvsParam_quant(param, TR_c, bins_RealT, n_thresh, regress_thresh_time, varargin)
%%
% this version divides data into quantiles - so each point has same # pts
%
%logITIvsParam(param, {TR}, binsRealT, n_thresh, regress_thresh_time)
%
% computes a semilog regression of inter-trial-interval against variable
% passed in as 'param', in each vector of trial objects {TR}. Returns
% fitting coefficients, stats, and data points.
%
% INPUTS:
%       param = variable name for property to regress against ITI
%       {TR}  = cell array of vector of trial objects; 
%       bins_RealT = time bins for computation in ms
%       n_thresh = threshold number of data points to use for computing averages
%       regress_thresh_time = upper limit in time for regression
%
%       options:
%
%           -QuantBins      <n> divides distribution for each vector of
%                           objects in TR into n equal percentiles. Ignores
%                           the vecotr 'binsRealT'.
%
%           -GroupingProp   <none (default) | Rat_id | Group_id | Drug | or other property> 
%                           color code each TR_c group by field. Has to be consistent
%                           in every trial in each TR_c passed in. 
%
%           -CumSumPltVec   Vector of bins for cumultive sum plot (e.g. [100:100:50000])
%
%           -PolyOrder      order of polynomial for fitting; if not set, then linear regression   
%
%           -PlotEachGrp    plot distributions and fits for each group
%                           (e.g. animal) in the input cell
%
%
%
% OUTPUTS:
%       d: structure with one instance for each group of TR in input cell
%       array
%
%
%  2015_1_20;    Aaron Gruber

%% chagelog
%
%  2015_1_23 added options for producing group-specific colors
%  2015_1_23 added option for using equal time bins (default; has unequal samples/bin) 
%            or variable bins with euqual numbers of samples (with -QantBins)

%% parse options
grp_prop = 'none';  % default property for (optional) grouping by color in plots 
nBinsQuant = 0;     % parameter used for -QuantBins that determines the number of bins with equal samples.
CumSumPltVec = [100:100:50000]; % vector of bins for cumultive sum distribtion plot
PolyOrder = 0;  % specify the polynomial model for fitting data; if = 0, then do linear regression
Plt_each_flg = 0;

if nargin > 5; 
    varargin_txt = varargin;        
    varargin_txt(~cellfun(@ischar, varargin)) = {'placeholder'};
   optIndx = find(cellfun(@isempty, strfind(varargin_txt,'-'))==0);
   opts_cell = varargin_txt(optIndx);
   for opt=opts_cell
       switch opt{:}
           case '-GroupingProp'
               indx = find(strcmp('-GroupingProp',varargin_txt)==1);
               grp_prop = varargin{indx+1};
           case '-CumSumPltVec'
               indx = find(strcmp('-CumSumPltVec',varargin_txt)==1);
               CumSumPltVec = varargin{indx+1};
           case '-PolyOrder'
               indx = find(strcmp('-PolyOrder',varargin_txt)==1);
               PolyOrder = varargin{indx+1};
           case '-PlotEachGrp'
               Plt_each_flg = 1;    
           case '-QuantBins'
               indx = find(strcmp('-QuantBins',varargin_txt)==1);
               if(indx<numel(varargin_txt));
                    nBinsQuant = varargin{indx+1};
               else
                    nBinsQuant = 0; % if none passed in, set to 0, will go to devault in next line
               end
               if(nBinsQuant < 1 || ~isnumeric(nBinsQuant)) % this takes care of invalid values (e.g. <1)
                   warning(['No # of bins specified in -QuantBins option, using 5 as default'])
                   nBinsQuant = 5;
               end
            otherwise
               if(strcmp(opt{:}(1),'-'))
                   error([opt{:}, ' is not a valid option']);
               end
       end 
   end
end


if (strcmpi(grp_prop, 'none'))
    grp_indx = ones(1,numel(TR_c));
    grp = {};
else
    for i=1:numel(TR_c)
        grp{i} = char(TR_c{i}(1).(grp_prop)); % get the value of the group prop for the first object of each TR_c group of objects
        if(numel(unique({char(TR_c{i}(1).(grp_prop))})) > 1) % throw an error if there is more than one level of the property
            error(['GroupingProp ',grp_prop, ' is invalid for inpu TR num: ', num2str(i), 'because there are more than one level present in the object array'])
        end
    end
    [~,~,grp_indx] = unique(grp); % get the index of the unique values
end

%% do analysis
figure
ax(1) = subplot(2,2,1); ylabel('count'); xlabel('ITI')
ax(2) = subplot(2,2,2); ylabel('cum PDF'); xlabel('ITI');
ax(3) = subplot(2,2,3); ylabel(param); xlabel('ITI'); %
ax(4) = subplot(2,2,4); ylabel('# trials'); xlabel('ITI')
ax(5) = axes('position', [0.8 0.25 0.15 0.15]); ylabel('slope'); xlabel('ITI');
set(ax(5),'NextPlot','add')
set(ax([1:4]),'XScale','log','XLim',[1,100],'NextPlot', 'add');

for i=1:numel(TR_c) % loop over all arrays of objects passed in
    TR = TR_c{i};
    
    % create counts of ITI accoridng to either bins passed in or quantize them
    
    if(nBinsQuant>0)
        bins_RealT = quantile([TR.InterTrialInterval], nBinsQuant);
        bins_RealT = [bins_RealT, bins_RealT(end)+bins_RealT(1)];
    end
    
    [nITI, indx_bin] = histc([TR.InterTrialInterval], bins_RealT);

    stairs(ax(1), bins_RealT(1:end-1)./1000, nITI(1:end-1),'Color',getColor(grp_indx(i),'num')); hold on
    
    tmp = num2cell(indx_bin);
    [TR.Grouped_trial_indx] = tmp{:}; % set the 'Grouped_trial' property to the index of the hist bin; this controls which elements are averaged (i.e. those in each bar of histogram)

    groups = {'Grouped_trial_indx'};

    % make a dataset to produce the averaging over the grouped trial (hist bins)
    ds = makeDS(param, TR, groups); 
    ds_ok = ds(ds.GroupCount > n_thresh , :); % get the bins with samples > threshold
    [numOK,~] = size(ds_ok);
    if(numOK <4 )
        warning('not enough samples to compute regression')
        [B,BINT,R,RINT,STATS] = deal([nan;nan],nan,nan,nan,nan); % have to put nan so there is no indexing problem later that comes with concatonation with an empty cell
    else
        %xvals = bins_RealT(double(ds_ok.Grouped_trial_indx) + 1); % since we had indicies, find corresponding time bins
        xvals = bins_RealT(str2num(char(ds_ok.Grouped_trial_indx)) + 1); % due to some bug, this works; str2num(char(ds_ok.Grouped_trial_indx(:)))+1

        plot(ax(3), xvals./1000, ds_ok.(param), 'o','Color',getColor(grp_indx(i),'num')); hold on
        % regression

        % transform to log10(seconds) basis
        x_o = log10(xvals./1000)'; % in log base
        y_o = ds_ok.(param); % the y vals
        
        % restrict time for regression
        indx = find(x_o > 0 & x_o<log10(regress_thresh_time/1000));
        x = x_o(indx);
        y = y_o(indx); 

        if(PolyOrder>0) % do polynomial fit if given order
            % first an unweighted polynomial regression
            p = polyfit(x,y,PolyOrder);
            x1 = linspace(x(1),x(end));
            f1 = polyval(p,x1);
            
           % plot(ax(3),10.^x1,f1,'y:')
            
            
            % weighted nonlinear regression
            modelFun = @(p,x) polyval(p,x);
            nlm = fitnlm(x,y,modelFun,p); % this gives the same as polyfit
            w = ds_ok.GroupCount(indx)./10.^x; % weighting function - number of samples per point
           % w = (ds_ok.GroupCount(indx)).^1.5./x.^2; 
            %w = ds_ok.GroupCount(indx)./10; 
            nlm_w = fitnlm(x,y,modelFun,p,'Weight',w); 
            p_w = nlm_w.Coefficients.Estimate;
                %f2 = polyval(p_w,x1);
                %plot(ax(3),10.^x1,f2,'Color',[0 0 0])
            [ypred,ypredci] = predict(nlm_w,x1','Simultaneous',true);
            plot(ax(3),10.^x1',ypred,'Color',getColor(grp_indx(i),'num')); 
            plot(ax(3),10.^x1',ypredci,'Color',getColor(grp_indx(i),'num'),'LineStyle',':'); 
            
            d(i).B = p_w';
            
            if(Plt_each_flg)  % if input flag thrown, plot each group independenty in a new fig
                figure; axx(1) = subplot(2,1,1); axx(2) = subplot(2,1,2); set(axx,'XScale','log','XLim',[1,100],'NextPlot', 'add')
                stairs(axx(1), bins_RealT(1:end-1)./1000, nITI(1:end-1)); % plot histogram
                xlabel(axx(1),'ITI'); ylabel('count')
                plot(axx(2), xvals./1000, ds_ok.(param), 'o'); hold on    % plot mean for each bin  
                plot(axx(2), 10.^x1',ypred,'b-', 10.^x1',ypredci,'b:');   % plot the fit
                xlabel(axx(2),'ITI'); ylabel(axx(2),param)
                if(strcmpi(grp_prop,'none'))
                    title(['data and fit for Rat =', char(TR_c{i}(1).Rat_id), '; coefs = ',num2str(p_w')])
                else
                    title(['data and fit for ',grp_prop,' = ', char(TR_c{i}(1).(grp_prop)), '; coefs=',num2str(p_w')])
                end
            end
            % stats
            %[p_w, F_w] = coefTest(nlm_w);
            
            
            if(i<2)
                title(ax(3),{['Polynomial fit (order=',num2str(PolyOrder),'); df=',num2str(nlm_w.DFE),'; R^2 adjusted=',num2str(nlm_w.Rsquared.Adjusted)];...
                      ['for F and p see the screen output']}); % ['F=',num2str(F_w),'; p=',num2str(p_w)]})
                   
                nlm_w.disp % show that stats on screen
            end
            
            
            
            %plot(ax(3),10.^x1,f1,'Color',getColor(grp_indx(i),'num'),'LineStyle','-.')
            
        else            % otherwise do a linear regression
            % scatter plot & fit
            [B,BINT,R,RINT,STATS] = regress(y, [ones(numel(x),1), x]);
            xlims = [max(x), min(x)];
            d(i).B = B; 
            d(i).STATS = STATS;
            plot(ax(3),[10.^xlims], [B(1)+xlims(1)*B(2), B(1)+xlims(2)*B(2)],'Color',getColor(grp_indx(i),'num'),'LineStyle',':'); hold on% put the x vals back to regular space so in log axis they match data points
        end
        
        if(~isempty(CumSumPltVec)) % make cumulitive sum plot
                [nITI_cum, indx_bin_cum] = histc([TR.InterTrialInterval], CumSumPltVec); % use the the other vector for cumultive sum
                nITI_CDF = cumsum(nITI_cum)./sum(nITI_cum); % make CDF
                stairs(ax(2), CumSumPltVec(1:end-1)./1000,nITI_CDF(1:end-1),'Color',getColor(grp_indx(i),'num')); hold on
        end
    end
     
    d(i).data_var = param;
    d(i).Rat_id = TR(1).Rat_id;
    d(i).Group_id = TR(1).Group_id;
    d(i).Drug = TR(1).Drug;
    d(i).Grp_prop = grp_prop; 
    d(i).ds = ds; d(i).xvals = xvals;
    d(i).nITI_CDF =nITI_CDF ;
    d(i).paramData = ds_ok.(param);
    d(i).LogBinedTime = x_o;
    d(i).grp_indx = grp_indx(i);
    d(i).color = getColor(grp_indx(i),'num');
    if(~strcmpi(grp_prop,'none'))
        d(i).Grp_prop_val = d(i).(grp_prop);
    end
end

% do KS test on CDFs and print results
ik = 1; kstxt{ik} = ['2 sample KStest'];
if(numel(d)>1 & numel(d)<8)
    for i=1:numel(d)-1
        for k=i+1:numel(d)
            [H,CDF_p(ik),KSSTAT] = kstest2([d(i).nITI_CDF], [d(k).nITI_CDF]);
            kstxt{ik} = [char(d(i).Grp_prop_val), '-', char(d(k).Grp_prop_val),'; P(ks) = ',num2str(CDF_p(ik))];            
            ik = ik+1;

        end
    end
end


axes(ax(2))
kstxt_h = text(11, 0.6, kstxt,'FontSize',7);


% put a legend if more than one group
if(numel(grp)>0)
    l_h = legend(ax(1),grp); set(l_h,'FontSize',7);
end

if(~isnan([d.B]))
    slope = [d.B];
    [H,P,CI,stats] = ttest(slope(2,:)); 
end

if(nBinsQuant > 0)
    title(ax(1),['quant method for bins', '   ; nSamples Threshold: ',num2str(n_thresh),  ])
elseif(PolyOrder==0)
    title(ax(1),['delta t: ', num2str(bins_RealT(2)-bins_RealT(1)), '   ; nSamples Threshold: ',num2str(n_thresh), '; tstat (all pts): ',num2str(stats.tstat), '; p: ',num2str(P), '; df: ', num2str(stats.df)  ])
else
    title(ax(1),['delta t: ', num2str(bins_RealT(2)-bins_RealT(1)), '   ; nSamples Threshold: ',num2str(n_thresh), '; Polynomial Order: ',num2str(PolyOrder)])
end

if(PolyOrder > 0)
    BforHterm = slope(1,:); % beta for the highest order term (e.g. x^2 for parabolic - dign determines if it is pointing up or down.
    bins = linspace(min(BforHterm), max(BforHterm), 20);
else
    bins = [-1:0.05:1];                   % bins for plotting dist of slopes
end

if( numel(grp)<2) % if 0 or 1 group, plot all of the coefficients
    [nSlope, bins]=hist(BforHterm,bins);
    stairs(ax(5),bins,nSlope); hold on; set(ax(5),'XLim',[min(BforHterm), max(BforHterm)])
    title(ax(5),['tstat: ',num2str(stats.tstat), '; p: ',num2str(P), '; df: ', num2str(stats.df) ])
else              % plot hist for each group in a different colour
    for i=1:numel(unique(grp))            % loop over the number of groups
        indxGpIndx = find(grp_indx == i); % find indicies where group is the same
        [nSlope, bins]=hist(slope(2,indxGpIndx),bins);
        stairs(ax(5),bins,nSlope,'Color',getColor(grp_indx(indxGpIndx(1)),'num')); hold on
        title(ax(5),['tstat: ',num2str(stats.tstat), '; p: ',num2str(P), '; df: ', num2str(stats.df) ])
    end
end


%set(ax(5),'XTick',[-1:0.5:1])
xlabel(ax(5),'regression coefficient')


for i=1:numel(TR_c)
    edges = quantile([TR_c{i}.InterTrialInterval]./1000, [0.20 0.8]);
    median_vec(i) = median([TR_c{i}.InterTrialInterval])./1000;
    numTR_vec(i) = numel(TR_c{i});
    plot(ax(4), edges, [numel(TR_c{i}), numel(TR_c{i})],'Color',getColor(grp_indx(i),'num')); hold on
    plot(ax(4), median([TR_c{i}.InterTrialInterval])./1000, numel(TR_c{i}), 'o','Color',getColor(grp_indx(i),'num'))
end

if(numel(median_vec)>1) % do a regression on the # trials vs iti if more than one point
    y = log10(median_vec)'; 
    x = numTR_vec';

    [B,BINT,R,RINT,STATS] = regress(y, [ones(numel(x),1), x]);
    xlims = [max(x), min(x)];
    plot(ax(4), 10.^[B(1)+xlims(1)*B(2), B(1)+xlims(2)*B(2)], [xlims],'m:'); hold on
    title(ax(4),['R=',num2str(STATS(1)),'; F=',num2str(STATS(2)),'; P=',num2str(STATS(3))])
end


xlabel('ITI')

%'hi'





