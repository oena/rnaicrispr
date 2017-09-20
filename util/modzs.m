function [czs, norm_wt, cc] = modzs(zs, ridx, varargin)
% MODZS Compute a moderated Z-score.
% CZS = MODZS(ZS, RIDX)

pnames = {'clip_low_wt', 'clip_low_cc',...
          'low_thresh_wt', 'low_thresh_cc', 'metric'};
dflts = {false, false,...
         0.01, 0, 'wt_avg'};
args = parse_args(pnames, dflts, varargin{:});

[~, nc] = size(zs);
if nc > 1    
    cc = fastcorr(zs(ridx, :), 'type', 'spearman');    
%     if args.weighted_cc
%         % apply a sigmoidal weighting function
%         cc = sigmoid(cc, -1, 0.4, 0.1);                
%     end
    cc = cc - diag(diag(cc));    
    if args.clip_low_cc
        % clip low cc's
        cc = clip(cc, args.low_thresh_cc, inf);
    end
    
    wt = 0.5*sum(cc, 2);
    if args.clip_low_wt
        % clip low weights
        wt = clip(wt, args.low_thresh_wt, inf);
    end
    
    switch(lower(args.metric))
        case 'wt_avg'
            sum_wt = sum(abs(wt));
        case 'wt_stouffer'
            sum_wt = sqrt(sum(wt.^2));
    end
    norm_wt = wt / sum_wt;
    czs = zs * norm_wt;
else
    czs = zs;
    norm_wt = 1;
    cc = 1;
end
end