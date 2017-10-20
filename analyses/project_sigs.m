function [out, brewresids] = project_sigs(prefs, ds, varargin)

% [out, brewresids] = project_sigs(ds, varargin)
%    'groupvars' - cell array containing elements of ds.chd with which to group signatures
%    'outdir' - string, output directory
%    'dataname' - identifier to name the output files
%    'mkfigs' - bool, whether to make figures
%    'n_random' - number of random permutations to run for elt_proj; set to at least 200
%    'rds' - roast dataset
%    'rcol' - identifier in ds.chd for corresponding roast signatures

params = {'groupvars', ...
    'outdir', ...
    'dataname', ...
    'mkfigs', ...
    'n_random', ...
    'rds', ...
    'rcol'};
dflts = {'', ...
    prefs.outdir, ...
    '', ...
    1, ...
    200, ...
    '', ...
    ''};
args = parse_args(params, dflts, varargin{:});

if isempty(args.dataname)
  error('dataname not supplied');
end

if ~exist(fullfile(args.outdir, args.dataname), 'dir')
  mkdir(fullfile(args.outdir, args.dataname));
end

[out, brewresids] = projection_iterator(ds, args);

if args.mkfigs
  mkfigs(out, args);
end

end


function mkfigs(out, args)
  f = out.elt_proj > -2;

  figure('Position', [10 10 1000 800]);
  hist(out.elt_rank(f), 100); grid on;
  xlabel('On-target CGS projection rank');
  ylabel('Count');
  title({sprintf('Distribution of on-target projection ranks relative to null, %s', args.dataname); ...
      sprintf('N sigs = %d, Null perms = %d', sum(~isnan(out.elt_rank)), args.n_random)});
  print(gcf, '-dpng', '-r250', fullfile(args.outdir, args.dataname, ...
    sprintf('%s_projection_rank_hist.png', args.dataname)));

  clf; 
  bm = cellfun(@(x) mean(x), out.null_proj);
  scatter(out.elt_proj(f), bm(f), 100, '.');
  hold on; grid on; 
  plot([-1 1], [-1 1], 'r--', 'LineWidth', 2);
  xlim([-0.6, 1]); ylim([-0.4, 0.4]);
  xlabel('On-target projection magnitude');
  ylabel('Mean Null projection magnitude');
  title({sprintf('Scatter of on-target projections vs mean null, %s', args.dataname); ...
      sprintf('N sigs = %d, Null perms = %d', sum(~isnan(out.elt_rank)), args.n_random)});
  print(gcf, '-dpng', '-r250', fullfile(args.outdir, args.dataname, ...
    sprintf('%s_projection_scatter_vsnull.png', args.dataname)));

  clf;
  [zy, zx] = ksdensity(out.elt_zscore(f), 'bandwidth', 0.02);
  plot(zx, zy, 'r', 'LineWidth', 2);
  grid on; xlim([-10 10]);
  xlabel('Z-score of on-target projection relative to null');
  ylabel('Density');
  title({sprintf('PDF of on-target projection z-score (vs null), %s', args.dataname), ...
      sprintf('N sigs = %d, Null perms = %d', sum(~isnan(out.elt_rank)), args.n_random)});
  print(gcf, '-dpng', '-r250', fullfile(args.outdir, args.dataname, ...
    sprintf('%s_proj_zscore_pdf.png', args.dataname)));

  clf;
  scatter(out.elt_proj(f), out.elt_rank(f), 100, '.');
  axis([-1 1 0 1]); grid on;
  xlabel('On-target projection magnitude');
  ylabel('Rank of CGS compared to null');
  title({sprintf('Scatter of projection magnitude vs rank compared to null, %s', args.dataname), ...
      sprintf('N sigs = %d, Null perms = %d', sum(~isnan(out.elt_rank)), args.n_random)});
  print(gcf, '-dpng', '-r250', fullfile(args.outdir, args.dataname, ...
    sprintf('%s_proj_rank_scatter.png', args.dataname)));
end


function [out, brewresids] = projection_iterator(ds, args)

% for now, I assume the variables are strings
grpvar = ds.cdesc(:, ds.cdict(args.groupvars{1}));

for k = 2:numel(args.groupvars)
  grpvar = cellfun(@(x,y) sprintf('%s_%s', x, y), grpvar, ds.cdesc(:, ds.cdict(args.groupvars{k})), ...
      'UniformOutput', 0);
end

[u, c, g] = cellcount(grpvar);
grpmap = containers.Map(u, g);
magns = sqrt(sum(ds.mat.^2));
elt_proj = -2*ones(size(ds.cid));
rand_proj = cell(size(ds.cid));

for ii = 1:numel(g)
  nelts(g{ii},1) = c(ii);
end

roast_mean_proj = -2*ones(size(ds.cid));
roast_resid_proj = -2*ones(size(ds.cid));
roast_resid_corr = -2*ones(size(ds.cid));

roast_null_proj = cell(size(ds.cid));
roast_null_resid_corr = cell(size(ds.cid));
roast_null_resid_proj = cell(size(ds.cid));

% For each element, generate one (or many?) CGSs from its sisters, etc
% Calculate the projection of the element onto the CGS.
brewresids = ds;
brewresids.cid = cellfun(@(x) sprintf('%s_residual', x), brewresids.cid, 'UniformOutput', 0);
brewresids.mat = zeros(size(brewresids.mat));

for k = 1:numel(ds.cid)
    if mod(k, 100) == 0
        disp(k);
    end

    ix = setdiff(grpmap(grpvar{k}), k);

    if ~isempty(ix)
      tcgs = modzs(ds.mat(:, ix), 1:978, ...
          'clip_low_wt', true, ...
          'clip_low_cc', true, ...
          'low_thresh_wt', 0.01, 'low_thresh_cc', 0, ...
          'metric', 'wt_avg');
    
      mag_tcgs = sqrt(sum(tcgs.^2));
      elt_proj(k) = (ds.mat(:, k)' * tcgs)/(magns(k) * mag_tcgs);

      % Grab roast replicates for residuals
      if ~isempty(args.rds)
        repids = find(ismember(args.rds.cid, strsplit(ds.cdesc{k, ds.cdict(args.rcol)}, '|')));
        repds = gctsubset(args.rds, 'csubset', repids);
        mag_rep = sqrt(sum(repds.mat.^2));
        roast_mean_proj(k) = mean((repds.mat' * tcgs)./(mag_rep' * mag_tcgs));

        if size(repds.mat, 2) > 1
          rep_resid = repds.mat - tcgs*(tcgs' * repds.mat)/(mag_tcgs^2);
          rep_resid_corrs = fastcorr(rep_resid, 'type', 'Spearman');
          roast_resid_corr(k) = mean(rep_resid_corrs(triu(true(size(rep_resid_corrs)), 1)));

          rtmean = zeros(size(rep_resid, 2), 1);
          for rk = 1:numel(rtmean)
            rtcgs = modzs(rep_resid(:, setdiff(1:numel(rtmean), rk)), 1:978, ...
                 'clip_low_wt', true, ...
                 'clip_low_cc', true, ...
                 'low_thresh_wt', 0.01, 'low_thresh_cc', 0, 'metric', 'wt_avg');
            rtmean(rk) = (rep_resid(:,rk)' * rtcgs)/(mag_rep(rk) * sqrt(sum(rtcgs.^2)));
                 %(sqrt(sum(rep_resid(:,rk).^2)) * sqrt(sum(rtcgs.^2)));
          end
          roast_resid_proj(k) = mean(rtmean);

          rtcgs = modzs(rep_resid, 1:978, 'clip_low_wt', true, 'clip_low_cc', true, ...
                'low_thresh_wt', 0.01, 'low_thresh_cc', 0, 'metric', 'wt_avg');
          brewresids.mat(:,k) = rtcgs;

        end
      end

      n_random = args.n_random;
      remix = setdiff(1:numel(ds.cid), grpmap(grpvar{k}));
      ncomp = -2*ones(n_random,1);
      nrcomp = -2*ones(n_random,1);
      nrcorrs = -2*ones(n_random,1);
      nrproj = -2*ones(n_random,1);

      for j = 1:n_random
        nix = remix(randperm(numel(remix), numel(ix)));
        ncgs = modzs(ds.mat(:, nix), 1:978, ...
          'clip_low_wt', true, ...
          'clip_low_cc', true, ...
          'low_thresh_wt', 0.01, 'low_thresh_cc', 0, ...
          'metric', 'wt_avg');

        mag_ncgs = sqrt(sum(ncgs.^2));
        ncomp(j) = (ds.mat(:,k)' * ncgs)/(magns(k) * mag_ncgs);

        % Computing the null residuals as below is very expensive and isn't used in subsequent
        % analyses.  Moreover, the residuals when projected onto null CGSs should be quite large;
        % is the magnitude of the projection onto the null CGS which is most interesting

        %if ~isempty(args.rds)
        %  nrcomp(j) = mean((repds.mat' * ncgs)./(mag_rep' * mag_ncgs));
        %  if size(repds.mat, 2) > 1
        %    q_resids = repds.mat - ncgs*(ncgs' * repds.mat)/(mag_ncgs^2);
        %    q_resid_corrs = fastcorr(q_resids, 'type', 'Spearman');
        %    nrcorrs(j) = mean(q_resid_corrs(triu(true(size(q_resid_corrs)), 1)));
        %    nrtmean = zeros(size(q_resids, 2), 1);
        %    for rk = 1:numel(nrtmean)
        %      nrtcgs = modzs(q_resids(:, setdiff(1:numel(nrtmean), rk)), 1:978, ...
        %            'clip_low_wt', true, 'clip_low_cc', true, ...
        %            'low_thresh_wt', 0.01, 'low_thresh_cc', 0, 'metric', 'wt_avg');
        %      nrtmean(rk) = (q_resids(:,rk)' * nrtcgs)/(mag_rep(rk) * sqrt(sum(rtcgs.^2)));
        %    end          
        %    nrproj(j) = mean(nrtmean);
        %  end
        %end 
      end
      
      rand_proj(k) = {ncomp};
      %if ~isempty(args.rds)
      %  roast_null_proj(k) = {nrcomp};
      %  roast_null_resid_corr(k) = {nrcorrs};
      %  roast_null_resid_proj(k) = {nrproj};
      %end
      
    end
end
        
out.cid = ds.cid;
for ii = 1:numel(args.groupvars)
  out.(args.groupvars{ii}) = ds.cdesc(:, ds.cdict(args.groupvars{ii}));
end
out.grpvar = grpvar;
out.nelts = nelts;
out.elt_proj = elt_proj;
out.elt_rank = -100*ones(size(out.elt_proj));
out.elt_zscore = -100*ones(size(out.elt_proj));

for k = 1:numel(out.elt_rank)
  out.elt_rank(k) = sum(out.elt_proj(k) < rand_proj{k})/numel(rand_proj{k});
  out.elt_zscore(k) = (out.elt_proj(k) - mean(rand_proj{k}))/std(rand_proj{k});
end
out.null_proj_mean = cellfun(@(x) mean(x), rand_proj);
out.null_proj_stdev = cellfun(@(x) std(x), rand_proj);

if ~isempty(args.rds)
  out.rep_mean_proj = roast_mean_proj;
  out.rep_resid_proj = roast_resid_proj;
  out.rep_resid_corr = roast_resid_corr;
end

mktbl(fullfile(args.outdir, args.dataname, sprintf('%s_cgs_projection.txt', args.dataname)), ...
    struct_cellarray(out), 'header', fieldnames(out));
out.null_proj = rand_proj;
mkgctx(fullfile(args.outdir, args.dataname, sprintf('%s_proj_residuals.gctx', args.dataname)), ...
    brewresids);

%out.rep_null_proj = roast_null_proj;
%out.rep_null_resid_corr = roast_null_resid_corr;
%out.rep_null_resid_proj = roast_null_resid_proj;

end
