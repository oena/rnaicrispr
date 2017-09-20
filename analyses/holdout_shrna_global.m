function [ret, nulls] = holdout_shrna_global(prefs, varargin)

pnames = {'sig_id', ...
    'pert_iname', ...
    'cell_ids', ...
    'groupvar', ...
    'outdir', ...
    'mkfigs', ...
    'dataname', ...
    'null_iters', ...
    'flatpc1', ...
    'datadir', ...
    'lms'};
dflts = {'', ...
    '', ...
    {'A375','A549','HA1E','HCC515','HEPG2','HT29','MCF7','PC3','VCAP'}, ...
    '', ...
    prefs.outdir, ... 
    0, ...
    'myholdout', ...
    1000, ...
    1, ...
    prefs.datadir, ...
    prefs.lms};

args = parse_args(pnames, dflts, varargin{:});
if ~exist(args.outdir, 'dir')
  mkdir(args.outdir);
end
if ~exist(fullfile(args.outdir, 'bycontext'), 'dir')
  mkdir(args.outdir);
end

%%datadir = prefs.datadir;
%%lms = prefs.lms;

% Make unique by cell line
if and(isempty(args.sig_id), and(numel(args.cell_ids) > 1, ~isempty(args.pert_iname)))
  ret = [];
  nulls = [];
  for k = 1:numel(args.cell_ids)
    [r,n] = holdout_shrna_global('pert_iname', args.pert_iname, 'cell_ids', args.cell_ids(k), ...
        'groupvar', args.groupvar, 'outdir', args.outdir, 'mkfigs', args.mkfigs, ...
        'dataname', sprintf('%s_%s', args.dataname, args.cell_ids{k}), ...
        'null_iters', args.null_iters, 'flatpc1', args.flatpc1);
    ret = struct_append(ret, r);
    nulls = struct_append(nulls, n);
  end
  
  return;
end

%%if and(isempty(args.sig_id), ~isempty(args.pert_iname))
%%  % I don't really want to do this
%%  % Get the sigs corresponding to the given pert_inames in the LINCS lines
%%  sigs = sig_info(sprintf('{"pert_type":"trt_sh", "cell_id":{$in:[%s]}, "pert_itime":"96 h"%s}', ...
%%      cell2str(args.cell_ids, 'quotes', '"'), ...
%%      sprintf(', "pert_iname":{$in:[%s]}', cell2str(args.pert_iname, 'quotes', '"'))));
%%  f = fieldnames(sigs);
%%  sc = struct2cell(sigs)';
%%else
%%  sigs = sig_info(args.sig_id);
%%  f = fieldnames(sigs);
%%  sc = struct2cell(sigs)';
%%end

sigs = parse_tbl(fullfile(args.datadir, 'rnai_annot', sprintf('%s_hp.txt', args.cell_ids{1})));
sc = struct_cellarray(sigs);
f = fieldnames(sigs);

fpertid = cellstrfind(f, 'pert_id');
fcellid = cellstrfind(f, 'cell_id');
%%[~, ix] = unique(cellfun(@(x,y) horzcat(x,y), sc(:,fpertid), sc(:,fcellid), 'UniformOutput', 0));
%%sc = sc(ix,:);

%%lms = parse_grp('/xchip/cogs/data/vdb/spaces/lm_epsilon_n978.grp');
%%ds = parse_gctx('/xchip/cogs/data/build/a2y13q1/modzs.gctx', 'cid', sc(:,cellstrfind(f, 'sig_id')), ...
%%        'rid', lms);

% A bit unnecessary, but it's nice to be explicit
if args.flatpc1
  ds = get_data(args, 'rnai_lvl5', 'cid', sc(:, cellstrfind(f, 'sig_id')), 'rid', args.lms, 'flatpc1', 1);
else
  ds = get_data(args, 'rnai_lvl5', 'cid', sc(:, cellstrfind(f, 'sig_id')), 'rid', args.lms, 'flatpc1', 0);
end


%%if args.flatpc1
%%  [coeff, score, latent, tsq, expl] = pca(ds.mat');
%%  score(:,1) = 0;
%%  mumat = mean(ds.mat, 2);
%%  mymat = (score*coeff' + repmat(mumat', size(ds.mat, 2), 1))';
%%  ds.mat = mymat;
%%end


idx = match_vectors(ds.cid, sc(:, cellstrfind(f, 'sig_id')));
sc = sc(idx,:);

ds.cdesc = sc;
ds.chd = f;

nulls = compute_null_holdout(ds, 30, 3:6, args);

ret = get_holdout_data(ds, nulls, args);

mktbl(fullfile(args.outdir, sprintf('%s_stats%s.txt', args.dataname, ifelse(args.flatpc1,'_pc1=0',''))), ...
    ret(2:end, 1:end-2), 'header', ret(1, 1:end-2));

ret = cell2structarray(ret(2:end,:), ret(1,:)');

end


function ret = get_holdout_data(ds, nulls, args)

%%lms = parse_grp('/xchip/cogs/data/vdb/spaces/lm_epsilon_n978.grp');
ret = {};
retfields = {'pert_iname', 'cell_id', 'n_perts', 'groupvar', 'time_pt', ...
    'holdout_median', 'null_median', 'holdout_std', 'null_std', ...
    'empirical_pval', 'qval', 'delta', 'tdelta', 'tmeddelta', ...
    'holdout_data', 'null_data'};

for k = 3:6
  [t.(sprintf('bn%d',k)), t.(sprintf('an%d',k))] = ksdensity(nulls.(sprintf('n%d',k)), 'bandwidth', 0.01);
end

% Generate null distribution with same sizes of CGSs
%nd.ds3 = parse_gctx('../data/cgs/nullshrna_cgs3_randperm_bycontext_ngenes=3945_n66962x978.gctx', 'rid', lms);
%nd.ds4 = parse_gctx('../data/cgs/nullshrna_cgs4_randperm_bycontext_ngenes=3945_n50218x978.gctx', 'rid', lms);
%nd.ds5 = parse_gctx('../data/cgs/nullshrna_cgs5_randperm_bycontext_ngenes=3945_n40176x978.gctx', 'rid', lms);
%nd.ds6 = parse_gctx('../data/cgs/nullshrna_cgs6_randperm_bycontext_ngenes=3945_n33474x978.gctx', 'rid', lms);

%nd.ds3 = gctsubset(nd.ds3, 'csubset', randperm(66962, 2000));
%nd.ds4 = gctsubset(nd.ds4, 'csubset', randperm(50218, 2000));
%nd.ds5 = gctsubset(nd.ds5, 'csubset', randperm(40176, 2000));
%nd.ds6 = gctsubset(nd.ds6, 'csubset', randperm(33474, 2000));

col_piname = cellstrfind(ds.chd, 'pert_iname');
col_pid = cellstrfind(ds.chd, 'pert_id');

% generate cgs
[u_cell, ~, g_cell] = cellcount(ds.cdesc(:, cellstrfind(ds.chd, 'cell_id')));
if args.mkfigs
  figure('Position', [100 100 1200 1000]);
end

for k_cell = 1:numel(u_cell)
    cell_ix = g_cell{k_cell};
    disp(u_cell{k_cell});
    sub = gctsubset(ds, 'csubset', cell_ix);

    [u, ~, g] = cellcount(sub.cdesc(:, col_piname));
    groupvar = u;

    for k_pert = 1:numel(u)
        pertcount = numel(unique(sub.cdesc(g{k_pert}, col_pid)));
        if mod(k_pert, 1000) == 0
           disp(k_pert);
        end

        if pertcount >= 6
          cgssize = min(6, floor(pertcount/2));
          sims = [];
          pertsigs = [];
          for ii = 1:30
            idx = randperm(numel(g{k_pert}), 2*cgssize);
            seta = idx(1:cgssize); setb = idx(cgssize+1:end);

            cgs_a = modzs(sub.mat(:, g{k_pert}(seta)), 1:978, ...
                'clip_low_wt', true, ...
                'clip_low_cc', true, ...
                'low_thresh_wt', 0.01, 'low_thresh_cc', 0, ...
                'metric', 'wt_avg');
            cgs_b = modzs(sub.mat(:, g{k_pert}(setb)), 1:978, ...
                'clip_low_wt', true, ...
                'clip_low_cc', true, ...
                'low_thresh_wt', 0.01, 'low_thresh_cc', 0, ...
                'metric', 'wt_avg');

            sims(ii) = fastcorr(cgs_a, cgs_b, 'type', 'Spearman');
            pertsigs = horzcat(pertsigs, cgs_a, cgs_b);
          end
          
          %nullsims = fastcorr(nd.(sprintf('ds%d', cgssize)).mat, pertsigs, 'type', 'Spearman');
          nullsims = nulls.(sprintf('n%d', cgssize));
          
          % Obsolete Wilcoxon method
          %[p, h, stats] = ranksum(sims, mat2vec(nullsims), 'tail', 'right');

          % Empirical null - permute within each column to generalize, to avoid comparing each 
          % individual null CGS with the N holdout CGSs
          empp = (sum(nullsims > median(sims(:))) + 1)/(numel(nullsims)+1);
          delta = (median(sims(:)) - median(nullsims(:)))/sqrt(std(mat2vec(sims))^2 + std(mat2vec(nullsims))^2);
          tmeddelta = (median(sims(:)) - median(nullsims(:)))/sqrt(std(mat2vec(sims))^2/numel(sims) ...
                 + std(mat2vec(nullsims))^2/numel(nullsims));
          tdelta = (mean(sims(:)) - mean(nullsims(:)))/sqrt(std(mat2vec(sims))^2/numel(sims) ...
                 + std(mat2vec(nullsims))^2/numel(nullsims));
          

          mytimept = sub.cdesc{g{k_pert}(1), cellstrfind(ds.chd, 'pert_itime')};
          ret = vertcat(ret, {u{k_pert}, u_cell{k_cell}, pertcount, u{k_pert}, mytimept, ...
                  median(sims(:)), median(nullsims(:)), std(mat2vec(sims)), std(mat2vec(nullsims)), ...
                  empp, 1, delta, tdelta, tmeddelta, ...
                  mat2vec(sims), nullsims});

          if args.mkfigs
            clf; hold on; grid on;
            [b,a] = ksdensity(sims(:), 'bandwidth', 0.01); plot(a,b, 'r', 'LineWidth', 2);
            [b,a] = ksdensity(mat2vec(nullsims), 'bandwidth', 0.01); plot(a,b, 'b', 'LineWidth', 2);
            switch cgssize
              case 3; plot(t.an3, t.bn3, 'g', 'LineWidth', 2);
              case 4; plot(t.an4, t.bn4, 'g', 'LineWidth', 2);
              case 5; plot(t.an5, t.bn5, 'g', 'LineWidth', 2);
              case 6; plot(t.an6, t.bn6, 'g', 'LineWidth', 2);
            end
            xlim([-1, 1]);
            xlabel('shRNA CGS Spearman Similarity'); ylabel('Density');
            legend('shRNA CGS - Holdout', 'shRNA CGS - Null', 'Null - Null', 'Location', 'NorthWest');
            title({sprintf('shRNA CGS Similarity for %s in %s at %s H', u{k_pert}, u_cell{k_cell}, mytimept); ...
                sprintf('Metric: Spearman, median holdout sim: %0.3f, median holdout-null: %0.3f', ...
                    median(sims(:)), median(mat2vec(nullsims)))});
            print(gcf, '-dpng', '-r250', fullfile(args.outdir, 'bycontext', ... 
                   sprintf('%s_%s_%dH_shrnacgssim.png', u{k_pert}, u_cell{k_cell}, 96)));
          end
        else
          ret = vertcat(ret, {u{k_pert}, u_cell{k_cell}, pertcount, u{k_pert}, ...
                  sub.cdesc{g{k_pert}(1), cellstrfind(ds.chd, 'pert_itime')}, ...
                  -1, -1, -1, -1, ...
                  2, 1, 0, 0, 0, [], []});
        end
    end
end

pvals = cell2mat(ret(:, cellstrfind(retfields, 'empirical_pval')));
qvals = ones(size(pvals));
qvals(pvals < 2) = fdr_calc(pvals(pvals < 2));

ret(:, cellstrfind(retfields, 'qval')) = num2cell(qvals);
ret = vertcat(retfields, ret);

end


function nulls = compute_null_holdout(ds, n_perms, cgs_size, args)
  nulls = [];
  niter = args.null_iters; %max(10, args.null_iters);

  for ii = 1:numel(cgs_size)
    csz = cgs_size(ii);
    nullsims = zeros(niter,1);

    for k = 1:niter
      if mod(k, 200) == 0
        disp(k);
      end
      mycids = randperm(numel(ds.cid), 2*csz);
      sims = zeros(n_perms, 1);
     
      for iter = 1:n_perms
        seta = randperm(numel(mycids), csz);
        setb = setdiff(1:numel(mycids), seta);

        cgs_a = modzs(ds.mat(:, mycids(seta)), 1:978, ...
          'clip_low_wt', true, 'clip_low_cc', true, ...
          'low_thresh_wt', 0.01, 'low_thresh_cc', 0, 'metric', 'wt_avg');
        cgs_b = modzs(ds.mat(:, mycids(setb)), 1:978, ...
          'clip_low_wt', true, 'clip_low_cc', true, ...
          'low_thresh_wt', 0.01, 'low_thresh_cc', 0, 'metric', 'wt_avg');

        sims(iter) = fastcorr(cgs_a, cgs_b, 'type', 'Spearman');
      end
      nullsims(k) = median(sims);
    end
    nulls.(sprintf('n%d', csz)) = nullsims;
  end

%mktbl(fullfile(args.outdir, sprintf('%s_null_holdoutsims.txt', args.dataname)), struct_cellarray(nulls), ...
%    'header', fieldnames(nulls));
end
