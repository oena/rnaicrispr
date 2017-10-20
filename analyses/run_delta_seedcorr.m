function run_delta_seedcorr(varargin)

% For all CGSs in the LINCS shRNA dataset, this function calculates the 
% mean and median correlation to shRNAs that are seedologs of shRNAs in the CGS
% and the mean and median correlation of shRNAs that are seedologs to their 
% seedolog shRNAs.  This is for calculating the change in seed correlation due to 
% the CGS
%
% parameters: 
%    flatpc1: boolean, flatten pc1.  1 by default

params = {'flatpc1'};
dflts = {1};
args = parse_args(params, dflts, varargin{:});

label = ifelse(args.flatpc1, 'pc1=0_', '');

topdir = '/xchip/cogs/projects/rnai_analysis/xpr_analysis2';
outdir = fullfile(topdir, 'shrna/delta_seedcorr');
lms = parse_grp('/xchip/cogs/data/vdb/spaces/lm_epsilon_n978.grp');
datapath = '/xchip/cogs/data/build/a2y13q1/modzs.gctx';

f = dir(fullfile(topdir, 'shrna/annotations/*.txt'));
f = struct2cell(f);
f = f(1,:)';

for k = 1:numel(f)
  out = [];
  disp(k);
  annots = parse_tbl(fullfile(topdir, sprintf('shrna/annotations/%s', f{k})));
  annots.distil_id = cellfun(@(x) strrep(x(4:end-2), ''', u''', '|'), annots.distil_id, 'UniformOutput', 0);

  ds = parse_gctx(datapath, 'rid', lms, 'cid', annots.sig_id);
  if args.flatpc1
    ds = project_pc1(ds, 'use_ref', 1, 'zeromean', 0);
  end

  [ug,cg,gg] = cellcount(annots.pert_iname);
  [us,cs,gs] = cellcount(annots.seed_seq_6mer);
  genemap = containers.Map(ug, 1:numel(ug));
  seedmap = containers.Map(us, gs);
  
  seedolog_mean_corrs = -2*ones(size(annots.sig_id));
  seedolog_median_corrs = -2*ones(size(annots.sig_id));
  cgs_mean_corrs = -2*ones(size(annots.sig_id));
  cgs_median_corrs = -2*ones(size(annots.sig_id));

  sister_count = zeros(size(annots.sig_id));
  seedolog_count = zeros(size(annots.sig_id));

  mat_cgs = [];
  for ii = 1:numel(ug)
    mat_cgs(:,ii) = modzs(ds.mat(:,gg{ii}), 1:978, 'clip_low_wt', true, ...
        'clip_low_cc', true, 'low_thresh_wt', 0.01, 'low_thresh_cc', 0, ...
        'metric', 'wt_avg');
    sister_count(gg{ii}) = cg(ii);
  end

  for ii = 1:numel(gs)
    seedolog_count(gs{ii}) = cs(ii);
  end

  for ii = 1:numel(annots.sig_id)
    seedix = setdiff(seedmap(annots.seed_seq_6mer{ii}), ii);
    cgsix = genemap(annots.pert_iname{ii});
    
    if ~isempty(seedix)
      seedcorrs = fastcorr(ds.mat(:,seedix), ds.mat(:,ii),'type','Spearman');
      cgscorrs = fastcorr(ds.mat(:,seedix), mat_cgs(:,cgsix),'type','Spearman');
      
      seedolog_mean_corrs(ii) = mean(seedcorrs);
      seedolog_median_corrs(ii) = median(seedcorrs);
      cgs_mean_corrs(ii) = mean(cgscorrs);
      cgs_median_corrs(ii) = median(cgscorrs);
    end
  end

  out.sig_id = annots.sig_id;
  out.cell_id = annots.cell_id;
  out.is_gold = double(ismember(annots.is_gold, {'True'}));
  out.pert_iname = annots.pert_iname;
  out.seed_seq_6mer = annots.seed_seq_6mer;
  out.seed_seq_7mer = annots.seed_seq_7mer;
  out.sister_count = sister_count;
  out.seedolog_count = seedolog_count;
  out.seedolog_mean_corrs = seedolog_mean_corrs;
  out.seedolog_median_corrs = seedolog_median_corrs;
  out.cgs_mean_corrs = cgs_mean_corrs;
  out.cgs_median_corrs = cgs_median_corrs;

  mktbl(fullfile(outdir, sprintf('%s_%sdelta_seedcorr.txt', strtok(f{k}, '.'), label)), ...
    struct_cellarray(out), 'header', fieldnames(out));
end


end
