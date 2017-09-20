function out = fig1b(prefs, varargin)

pnames = {'cls', ...
    'flatpc1', ...
    'cdf', ...
    'outdir', ...
    'datadir', ...
    'annotdir'};
dflts = { {'a375', 'a549','ha1e','hcc515','hepg2','ht29','mcf7','pc3','vcap'}, ...
    0, ...
    0, ...
    fullfile(prefs.outdir, 'figures'), ...
    prefs.datadir, ...
    fullfile(prefs.datadir, 'rnai_annot')};
args = parse_args(pnames, dflts, varargin{:});

label = ifelse(args.flatpc1, 'pc1=0', '');
cdflabel = ifelse(args.cdf, '_cdf', '_pdf');
cls = args.cls;
celllabel = ifelse(numel(cls) == 1, upper(cls{1}), 'all_lincs');

%%topdir = '/cmap/projects/rnai_analysis/xpr_analysis2/';
%%outdir = '/cmap/users/ismith/projects/seedpaper_misc/figures/fig1b_cell';
%%annotdir = fullfile(args.datadir, 'rnai_annot');   %%fullfile(topdir, 'shrna/annotations');
lms = prefs.lms;

if ~exist(args.outdir)
  mkdir(args.outdir);
end

%%datapth = '/cmap/data/build/a2y13q1/modzs.gctx';

% get corrs
sistercorr = [];
seed6corr = [];
seed7corr = [];
allcorrs = [];
pertcount = 0;

% generate meaningful statistics
out.cell_id = [];
out.ss6_fdr25_frac = [];
out.ss7_fdr25_frac = [];
out.gene_fdr25_frac = [];
out.ss6_mean = [];  out.ss6_std = [];
out.ss7_mean = [];  out.ss7_std = [];
out.gene_mean = []; out.gene_std = [];

for ii = 1:numel(cls)
  cls_6c = [];
  cls_7c = [];
  cls_genec = [];
  cls_allc = [];

  annots = [];  
  t = get_data(args, 'rnai_annot', 'cell_id', cls{ii});   %%parse_tbl(fullfile(annotdir, sprintf('%s_hp.txt', cls{ii})));
  annots = struct_append(annots, t);
  pertcount = pertcount + numel(annots.pert_iname);

  ds = get_data(args, 'rnai_lvl5', 'rid', lms, 'cid', annots.sig_id);   %%parse_gctx(datapth, 'rid', lms, 'cid', annots.sig_id);

  if args.flatpc1
    ds = project_pc1(ds, args);
  end

  [u,c,g] = cellcount(annots.pert_iname);
  for k = 1:numel(u)
    if c(k) > 1
      cls_genec = vertcat(cls_genec, tri2vec(fastcorr(ds.mat(:, g{k}), 'type', 'Spearman')));
    end
  end

  [u6,c6,g6] = cellcount(annots.seed_seq_6mer);
  for k = 1:numel(u6)
    if c6(k) > 1
      cls_6c = vertcat(cls_6c, tri2vec(fastcorr(ds.mat(:, g6{k}), 'type', 'Spearman')));
    end
  end

  [u7,c7,g7] = cellcount(annots.seed_seq_7mer);
  for k = 1:numel(u7)
    if c7(k) > 1
      cls_7c = vertcat(cls_7c, tri2vec(fastcorr(ds.mat(:, g7{k}), 'type', 'Spearman')));
    end
  end

  cls_allc = tri2vec(fastcorr(ds.mat, 'type', 'Spearman'));

  % Summarize cell line
  % Using tenth of percentiles as below is faster than ranking relative to all values and precise to 0.001. 
  q = quantile(cls_allc, 1000);
  
  out.cell_id = vertcat(out.cell_id, cls(ii)); 
  out.gene_mean = vertcat(out.gene_mean, mean(cls_genec));  
  out.gene_std = vertcat(out.gene_std, std(cls_genec));
  out.ss6_mean = vertcat(out.ss6_mean, mean(cls_6c));
  out.ss6_std = vertcat(out.ss6_std, std(cls_6c));
  out.ss7_mean = vertcat(out.ss7_mean, mean(cls_7c));
  out.ss7_std = vertcat(out.ss7_std, std(cls_7c));
  
  out.gene_fdr25_frac = vertcat(out.gene_fdr25_frac, mean(fdr_calc(arrayfun(@(x) mean(x < q), cls_genec)) < 0.25));
  out.ss6_fdr25_frac = vertcat(out.ss6_fdr25_frac, mean(fdr_calc(arrayfun(@(x) mean(x < q), cls_6c)) < 0.25));
  out.ss7_fdr25_frac = vertcat(out.ss7_fdr25_frac, mean(fdr_calc(arrayfun(@(x) mean(x < q), cls_7c)) < 0.25));

  % Append data
  sistercorr = vertcat(sistercorr, cls_genec);
  seed6corr = vertcat(seed6corr, cls_6c);
  seed7corr = vertcat(seed7corr, cls_7c);
  allcorrs = vertcat(allcorrs, cls_allc);
end

if strcmpi(celllabel, 'all_lincs')
  mktbl(fullfile(args.outdir, sprintf('fig1b_stats_%s%s.txt', celllabel, label)), struct_cellarray(out), ...
    'header', fieldnames(out));
end

%[h,p,ksstat] = kstest2(sistercorr, seed6corr);
%mykstat = ksstat*sqrt(numel(sistercorr)*numel(seed6corr)/(numel(sistercorr)+numel(seed6corr)));
%critval = sqrt(-0.5 * log(myalpha/2));
%http://www.mathematik.uni-kl.de/~schwaar/Exercises/Tabellen/table_kolmogorov.pdf
%disp('KS Test parameters available in the code, including critical values');


% Plot
figure('Position', [10 10 1200 1000]); hold on; %grid on;
[b,a] = ksdensity(sistercorr, 'bandwidth', 0.02, 'function', ifelse(args.cdf, 'cdf', 'pdf'));
plot(a,b,'r','LineWidth',2);
[b,a] = ksdensity(seed6corr, 'bandwidth', 0.02, 'function', ifelse(args.cdf, 'cdf', 'pdf'));
plot(a,b,'b','LineWidth',2);
[b,a] = ksdensity(seed7corr, 'bandwidth', 0.02, 'function', ifelse(args.cdf, 'cdf', 'pdf'));
plot(a,b,'LineWidth',2,'Color',[0 0.5 0]);
[b,a] = ksdensity(allcorrs(randperm(numel(allcorrs), min(numel(allcorrs), 1e6))), ...
        'bandwidth', 0.01, 'function', ifelse(args.cdf, 'cdf', 'pdf'));
plot(a,b,'k','LineWidth',2);
xlabel('Spearman Correlation');
ylabel(ifelse(args.cdf, 'Cumulative Density', 'Density'));
xlim([-0.8, 0.8]);
legend('Same Gene', 'Same 6-mer Seed','Same 7-mer Seed', 'All pairs');
title({sprintf('Distribution of all pairwise correlations for shRNA signatures in %s', celllabel); ...
    sprintf('N = %d, %s', pertcount, label)});
print(gcf, '-dpdf', '-r250', '-bestfit', fullfile(args.outdir, sprintf('fig1b_%s_%s%s.pdf', celllabel, label, cdflabel)));

end
