function fig6a_holdout_plots(prefs, flatpc1)

topdir = fullfile(prefs.tempdir); 	      %%'/cmap/projects/rnai_analysis/xpr_analysis2/output2/';
outdir = fullfile(prefs.outdir, 'figures');   %%'/cmap/users/ismith/projects/seedpaper_misc/figures';


if flatpc1
  pcstr = '_flatpc1';
  hxpr = parse_tbl(fullfile(topdir, 'cgs_pc1=0/xpr_cgs_holdout_stats_cellnull_pc1=0.txt'));
  hsh = parse_tbl(fullfile(topdir, 'cgs_shrna_pilot2_pc1=0/shrna_pilot2_combined_stats_pc1=0.txt'));
  hsh.qval(hsh.empirical_pval == 2) = 1.01;

  %% Comparison Data
  t = parse_tbl(fullfile(topdir, 'cgs_comparison/cgs_query_compare_augmented_xpr2_pc1=0.txt'));

else
  pcstr = '_stddata';

  hxpr = parse_tbl(fullfile(topdir, 'cgs/xpr_cgs_holdout_stats_cellnull_old_variablecgssize.txt'));
  hxpr.qval = fdr_calc(hxpr.pval_cell);
  hsh = parse_tbl(fullfile(topdir, 'cgs_shrna_pilot2/shrna_pilot2_combined_stats_noha1e.txt'));
  hsh.qval(hsh.empirical_pval == 2) = 1.01;

  t = parse_tbl(fullfile(topdir, 'cgs_comparison/cgs_query_compare_augmented_xpr2.txt'));
end

cmap = colormap('bone');  cmap(:,2:3) = 0.4+0.5*cmap(:,2:3); 
cmap(:,1) = 0.9*cmap(:,1);
cmap = vertcat(cmap, repmat([1 1 1],192,1), [0 0 0]);

xprmat = 1.01*ones(numel(unique(setdiff(hxpr.pert_iname, 'EGFP'))), numel(unique(hxpr.cell_id)));
shmat = 1.01*ones(numel(unique(hsh.pert_iname)), numel(unique(hsh.cell_id)));
combmat = 1.01*ones(numel(unique(setdiff(hxpr.pert_iname, 'EGFP'))), numel(unique(t.cell_id)));

cellmap = containers.Map(unique(hxpr.cell_id), 1:numel(unique(hxpr.cell_id)));
cshmap = containers.Map(unique(hsh.cell_id), 1:numel(unique(hsh.cell_id)));
xgmap = containers.Map(setdiff(hxpr.pert_iname, 'EGFP'), 1:numel(setdiff(hxpr.pert_iname, 'EGFP')));
sgmap = containers.Map(unique(hsh.pert_iname), 1:numel(unique(hsh.pert_iname)));

for k = 1:numel(hxpr.pert_iname)
  if ~strcmpi(hxpr.pert_iname{k}, 'EGFP')
    xprmat(xgmap(hxpr.pert_iname{k}), cellmap(hxpr.cell_id{k})) = hxpr.qval(k);
  end
end

for k = 1: numel(hsh.pert_iname)
  shmat(sgmap(hsh.pert_iname{k}), cshmap(hsh.cell_id{k})) = hsh.qval(k);
end


jx = find(sum(shmat > 1, 2) < numel(unique(hsh.cell_id)));
shperts = unique(hsh.pert_iname);
heatmapi(shmat(jx,:), 'map', cmap, 'colorrange', [0 1.01], 'xannot', unique(hsh.cell_id), ...
    'yannot', flipud(shperts(jx)));
set(gca, 'PlotBoxAspectRatioMode', 'manual');
set(gca, 'PlotBoxAspectRatio', [1 3.113 1]);
set(gca, 'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 4 6]);
print(gcf, '-dpdf', '-r250', fullfile(outdir, strcat('fig6a_shrna_holdout', pcstr, '.pdf')));
mktbl(fullfile(outdir, strcat('fig6a_shrna_holdout', pcstr, '_data.txt')), vertcat(horzcat(cell(1), unique(hsh.cell_id)'), horzcat(shperts, num2cell(shmat))));


heatmapi(xprmat, 'map', cmap, 'colorrange', [0 1.01], 'xannot', ...
    unique(hsh.cell_id), 'yannot', flipud(setdiff(unique(hxpr.pert_iname), 'EGFP')));
set(gca, 'PlotBoxAspectRatioMode', 'manual');
set(gca, 'PlotBoxAspectRatio', [1 5 1]);
set(gca, 'FontSize', 8);
set(gcf, 'PaperPosition', [0.25 2.5 3 6]);
print(gcf, '-dpdf', '-r250', fullfile(outdir, strcat('fig6b_crispr_holdout', pcstr, '.pdf')));
mktbl(fullfile(outdir, strcat('fig6b_crispr_holdout', pcstr, '_data.txt')), vertcat(horzcat(cell(1), unique(hxpr.cell_id)'), ...
    horzcat(setdiff(unique(hxpr.pert_iname), 'EGFP'), num2cell(xprmat))));


%% Comparison - data read is now dependent on pc1, see above
%t = parse_tbl(fullfile(topdir, 'cgs_comparison/cgs_query_compare_augmented_xpr2_pc1=0.txt'));
pertmap = containers.Map(unique(t.pert_iname), 1:56);
cellmap = containers.Map(unique(t.cell_id), 1:6);

compmat = 2*ones(56,6);
for ii = 1:numel(t.sig_id)
  compmat(pertmap(t.pert_iname{ii}), cellmap(t.cell_id{ii})) = t.qval(ii);
end
jx = find(sum(compmat > 1, 2) < 6);
compperts = unique(t.pert_iname);
heatmapi(compmat(jx,:), 'map', cmap, 'colorrange', [0 1.01], 'xannot', unique(t.cell_id), ...
    'yannot', flipud(compperts(jx)));
set(gca, 'PlotBoxAspectRatioMode', 'manual');
set(gca, 'PlotBoxAspectRatio', [1 5 1]);
set(gca, 'FontSize', 8);
set(gcf, 'PaperPosition', [0.25 2.5 3 6]);
print(gcf, '-dpdf', '-r250', fullfile(outdir, strcat('fig6c_querycomp', pcstr, '.pdf')));
mktbl(fullfile(outdir, strcat('fig6c_querycomp', pcstr, '_data.txt')), vertcat(horzcat(cell(1), unique(t.cell_id)'), horzcat(compperts(jx), num2cell(compmat(jx,:)))));


% Notables
t.grp = cellfun(@(x,y) sprintf('%s, %s',x,y), t.pert_iname, t.cell_id, 'UniformOutput', 0);

ix = find(and(t.holdout_shrna_qval < 0.25, t.holdout_xpr_qval < 0.25));
[~,jx] = sort(t.grp(ix));
ix = ix(jx);
ixmat = horzcat(t.holdout_shrna_qval(ix), t.holdout_xpr_qval(ix), t.qval(ix));

if size(ixmat,1) > 25
  heatmapi(ixmat(1:25,:), 'map', cmap, 'colorrange', [0 1], 'yannot', flipud(t.grp(ix(1:25))), ...
      'xannot', {'RNAi q-val', 'CRISPR q-val', 'Connectivity q-val'});
  set(gca, 'PlotBoxAspectRatioMode', 'manual');
  set(gca, 'PlotBoxAspectRatio', [1 5 1]);
  set(gcf, 'PaperPosition', [0.25 2.5 4 6]);
  title('Notables fig1');
  print(gcf, '-dpdf', '-r250', fullfile(outdir, strcat('fig6b_holdout_notables1', pcstr, '.pdf')));
  mktbl(fullfile(outdir, strcat('fig6b_holdout_notables1', pcstr, '_data.txt')), vertcat({'', 'RNAi q-val', 'CRISPR q-val', 'Connectivity q-val'}, ...
     horzcat(t.grp(ix(1:25)), num2cell(ixmat(1:25,:)))));

  heatmapi(ixmat(26:end,:), 'map', cmap, 'colorrange', [0 1], 'yannot', flipud(t.grp(ix(26:end))), ...
      'xannot', {'RNAi q-val', 'CRISPR q-val', 'Connectivity q-val'});
  set(gca, 'PlotBoxAspectRatioMode', 'manual');
  set(gca, 'PlotBoxAspectRatio', [1 5 1]);
  set(gcf, 'PaperPosition', [0.25 2.5 4 6]);
  title('Notables fig2');
  print(gcf, '-dpdf', '-r250', fullfile(outdir, strcat('fig6b_holdout_notables2', pcstr, '.pdf')));
  mktbl(fullfile(outdir, strcat('fig6b_holdout_notables2', pcstr, '_data.txt')), vertcat({'', 'RNAi q-val', 'CRISPR q-val', 'Connectivity q-val'}, ...
      horzcat(t.grp(ix(26:end)), num2cell(ixmat(26:end,:)))));
else
  heatmapi(ixmat, 'map', cmap, 'colorrange', [0 1], 'yannot', flipud(t.grp(ix)), ...
      'xannot', {'RNAi q-val', 'CRISPR q-val', 'Connectivity q-val'});
  set(gca, 'PlotBoxAspectRatioMode', 'manual');
  set(gca, 'PlotBoxAspectRatio', [1 5 1]);
  set(gcf, 'PaperPosition', [0.25 2.5 4 6]);
  title('Notables');
  print(gcf, '-dpdf', '-r250', fullfile(outdir, strcat('fig6b_holdout_notables', pcstr, '.pdf')));
  mktbl(fullfile(outdir, strcat('fig6b_holdout_notables', pcstr, '_data.txt')), vertcat({'', 'RNAi q-val', 'CRISPR q-val', 'Connectivity q-val'}, ...
     horzcat(t.grp(ix), num2cell(ixmat))));
end
end
