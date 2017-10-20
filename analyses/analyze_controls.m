function retds = analyze_controls()

% Examine the replicate correlations of empty_t, untrt, empty_vector, and EGFP controls of size n = 4
% Create controls of size four from roast data, examine correlations
% Also consider the distribution of all landmark z-scores for brew controls vs BRAF

args.topdir = '/xchip/cogs/projects/rnai_analysis/xpr_analysis2/';
args.outdir = '/xchip/cogs/projects/rnai_analysis/xpr_analysis2/output/controls';
brewds = parse_gctx('../data/XPR.XPR001_Aggregate_COMPZ.MODZ_SCORE_LM_n2715x978.gctx');
roastds = parse_gctx('../data/XPR.XPR001_Aggregate_ZSPCQNORM_n11144x978.gctx');


retds = generate_controls(roastds, args);
mk_control_figs(retds, args);
%analyze_brew_controls(brewds, args);
%analyze_roast_controls(roastds, args);

end


function analyze_roast_controls(roastds, args)

% plot rms (i.e. length) of each signature vector
absmag = sqrt(sum(roastds.mat.^2)/978);
[u,~,g] = cellcount(roastds.cdesc(:, roastds.cdict('pert_mfc_desc')));
[~,ix] = sort(cellfun(@(x) median(absmag(x)), g));
boxplot(absmag, roastds.cdesc(:, roastds.cdict('pert_mfc_desc')), ...
    'plotstyle', 'compact', 'symbol', '.', 'grouporder', u(ix));
grid on; ylim([0 3]);
ylabel('RMS Roast Z-score Signature');
title({'Distribution of XPR roast RMS signature z-scores from XPR pilot by target gene'; ...
    sprintf('N signatures: %d; %d pert_inames; %d contexts', numel(absmag), numel(u), 15)});
print(gcf, '-dpng', '-r250', fullfile(args.outdir, 'xpr_roast_control_rms_boxplot_bygene.png'));

end


function mk_control_figs(ds, args)

[u,~,g] = cellcount(ds.cdesc(:,8));
corrs = cell2mat(ds.cdesc(:,2));

figure('Position', [100 100 1000 800]); hold on; grid on;
[b,a] = ksdensity(corrs(horzcat(g{setdiff(1:20, [7,9,20])})), 'bandwidth', 0.03);
plot(a,b, 'r', 'LineWidth', 2);
[b,a] = ksdensity(corrs(g{7}), 'bandwidth', 0.03);
plot(a,b, 'g', 'LineWidth', 2);
[b,a] = ksdensity(corrs(g{9}), 'bandwidth', 0.03);
plot(a,b, 'b', 'LineWidth', 2);
[b,a] = ksdensity(corrs(g{20}), 'bandwidth', 0.03);
%plot(a,b, 'k', 'LineWidth', 2);
xlabel('Distil_cc_q75, distil_nsample = 4');
ylabel('Density');
%legend('All XPRs', u{7}, u{9}, u{20}, 'Location', 'NorthWest');
legend('All XPRs', u{7}, u{9}, 'Location', 'NorthWest');
xlim([-0.8, 0.8]);
title({sprintf('Distribution of XPR Roast Replicate Correlations for Brews of size <= 4'); ...
   sprintf('All contexts, n = %d', size(ds.mat, 2))});
print(gcf, '-dpng', '-r250', fullfile(args.outdir, 'brew_n=4', 'xpr_ccq75_brewn=4_densityplot.png'));


figure('Position', [100 100 1000 800]); hold on; grid on;
[b,a] = ksdensity(vertcat(ds.cdesc{horzcat(g{setdiff(1:20, [7,9,20])}),1}), 'bandwidth', 0.03);
plot(a,b, 'r', 'LineWidth', 2);
[b,a] = ksdensity(vertcat(ds.cdesc{horzcat(g{7}),1}), 'bandwidth', 0.03);
plot(a,b, 'g', 'LineWidth', 2);
[b,a] = ksdensity(vertcat(ds.cdesc{horzcat(g{9}),1}), 'bandwidth', 0.03);
plot(a,b, 'b', 'LineWidth', 2);
[b,a] = ksdensity(vertcat(ds.cdesc{horzcat(g{20}),1}), 'bandwidth', 0.03);
%plot(a,b, 'k', 'LineWidth', 2);
xlabel('All Spearman correlations, distil_nsample = 4');
ylabel('Density');
%legend('All XPRs', u{7}, u{9}, u{20}, 'Location', 'NorthWest');
legend('All XPRs', u{7}, u{9}, 'Location', 'NorthWest');
xlim([-0.8, 0.8]);
title({sprintf('Distribution of XPR Roast Replicate Correlations for Brews of size <= 4'); ...
   sprintf('All contexts, n = %d', size(ds.mat, 2))});
print(gcf, '-dpng', '-r250', fullfile(args.outdir, 'brew_n=4', 'xpr_allcc_brewn=4_densityplot.png'));

% Confine analysis to A375, 96H
uix = intersect(cellstrfind(ds.cdesc(:,4), 'A375'), find(cell2mat(ds.cdesc(:,9)) == 96));
sub = gctsubset(ds, 'csubset', uix);
ix_egfp = cellstrfind(sub.cdesc(:,8), 'EGFP');
perts_egfp = unique(sub.cdesc(ix_egfp, 6));

mydata = [];
labels = cell(0,1);

for k = 1:numel(perts_egfp)
  ix = cellstrfind(sub.cdesc(:,6), perts_egfp{k});
  mydata = vertcat(mydata, mat2vec(sub.mat(:,ix)));
  labels = vertcat(labels, repmat({sprintf('EGFP_%d', k)}, numel(sub.mat(:,ix)),1));
end

emptyv_data = mat2vec(sub.mat(:, cellstrfind(sub.cdesc(:,8), 'EMPTY_VECTOR')));
mydata = vertcat(mydata, emptyv_data);
labels = vertcat(labels, repmat({'EmptyVector'}, numel(emptyv_data),1));

braf_data = mat2vec(sub.mat(:, cellstrfind(sub.cdesc(:,8), 'BRAF')));
mydata = vertcat(mydata, braf_data);
labels = vertcat(labels, repmat({'BRAF'}, numel(braf_data),1));

figure('Position', [10 10 1000 800]);
boxplot(mydata, labels, 'jitter', 0.10, 'symbol', 'b.');
set(findobj(gca, 'type', 'line'), 'linew', 2);
h = findobj(gca, 'Type', 'text');
grid on;

for ii = 1:numel(h)
  set(h(ii), 'Rotation', 90, 'HorizontalAlignment', 'right', 'FontSize', 10);
end
set(gca, 'OuterPosition', [0.02 0.10 1 0.85]);
ylim([-8, 8]);
ylabel('Z-score distribution');
title({'Distribution of Z-scores from control XPR signatures, A375 96H'; ...
    sprintf('Brews of size 4, n = %d', size(sub.mat, 2))});
print(gcf, '-dpng', '-r250', fullfile(args.outdir, 'brew_n=4', 'xpr_control_zscore_boxplot_A375_brew_n=4.png'));

% Density plot
figure('Position', [10 10 1000 800]); hold on; grid on;
mycolmap = varycolor(13);
for k = 1:numel(perts_egfp)
    ix = cellstrfind(sub.cdesc(:,6), perts_egfp{k});
    [b,a] = ksdensity(mat2vec(sub.mat(:,ix)), 'function', 'cdf', 'bandwidth', 0.1);
    plot(a,b, ':', 'Color', mycolmap(k,:), 'LineWidth', 2);
end
[b,a] = ksdensity(emptyv_data, 'function', 'cdf', 'bandwidth', 0.1);
plot(a,b, 'Color', mycolmap(12,:), 'LineWidth', 2);
[b,a] = ksdensity(braf_data, 'function', 'cdf', 'bandwidth', 0.1);
plot(a,b, 'Color', mycolmap(13,:), 'LineWidth', 2);
legend(vertcat(perts_egfp, 'Empty_t', 'BRAF'), 'Location', 'SouthEast');
xlabel('Z-score, all landmarks');
ylabel('Cumulative Density');
xlim([-4 4]);
title('Cumulative density plot of z-score distributions for XPR control signatures in A375');
print(gcf, '-dpng', '-r250', fullfile(args.outdir, 'brew_n=4', 'xpr_control_zscore_cdf_A375_brew_n=4.png'));


end




function retds = generate_controls(ds, args)

grpvar = cellfun(@(x,y,z) sprintf('%s_%d_%s',x,y,z), ds.cdesc(:,ds.cdict('cell_id')), ds.cdesc(:,ds.cdict('pert_time')), ds.cdesc(:, ds.cdict('pert_id')),'UniformOutput', 0);

ds.mat = clip(ds.mat, -10, 10);
[u, c, g] = cellcount(grpvar);

cix = [4, 8, 13, 16, 17, 19, 21];
outmat = [];
outcid = {};
outchd = {'allcorrs', 'cc_q75', 'distil_nsample', ds.chd{cix}, 'roast_cix'};
outcell = {};

for k = 1:numel(u)
    b = make_random_groups(g{k}, 4);
    
    if or(numel(b) == 1, numel(b{end}) == 4)
        b = vertcat(b,{0});
    end

    for ii = 1:numel(b)-1
        cmat = fastcorr(ds.mat(:, b{ii}), 'type', 'Spearman');
        cmat = cmat(triu(true(size(cmat)),1));
        t = horzcat({cmat}, {quantile(cmat, 0.75)}, {numel(b{ii})});
        tcdesc = ds.cdesc(b{ii}(1), cix);

        outmat = horzcat(outmat, ...
            modzs(ds.mat(:, b{ii}), 1:numel(b{ii}), 'clip_low_wt', true, 'clip_low_cc', true));
        outcid = vertcat(outcid, sprintf('%s_%d_n=%d', u{k}, ii, numel(b{ii})));
        outcell = vertcat(outcell, horzcat(t, tcdesc, b(ii)));
    end
end

qscore = parse_tbl(fullfile(args.topdir, 'data', 'cmap_lims_expandedSpacers.txt'));
qmap = containers.Map(qscore.pert_id, qscore.score);
scoredata = -666*ones(size(outcid));
for k = 1:numel(outcid)
    if isKey(qmap, outcell{k,6})  % Good Lord, this code is brittle.  Need cix = 13.
        scoredata(k) = qmap(outcell{k,6});
    end
end
outcell = horzcat(outcell, num2cell(scoredata));
outchd = horzcat(outchd, 'xpr_quality_score');
retds = mkgctstruct(outmat, ...
    'rid', ds.rid, 'rhd', ds.rhd, 'rdesc', ds.rdesc, ...
    'cid', outcid, 'chd', outchd, 'cdesc', outcell);

end


function analyze_brew_controls(brewds, args)
%outdir = '/xchip/cogs/projects/rnai_analysis/xpr_analysis2/output/controls';
%ds = parse_gctx('../data/XPR.XPR001_Aggregate_COMPZ.MODZ_SCORE_LM_n2715x978.gctx');
ds = brewds;
outdir = args.outdir;

ix_egfp = cellstrfind(ds.cdesc(:, ds.cdict('pert_mfc_desc')), 'EGFP');
perts_egfp = unique(ds.cdesc(ix_egfp, ds.cdict('pert_id')));
[u_pert, ~, g_pert] = cellcount(ds.cdesc(:, ds.cdict('pert_id')));

ix = find(ismember(u_pert, perts_egfp));
mydata = [];
labels = cell(0,1);

for k = 1:numel(ix)
  mydata = vertcat(mydata, mat2vec(ds.mat(:, g_pert{ix(k)})));
  labels = vertcat(labels, repmat({sprintf('EGFP_%d', k)}, numel(ds.mat(:, g_pert{ix(k)})),1));
end

untrt_data = mat2vec(ds.mat(:, g_pert{cellstrfind(u_pert, 'CMAP-000')}));
mydata = vertcat(mydata, untrt_data);
labels = vertcat(labels, repmat({'UnTrt'}, numel(untrt_data), 1));

emptyv_data = mat2vec(ds.mat(:, g_pert{cellstrfind(u_pert, 'TRCN0000208001')}));
mydata = vertcat(mydata, emptyv_data);
labels = vertcat(labels, repmat({'EmptyVector'}, numel(emptyv_data),1));

% add pos con for reference
ix_braf = cellstrfind(ds.cdesc(:, ds.cdict('pert_mfc_desc')), 'BRAF');
braf_data = mat2vec(ds.mat(:, ix_braf));
mydata = vertcat(mydata, braf_data);
labels = vertcat(labels, repmat({'BRAF XPRs'}, numel(braf_data),1));

figure('Position', [10 10 1200 1000]);
boxplot(mydata, labels, 'jitter', 0.10, 'symbol', 'b.');
set(findobj(gca, 'type', 'line'), 'linew', 2);
h = findobj(gca, 'Type', 'text');
grid on;

for ii = 1:numel(h)
  set(h(ii), 'Rotation', 90, 'HorizontalAlignment', 'right', 'FontSize', 10);
end
set(gca, 'OuterPosition', [0.02 0.10 1 0.85]);

ylim([-5 5]);
ylabel('Z-score distribution');
title('Distribution of Z-scores from control XPR signatures aggregated across 5 cell lines, 3 time pts');
print(gcf, '-dpng', '-r250', fullfile(outdir, 'xpr_control_zscore_boxplot.png'));

end
