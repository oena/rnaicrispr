function [ds] = get_data(prefs, dataset, varargin)

% This function is used by a majority of the RNAi-CRISPR comparison code to 
% locate and parse data correctly.  If you have downloaded the data and code, 
% you should only need to modify the preferences.txt file indicating where the
% data is located and where the analyses should be output. 
%
% ds = get_data(prefs, dataset, varargin)
%    prefs - struct from setup.m that includes datadir, annotdir
%    dataset - one of {'rnai_lvl5', 'rnai_annot', 'refpc', 'refmean'}
%    varargin - 'cid', 'rid', 'flatpc1', 'cell_id'


pnames = {'cid', ...
    'rid', ...
    'flatpc1', ...
    'cell_id'};
dflts = {{}, ...
    {}, ...
    0, ...
    ''};
args = parse_args(pnames, dflts, varargin{:});


switch dataset
  case 'rnai_lvl5'
    if args.flatpc1
      ds = parse_gctx(fullfile(prefs.datadir, 'rnai', 'a2_shrna_lvl5_flatpc1_extract_lm_n116782x978.gctx'), 'rid', args.rid, 'cid', args.cid);
    else
      ds = parse_gctx(fullfile(prefs.datadir, 'rnai', 'a2_shrna_lvl5_extract_lm_n116782x978.gctx'), 'rid', args.rid, 'cid', args.cid);
    end
  case 'rnai_annot'
    ds = parse_tbl(fullfile(prefs.annotdir, sprintf('%s_hp.txt', args.cell_id)));
  case 'refpc'
    ds = parse_gctx(fullfile(prefs.datadir, 'pca', 'affogato_pcs_global_n978x978.gctx'), 'rid', args.rid, 'cid', args.cid);
  case 'refmean'
    ds = parse_gctx(fullfile(prefs.datadir, 'pca', 'affogato_means_global_n3x978.gctx'), 'rid', args.rid, 'cid', args.cid);
end

end
