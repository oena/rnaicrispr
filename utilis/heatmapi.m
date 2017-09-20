function [h, cmap] = heatmapi(data, varargin)

% h = heatmapi(data, varargin): Generate a heatmap from a data set.
%   data:   a square matrix representing intensity values on an nxn grid
% Parameters:
%   'desc': takes the descriptor value from histogram2 if that was
%       used to generate the data.  Use descriptor either/or xval, yval
%   'interp': Interpolate the heatmap; boolean.  Default is 0.
%   'colormap': colormap to use.  Default is 'cmap', alternatives are 'jet'
%               and 'hot'
%   'map'  : an Nx3 matrix specifying a custom colormap; overwrites colormap 
%       argument
%   'reverse': Invert the y-axis.  Default is 0.  1 forces the y-axis to
%       increase as it goes down.  
%   'colorrange': A vector of two elements indicating the minimum and
%       maximum values for the heat coloring.  Default is [-1 1];
%   'xval': the values of the bounds on the x-coordinates of the data.
%       Numel(xval) must equal size(data,2)+1.  If data is a histogram,
%       xval would correspond to the lower and upper bounds on each bin.
%   'yval': the values of the bounds on the y-coordinates of the data.
%       Numel(yval) must equal size(data,1)+1.
%   'xtick': a vector with the subset of 'xval' with which to annotate the x-axis
%       Defaults to one tick for each column of data when xval and xtick
%       are both unspecified.  (1:size(data,1))+0.5
%   'ytick': a vector with the subset of 'yval' with which to annotate the y-axis
%       Defaults to one tick for each column of data when yval and ytick
%       are both unspecified.  (1:size(data,2))+0.5
%   'xannot': a cell array with labels for xtick
%   'yannot': a cell array with labels for ytick

% valid options for colormap: jet, hot

pnames = {'desc', ...
    'interp', ...
    'colormap', ...
    'colorrange', ...
    'map', ...
    'reverse', ...
    'xval', ...
    'yval', ...
    'xtick', ...
    'ytick', ...
    'xannot', ...
    'yannot'};
dflts = {[], ...
    0, ...
    'cmap', ...
    [-1 1], ...
    [], ...
    0, ...    
    (1:size(data,2)+1)', ...
    (1:size(data,1)+1)', ...
    (1:size(data,2)) + 0.5, ...
    (1:size(data,1)) + 0.5, ...
    {}, ... 
    {}, ... 
    };
args = parse_args(pnames, dflts, varargin{:});

data = horzcat(data, zeros(size(data,1), 1));
data = vertcat(data, zeros(1, size(data, 2)));

if ~isa(data, 'double')
    data = double(data);
end

if ~isempty(args.desc)
    d = args.desc;
    args.xval = (d(1,1):(d(1,2)-d(1,1))/d(1,3):d(1,2))';
    args.yval = (d(2,1):(d(2,2)-d(2,1))/d(2,3):d(2,2))';
end

figure('Position', [500 500 1000 800]);
%pcolor requires args.yval be ordered descending for the graph to have the
%intuitive orientation
args.yval = flipud(args.yval);
h = pcolor(args.xval, args.yval, data);
if args.reverse
    set(gca, 'YDir', 'reverse');
end

% specify the color range
caxis(args.colorrange);

if args.interp
    shading interp;
end

if isempty(args.map)
    if strcmpi(args.colormap, 'cmap')
        r = min(0:1:99, 50)/50;
        b = r(end:-1:1);
        g = min(r, b);
        colormap([r' g' b']);
    else
        colormap(args.colormap);
    end
else
    colormap(args.map);
end

colorbar;

set(gca, 'XTick', args.xtick);
set(gca, 'YTick', args.ytick);
set(gca, 'FontSize', 7);
set(gca, 'XTickLabel', args.xannot);
set(gca, 'XTickLabelRotation', 90);
%xticklabel_rotate([], 90, args.xannot);
set(gca, 'YTickLabel', args.yannot);

end
