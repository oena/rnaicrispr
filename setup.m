% Run this matlab code to initialize the environment
% Run each time you are using this code base.
%
% In matlab, after changing to this directory, run
% setup()

% Customize your directory structure to point to the data directory and your desired output directory
prefs.datadir = '../../BroadInstitute/data';
prefs.outdir = '../../BroadInstitute/output';
prefs.annotdir = '../../BroadInstitute/data/rnai_annot';

% Add the code to your matlab path
addpath(genpath('.'));
if ~exist(prefs.outdir, 'dir')
  mkdir(prefs.outdir);
end

prefs.lms = parse_grp(fullfile(prefs.datadir, 'lm_epsilon_n978.grp'));
