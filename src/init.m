% Copyright (c) 2015, Philipp Simon Schmidt
% For more details see LICENSE.txt and AUTHORS.txt

function init()
    addpath('./helpers', '../lib','../lib/log4m', './test')
    
    global LIB_PATH
    LIB_PATH = fullfile(pwd,'../lib');

    global SAVE_FIGURES;
    if ~exist('SAVE_FIGURES','var') || isempty(SAVE_FIGURES)
        SAVE_FIGURES = 1;
    end

    global SHOW_FIGURES;
    if ~exist('SHOW_FIGURES','var') || isempty(SHOW_FIGURES)
        SHOW_FIGURES = 1;
    end

    global RESULTS_FOLDER;
    RESULTS_FOLDER = 'figures/';
end