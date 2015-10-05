function init()
    addpath('./helpers', '../lib','../lib/log4m')
    
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

    global results_folder;
    results_folder = 'figures/';
end