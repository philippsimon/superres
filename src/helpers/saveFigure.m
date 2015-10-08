% Copyright (c) 2015, Philipp Simon Schmidt
% For more details see LICENSE.txt and AUTHORS.txt

function saveFigure(fig, name, aspectRatio)
    if nargin < 3
        aspectRatio = [3 2];
    end
    WIDTH = 1200;
    size = (WIDTH / aspectRatio(1)) * aspectRatio;
    set(fig, 'Position', [0 0 size]);
    set(fig,'PaperPositionMode','auto');
    
    global RESULTS_FOLDER;
    folder = [RESULTS_FOLDER,'/figures'];
    mkdir(folder);
    
    D = dir([folder, '/*.png']);
    figureNum = length(D(not([D.isdir]))) + 1;
    
    fileName = [folder,'/',num2str(figureNum),'_',name];
    format = 'png';
    
    disp(['Write "',fileName,'.',format,'"']);
    print(fig,fileName,['-d',format],'-r0');
end