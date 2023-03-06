function pds =  ignit()
% initilize the all necessary 

tmpPath = which('ignit');
[phd, ~, ~] = fileparts(tmpPath); % phd : project home directory
                                  
% add necessary routines
addpath(genpath(fullfile(phd, 'src')));
% add path containg the main/index figure files
addpath(genpath(fullfile(phd, 'visualizations')));


%% creat project directory structure (pds)
global pds
pn = 'gpla_submission';

pds.pn               = pn;
pds.prj              = phd;

pds.src              = fullfile(pds.prj, 'src'); 
pds.ldat             = fullfile(pds.prj, 'localdata'); % local data
% pds.sim             = fullfile(pds.src, 'simulations'); 

% pds.datsrc.exp          = fullfile(pds.src, 'src', 'data', 'explorations'); 
% pds.datsrc.hnd          = fullfile(pds.src, 'src', 'data', 'handies'); 
% pds.datsrc.prc          = fullfile(pds.src, 'src', 'data', 'processed'); 
% pds.datsrc.sim          = fullfile(pds.src, 'src', 'data', 'simulations'); 

% pds.exp             = fullfile(pds.src, 'explorations');

% pds.viz             = fullfile(pds.src, 'visualizations');
% pds.vizabs             = fullfile(pds.viz, 'AbstractsPosters');
% pds.vizpap             = fullfile(pds.viz, 'papers');
% pds.vizpres            = fullfile(pds.viz, 'presentations');
% pds.vizconf            = fullfile(pds.viz, 'confPapers');

% pds.fig             = fullfile(pds.prj, 'figures');

% pds.sndb            = fullfile(pds.prj, 'sandbox');

clear phd tmpPath
