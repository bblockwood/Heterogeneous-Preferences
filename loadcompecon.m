% LOADCOMPECON Loads the CompEcon toolbox  
% Placed in own function because the local path is different from the
% server path.

function loadcompecon()
cepath='~/Dropbox/Programming/Matlab/utilities/compecon2011/';
path([cepath 'CEtools;' cepath 'CEdemos'],path);
