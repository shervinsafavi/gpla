function [lfpVec, spkVec, gPLV, pValue, varargout] = tngpla(spikeTrains_raw, lfpPhases_raw, varargin)

%% Statistics on gPLV
% [lfpVec, spkVec, gPLV, pValue, gPLV_stats, couplingMatrix, singularValues, SVs_stats] = tngpla(spikeTrains, lfpPhase, ...
%     nSpikeThreshold, unitSubset, temporalWindow)
%
%
% EXAMPLE:      
% an example syntax (if it's needed)
%
% ------
% Input:
% 1     spikeTrains: 
%       expnation
% 2     lfpPhase: 
%       expnation
% (3)   nSpikeThreshold: scalar
%       (default value = 1)
%       number of cells (i.e. generate Poisson spike rrain for a population
% (4)   unitSubset: scalar
%       (default value = 1)
%       number of cells (i.e. generate Poisson spike rrain for a population
% (5)   temporalWindow: scalar
%       (default value = 1)
%       number of cells (i.e. generate Poisson spike rrain for a population
% (6)   flag_origDimEigVec: scalar
%       (default value = 1)
%       number of cells (i.e. generate Poisson spike rrain for a population
% (7)   statTestInfo: structure
%       (default value = 1)
%       number of cells (i.e. generate Poisson spike rrain for a population
% 
% Output:
% 1     lfpVec: 
% 2     spkVec: 
% 3     gPLV: 
% (4)   complex_gPLV
% (5)   gPLV_stats:
%          CI: 
%          nullDistribution: 
%          pValue: 
% (6)   PLV 
% (7)   complex_PLV (couplingMatrix)
% (8)   PLV_stats: 
%          CI: 
%          nullDistribution: 
%          pValue: 
% (9)   SVs_stats
% ------
% potential improvments:
% (1) ???
% ------
% Code Info:
%   creation: 2019-06-04 by SS and MB (shervin.safavi@tuebingen.mpg.de)
%   modification:
%       $ 2018-??-?? ???
% ------
% see also pla, ppla


%%    
optionalVariables.flag_gPLVnrmlz      = [];     defaultValues{1} = 0;
optionalVariables.nSpikeThreshold     = [];     defaultValues{2} = NaN;
optionalVariables.unitSubset          = [];     defaultValues{3} = NaN;
optionalVariables.temporalWindow      = [];     defaultValues{4} = [];
optionalVariables.flag_origDimEigVec  = [];     defaultValues{5} = 0;
optionalVariables.statTestInfo        = [];     defaultValues{6} = [];
optionalVariables.iSV                 = [];     defaultValues{7} = 1;
optionalVariables.sameElecCheckInfo_r = [];     defaultValues{8} = [];
optionalVariables.plvNrmlzMethed      = [];     defaultValues{9} = 'nSpk-square-root';
optionalVariables.flag_whitening      = [];     defaultValues{10} = 0;
optionalVariables.flag_lfpNrmlz       = [];     defaultValues{11} = 0;

optionalVariables = handleVarargin(varargin, optionalVariables, defaultValues);

%%
[spikeTrains_allTrLong, lfpPhases_allTrLong, n, selectedUnits, unwhitenOpr] = ...
    prep_SpkLfpData(spikeTrains_raw, lfpPhases_raw, varargin{:});

% adjut the tabel of channels according to, selected units inside gpla
if ~isempty(optionalVariables.sameElecCheckInfo_r)
    sameElecCheckInfo = optionalVariables.sameElecCheckInfo_r;
    sameElecCheckInfo.spkU_lfpCh_cnvrtTabel = ...
        optionalVariables.sameElecCheckInfo_r.spkU_lfpCh_cnvrtTabel(selectedUnits, :);
else
    sameElecCheckInfo = optionalVariables.sameElecCheckInfo_r;
end

statTestInfo = optionalVariables.statTestInfo;
if ~isempty(optionalVariables.statTestInfo)
    [gPLV, pValue,  lfpVec, spkVec_raw,  couplingMatrix, singularValues, gPLV_nullHypoReject, gPLV_stats, PLV_stats, SV_stats] = ...
        tnstataliz_gPLV(spikeTrains_allTrLong, lfpPhases_allTrLong, ...
                        statTestInfo, optionalVariables.iSV, ...
                        sameElecCheckInfo, ...
                        optionalVariables.plvNrmlzMethed, unwhitenOpr, optionalVariables.flag_gPLVnrmlz);
    
    gPLV_stats.pValue           = pValue;
    gPLV_stats.nullHypoReject   = gPLV_nullHypoReject;

    % collect all the statistics within a structure
    stats.gPLV_stats = gPLV_stats;
    stats.PLV_stats = PLV_stats;
    stats.SV_stats = SV_stats;

    % varargout{2} = SV_stats;  
else
    [gPLV, pValue,  lfpVec, spkVec_raw,  couplingMatrix, singularValues] = ...
        tnstataliz_gPLV(spikeTrains_allTrLong, lfpPhases_allTrLong, ...
                        statTestInfo, optionalVariables.iSV, ...
                        sameElecCheckInfo, optionalVariables.plvNrmlzMethed, ...
                        unwhitenOpr, optionalVariables.flag_gPLVnrmlz);
    pValue = NaN;
    stats = NaN;   % (o1)   all statistics

end

%%
rawSvdStuff.singularValues = singularValues;
rawSvdStuff.couplingMatrix = couplingMatrix;

varargout{1} = stats;   % (o1)   all statistics
varargout{2} = rawSvdStuff;
varargout{3} = selectedUnits;       

%%
% *** should go into a function
% spkVec and PLV matrix 
if logical(optionalVariables.flag_origDimEigVec)
    % everything should be returned in the original dimention
    % spkVec      = nan(n.SpkUnit, 1);
    spkVec      = nan(n.SpkUnit, size(spkVec_raw, 2));
    PLV         = nan(n.LfpCh, n.SpkUnit); 
    % complex_PLV = nan(n.LfpCh, n.SpkUnit);
    
    % spkVec(selectedUnits, 1)        = spkVec_raw; 
    spkVec(selectedUnits, :)        = spkVec_raw; 
    % PLV(:, selectedUnits)           = abs(M);
    % complex_PLV(:, selectedUnits)   = M;

else % unchanged
    spkVec                          = spkVec_raw;
    % PLV                             = abs(M);
    % complex_PLV                     = M;
end

