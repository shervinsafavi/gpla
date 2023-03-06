function [spikeTrains_allTrLong, lfpPhases_allTrLong, n, selectedUnits, ...
          unwhitenOpr] = prep_SpkLfpData(spikeTrains_raw, lfpPhases_input, varargin)
% [spikeTrains_allTrLong, lfpPhases_allTrLong] == prep_SpkLfpData(spikeTrains_raw, lfpPhases_raw, varargin)

%% Handle optional inputs (varargin):
optionalVariables.flag_gPLVnrmlz      = [];     defaultValues{1} = 1;
optionalVariables.nSpikeThreshold     = [];     defaultValues{2} = NaN;
optionalVariables.unitSubset          = [];     defaultValues{3} = NaN;
optionalVariables.temporalWindow      = [];     defaultValues{4} = [];
optionalVariables.flag_origDimEigVec  = [];     defaultValues{5} = 0;   % not needed here
optionalVariables.statTestInfo        = [];     defaultValues{6} = [];  % not needed here
optionalVariables.iSV                 = [];     defaultValues{7} = 1;   % not needed here
optionalVariables.checkSameElecStuff_flag = [];     defaultValues{8} = 0;   % not needed here
optionalVariables.plvNrmlzMethed      = [];     defaultValues{9} = 'nSpk-square-root'; % not needed here
optionalVariables.flag_whitening      = [];     defaultValues{10} = 0;
optionalVariables.flag_lfpNrmlz       = [];     defaultValues{11} = 0;

optionalVariables = handleVarargin(varargin, optionalVariables, defaultValues);

%% Check the data format
% if length of spike and LFP vector are matching, otherwise inform user and
% exit the program

%%
n.Tr        = size(lfpPhases_input, 3);        % number of trials

%% whitening (if asked by the user)
% *** note: this is being done before epoch selection or any
% other preprocessing
if optionalVariables.flag_whitening == 1
    if isreal(lfpPhases_input)
        error('whitening is done on phase, u need to pass the analytical signal')
    else
        [tmpSig, witenOpr, unwhitenOpr] = ...
            whitenRed2(reshape(lfpPhases_input, size(lfpPhases_input, ...
                                                     1), []), NaN);
        % when is NaN, no rank reduction will be applied 
        % we needed to implemtn it like this to be compatible with
        % the rest of the codes
        lfpPhases_raw = reshape(tmpSig, size(tmpSig, 1), [], n.Tr);    
    end
elseif optionalVariables.flag_whitening == 2
    if isreal(lfpPhases_input)
        error('whitening is done on phase, u need to pass the analytical signal')
    else
        [lfpPhases_raw, witenOpr, unwhitenOpr] = whitenRed4(lfpPhases_input);
    end
else 
    lfpPhases_raw = lfpPhases_input;
    unwhitenOpr = [];

end

%% Data info
n.LfpCh     = size(lfpPhases_raw, 1);        % number of LFP channels
n.Sample    = size(lfpPhases_raw, 2);        % number of samples
n.SpkUnit   = size(spikeTrains_raw{1}, 1);   % number of spiking units


%% Specify the temporal window of analysis
if ~isempty(optionalVariables.temporalWindow) 
    % the specific window specified by used will be used for analysis
    % if there is start-end:
    if (...
        ~iscell(optionalVariables.temporalWindow) && ...
        numel(optionalVariables.temporalWindow) == 2) % assuming no one is interested in only 2 samples
        startInd = optionalVariables.temporalWindow(1);
        stopInd = optionalVariables.temporalWindow(2);
        selectedSamples = startInd : stopInd;
        %     elseif iscell(optionalVariables.temporalWindow)
        
        % if multiple time samples (similar across all trials or each trial sep within a cell array)
    else
        selectedSamples = optionalVariables.temporalWindow;
    end
else
    % *** need pnly 1 line !!
    startInd = 1;
    stopInd = n.Sample;
    selectedSamples = startInd : stopInd;
end

%% Specify the units to be included in the analysis

% pick units with spike count higher than 'nSpikeThreshold' (if asked by the user)
if ~isnan(optionalVariables.nSpikeThreshold)
    % find average spike count across trials 
    spikeCount = zeros(n.SpkUnit, n.Tr);
    for iTr = 1 : n.Tr
        spikeCount(:,iTr) = sum(spikeTrains_raw{iTr}, 2);
    end
    totSpikeCount = sum(spikeCount,2);
    selectedUnits_nSpkThrBased = find(totSpikeCount > optionalVariables.nSpikeThreshold);   
else
    selectedUnits_nSpkThrBased = (1:n.SpkUnit); % selecting all units
end

% pick units specified by the the user
if isnan(optionalVariables.unitSubset) 
    userSelectedUnits = (1:n.SpkUnit);
else
    userSelectedUnits = optionalVariables.unitSubset;
end

% selectedUnits = union of high spike guys and user selected units
selectedUnits = intersect(selectedUnits_nSpkThrBased, userSelectedUnits);

%% Prepare the data
% ~ allocate memory
spikeTrains{n.Tr} = [];
% F = exp(1i * lfpPhase_raw(:,startInd:stopInd, :));
% this 'if condition'  will help to save the some computation 
if ~isempty(optionalVariables.temporalWindow) 
    % preparing LFP phase matrix
    if iscell(selectedSamples)
        lfpPhases{n.Tr} = [];
        for iTr = 1 : n.Tr
            lfpPhases{iTr} = lfpPhases_raw(:,selectedSamples{iTr}, iTr);
        end
    else
        lfpPhases = lfpPhases_raw(:,selectedSamples, :);
    end
else
    % preparing LFP phase matrix
    lfpPhases = lfpPhases_raw;
end

% preparing spike matrix
if iscell(selectedSamples)
    for iTr = 1 : n.Tr
        spikeTrains{iTr} = spikeTrains_raw{iTr}(selectedUnits, selectedSamples{iTr});
    end    
else
    for iTr = 1 : n.Tr
        spikeTrains{iTr} = spikeTrains_raw{iTr}(selectedUnits, selectedSamples);
    end
end

%% Long spike trains and LFP 
spikeTrains_allTrLong = cell2mat(spikeTrains);


if iscell(selectedSamples)
    lfpPhases_allTrLong = cell2mat(lfpPhases);
else
    updated_nSample = numel(selectedSamples);
    lfpPhases_allTrLong = nan(n.LfpCh, n.Tr * updated_nSample);
    % tmpLfpData{n.Tr} = nan(n.LfpCh,
    for iTr = 1 : n.Tr
        tmpRange = (iTr-1)*updated_nSample +1 : (iTr-1)*updated_nSample + updated_nSample;
        lfpPhases_allTrLong(:,tmpRange) = lfpPhases(:,:, iTr);
    end
end

%% normalizng the signal if asked by user
if optionalVariables.flag_lfpNrmlz
    if isreal(lfpPhases_allTrLong)
        error('phase should not be normalized')
    else
        lfpPhases_allTrLong = bsxfun(@rdivide, lfpPhases_allTrLong, ...
                                     std(lfpPhases_allTrLong, [], 2));
    end
end