function [lfpVec, spkVec, gPLV, varargout] = gpla_core(spikeTrains_allTrLong, lfpPhases_allTrLong, varargin)
% [lfpVec, spkVec, gPLV, complex_gPLV, 
%               complex_PLV, singularValues, SVspectrum_recBsd, M_rankReduced_recBsd] = gpla(spikeTrains_allTrLong, lfpPhases_allTrLong, varargin)
% it takes directly trial concatenated data for both LFP and spikes
%
% EXAMPLE:      
% an example syntax (if it's needed)
%
% ------
% Input:
% 1     spikeTrains: trial concatenated spike trains: unit x time bin
% 2     lfpPhase: trial concatenated lfp : channel x time sample (same as
% time bin)
% (3) binary flag for normalization of gpla: if true normalize by square
% root of number of units
% (4) iSV: integer index of the sv to return, by default 1: returns the first
% (largest) singular vector/value
%
% Output:
% 1     lfpVec: 
% 2     spkVec: 
% 3     gPLV: 
% (4)   complex_gPLV
% (5)   complex_PLV (couplingMatrix)
% (6)   Singular values spectrum
% (7)  SVspectrum_recBsd ? spectrum when removed the top sv
% (8)  M_rankReduced_recBsd ? reduced rank matrix when removed the top sv
%
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

%% Handle optional inputs (varargin):
optionalVariables.flag_gPLVnrmlz      = [];     defaultValues{1} = 0;
optionalVariables.iSV                 = [];     defaultValues{2} = 1;
optionalVariables.sameElecCheckInfo   = [];     defaultValues{3} = [];
optionalVariables.plvNrmlzMethed      = [];     defaultValues{4} = 'nSpk-square-root';
optionalVariables.unwhitenOpr         = [];     defaultValues{5} = [];

optionalVariables = handleVarargin(varargin, optionalVariables, defaultValues);

%% prepare the data
% pick the relevant portion of the signal and units 
% [spikeTrains, lfpPhases, n] = prep_SpkLfpData(spikeTrains_raw, lfpPhases_raw, varargin{:});

% n.LfpCh     = size(lfpPhases_raw, 1);        % number of LFP channels
% n.Sample    = size(lfpPhases_raw, 2);        % number of samples
% n.SpkUnit   = size(spikeTrains_raw, 1);   % number of spiking units
% 
% [selectedUnits, spikeTrains_allTrLong, lfpPhases_allTrLong] = ...
%     prep_SpkLfpData(spikeTrains_raw, lfpPhases_raw, varargin{:});

n.LfpCh     = size(lfpPhases_allTrLong, 1);        % number of LFP channels
n.Sample    = size(lfpPhases_allTrLong, 2);        % number of samples
n.SpkUnit   = size(spikeTrains_allTrLong, 1);   % number of spiking units

% if sum(full(spikeTrains_allTrLong(:))) == 0
%     error('There is no spike!') 
% end

%% Compute the [complex] coupling matrix
% M = cmpt_couplingMatrixMB(full(spikeTrains_allTrLong), ...
%                           lfpPhases_allTrLong);

[M, totSpikeCounts] = tncmpt_couplingMatrix(full(spikeTrains_allTrLong), lfpPhases_allTrLong, ...
                                            optionalVariables.plvNrmlzMethed, optionalVariables.sameElecCheckInfo);


[singularLfpVecs_raw, singularSpkVecs_raw, singularValues] = fctrz_couplingMatrix(M);
% [SVspectrum_recBsd, M_rankReduced_recBsd] = cmpt_recBsdSVspectrum(M);

%% unwhiten the LFP if the operator is passed by user
if ~isempty(optionalVariables.unwhitenOpr)
    % unwhiten the signal
    tmpSingularLfpVecs = optionalVariables.unwhitenOpr * singularLfpVecs_raw;
    % and normalize it
    singularLfpVecs = bsxfun(@rdivide, tmpSingularLfpVecs, sqrt(sum(abs(tmpSingularLfpVecs).^2,1)));
else 
    % leave it as it was
    singularLfpVecs = singularLfpVecs_raw;
end

%% further normalizttion of spike vecttors
% spkVec 
% if we normalized the entries of M with square root of spike count
% we need the following additional normalization
switch optionalVariables.plvNrmlzMethed
  case 'nSpk-square-root'
    tmpSingularSpkVecs = bsxfun(@rdivide, singularSpkVecs_raw, ...
                                totSpikeCounts .^ .5) ;
    % singularSpkVecs = bsxfun(@rdivide,tmpSingularSpkVecs, sqrt(sum(abs(tmpSingularSpkVecs).^2,1)));
    % since we might have zero spikes, we should do nansum
    singularSpkVecs = bsxfun(@rdivide,tmpSingularSpkVecs, sqrt(nansum(abs(tmpSingularSpkVecs).^2,1)));
   % singularSpkVecs = tmpSingularSpkVecs;
otherwise
    singularSpkVecs = singularSpkVecs_raw;
end


%% phase rotation 
% phase normalization of spike and LFP eigen vector

if ~ischar(optionalVariables.iSV)

    if numel(optionalVariables.iSV) == 1
        meanPhase_lfpVec = angle(nanmean(singularLfpVecs(:, optionalVariables.iSV)));
        spkVec                          = singularSpkVecs(:, optionalVariables.iSV) * exp(-1i * meanPhase_lfpVec);
        lfpVec = singularLfpVecs(:, optionalVariables.iSV) * exp(-1i * meanPhase_lfpVec);

    elseif numel(optionalVariables.iSV) > 1
        % *** there is more effcient ways of implementing it, just for now
        for iSelSV = optionalVariables.iSV
            meanPhase_lfpVec = angle(nanmean(singularLfpVecs(:, iSelSV)));
            spkVec(:, iSelSV) = singularSpkVecs(:, iSelSV) * exp(-1i * meanPhase_lfpVec);
            lfpVec(:, iSelSV) = singularLfpVecs(:, iSelSV) * exp(-1i * meanPhase_lfpVec);
        end
    end
else % if is a char and is all
     % *** there is more effcient ways of implementing it, just for now
    for iSelSV = 1 : numel(singularValues)
        meanPhase_lfpVec = angle(nanmean(singularLfpVecs(:, iSelSV)));
        spkVec(:, iSelSV) = singularSpkVecs(:, iSelSV) * exp(-1i * meanPhase_lfpVec);
        lfpVec(:, iSelSV) = singularLfpVecs(:, iSelSV) * exp(-1i * meanPhase_lfpVec);
    end
    % should be completed later, all is passed all we passed     
end


%% PLV
PLV                             = abs(M);
complex_PLV                     = M;


% gPLV
% gPLV = singularValues(1); %

% *** I'm not sure if this normaliztion is any more applicable due
% to the different normaliztion of entries of coupling matrix
% in the case where whe normlaize by squre nSpk and without
% whitening is safe (-> no amplitude)

if optionalVariables.flag_gPLVnrmlz
        % gPLV = singularValues(1) * (n.LfpCh * n.SpkUnit)^-.5;
    % gPLV = singularValues(optionalVariables.iSV) * (n.LfpCh * n.SpkUnit)^-.5;
        gPLV = singularValues(1) * (size(M, 1) * size(M, 2))^-.5
else
    % gPLV = singularValues(optionalVariables.iSV);
    gPLV = singularValues(1);
end

% complex_gPLV
complex_gPLV = gPLV * exp(2i * meanPhase_lfpVec); 

%% Outputs 
% lfpVec is already assigned        % 1     LFP vector
% spkVec is already assigned        % 2     spike vector
% gPLV   is already assigned        % 3     gPLV
varargout{1} = complex_gPLV;        % (4)   complex gPLV
                                    %varargout{2} = PLV;                 % (5)   PLV absolute value: useless:
                                    %remove
varargout{2} = complex_PLV;         % (6)   complex PLV or coupling matrix 
varargout{3} = singularValues;      % (7)   Singular values spectrume
                                    % varargout{4} = SVspectrum_recBsd;   % (8)
                                    % varargout{5} = M_rankReduced_recBsd;% (9)
                                    % 
end

% function [SVspectrum_recBsd, varargout] = cmpt_recBsdSVspectrum(M)
% % substraction the first component to compute other svs: not used anymore
% iSV =  1;
% [singularLfpVecs, ~, SVs] = fctrz_couplingMatrix(M);
% SVspectrum_recBsd(iSV,1) = SVs(1);

% % u1 = singularLfpVecs(1,:);
% u1 = singularLfpVecs(:,1);
% projectorMat = (u1 * u1') / (norm(u1) ^ 2);
% M_wLSC = M;
% % M_xLSC = M_wLSC - (SVs(1) * projectorMat) * M_wLSC; 
% M_xLSC = M_wLSC - projectorMat * M_wLSC; 
% M_rankReduced_recBsd(:,:, 1) = M_xLSC;

% % for iSV = 2 : min(size(M))
% for iSV = 2 : 10
% %     [singularLfpVecs, ~, SVs] = fctrz_couplingMatrix(M_xLSC);    
% 	[singularLfpVecs, ~, SVs] = fctrz_couplingMatrix(M_rankReduced_recBsd(:,:, iSV-1));    
%     SVspectrum_recBsd(iSV,1) = SVs(1);

%     u1 = singularLfpVecs(:,1);
%     M_wLSC = M_xLSC; 
%     projectorMat = (u1 * u1') / (norm(u1) ^ 2);
% %     M_xLSC = M_wLSC - (SVs(1) * projectorMat) * M_wLSC; 
%     M_xLSC = M_wLSC - projectorMat * M_wLSC; 
%     M_rankReduced_recBsd(:,:, iSV) = M_xLSC;
% end
% varargout{1} = M_rankReduced_recBsd;
% end

