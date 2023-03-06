function [recentLfpPh, varargout] = tpp_filt_recenter(lfpLikeSig, freqCenter, halfFilterWidth, ...
                           SF, filterOrder, nIteCent)
% this function do the filtering and phase recursive recentering of
% LFP data

signalParams.nTr = size(lfpLikeSig, 3);
signalParams.nCh = size(lfpLikeSig, 1);

%% Filter parameters
% spike-lfp analysis
% freqBand = globalDynamicsParams.oscFreq(iFreq) + [-7.5 +7.5];
freqBand = freqCenter + [-halfFilterWidth +halfFilterWidth];
wn = freqBand / (SF/2);
% filterOrder = 2;
[b, a] = butter(filterOrder, wn, 'bandpass');


%% filter LFP
for iTr = 1 : signalParams.nTr
    for iCh = 1 : signalParams.nCh
        tmpLfp = lfpLikeSig(iCh,:, iTr);
        tmpLfp_filtered = filtfilt(b, a, tmpLfp);
        tmpAnlcLfp = hilbert(tmpLfp_filtered);
        
        analLfp(iCh,:, iTr) = tmpAnlcLfp;
        
        for kite = 1 :  nIteCent
            tmpAnlcLfp = tmpAnlcLfp./ abs(tmpAnlcLfp);
            tmpAnlcLfp = tmpAnlcLfp - mean(tmpAnlcLfp);
        end 
        
        recentLfpPh(iCh,:, iTr) = angle(tmpAnlcLfp);
        
    end
end


varargout{1} = analLfp;