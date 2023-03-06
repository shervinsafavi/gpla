function [stimVec, classIdx] = generWavStim(transientParams)
%GENERWAVSTIM - generate a time sequence of stimulations to generate
%various wave equation propagating patterns.
%GENERWAVSTIM(stimSequence,interStim,stimPer,stimRepet)
%
%  inputs
%
%	stimSequence - integer vector of stimulus sequence
%
%	interStim - number of time points between each stimulation periods
%
%	stimPer - period of repetitive stimulations
%
%	stimRepet - number of periodic repetitions of a single stimulation
%
%
%  outputs
%
%	stimVec - time sequence of stimulation
%
%
% Author : Michel Besserve, MPI for Intelligent Systems, MPI for Biological Cybernetics, Tuebingen, GERMANY
structunpack(transientParams)
pauseInter=zeros(interStim,1);
stimInter=zeros(stimPer,1);
stimVec=pauseInter;
classIdx=0*pauseInter;
for kstim=1:length(stimSequence)
    stimVec=cat(1,stimVec,repmat([stimSequence(kstim);stimInter],stimRepet,1),pauseInter);
    classIdx=cat(1,classIdx,repmat(stimSequence(kstim)*[1;1+0*stimInter],stimRepet,1),0*pauseInter);
end
