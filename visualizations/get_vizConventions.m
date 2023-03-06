function vc = get_vizConventions()
% get_vizConventions()
% this function give the parameters are needed only in visualization (visualization conventions)

%% size subplot labels
vc.sblfz = 15;

%% figure 2 of the manuscript
vc.f1.w = 17.4;
vc.f1.h = 21.75;


%% figure 4 of the manuscript
vc.f10.gfs = 7;
vc.f10.lfs = 8;
vc.f10.w = 17.4;
vc.f10.h = vc.f10.w * 874/1066;
vc.f10.sblfz = 14;

%% figure 6 of the manuscript
vc.f11.w = 17.4;
vc.f11.h = 23;


%% Visualization conventions

vc.yos = 1.2;

vc.lw_lfp = 1.55;
vc.lw_spk = 1.5;
vc.lfpv = [0 0 1];
vc.spkv = [1 0 0];
vc.illusCases = [1 1 .5];

vc.f7.barWidth = .7;

%% for figure 7 
vc.f7.barWidth = 0.7;
vc.f7.gtColor = [.5 .5 .5];
vc.f7.unitColor = [0 0 0];

vc.f7.bulletSize = 20;
vc.f7.offset = .065;
vc.f7.fontSize = 10;
vc.f7.alpha = 0.05;

%% 9 
vc.f9.w = 17.4;
vc.f9.h = 15;
vc.f9.gfs = 6; % global font size
vc.f9.bulletSize_tiny = 1;
vc.f9.bulletSize = 10;

% font size
vc.f9.couplingCmpr_fs = 7; % global font size

% tick size % ts = tick size
vc.f9.couplingCmpr_ts = .04; 
vc.f9.ts.oscDemo = 0.025;

vc.f9.lfs = 8;

%% 11
vc.f11.EIphaseShift.lw = 2;

%% 10
vc.f10.imts = 0.04;

%% figure 5
vc.f5.ppLW = 1.5;
end

