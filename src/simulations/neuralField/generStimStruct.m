function StimStruct = generStimStruct(szy,dsz)
% generate edp stimulation parameters
% inputs:
% szy: y dimension of the grid in meters

    dp = 0.0016; % droplet probability per one time step
    if nargin<2
        dsz = szy/5; % spatial scale of event
    end
    da = .6; % droplet amplitude
    ba = .15;%bar amplitude
    stimtype = 'droplet';
    x0d = szy/2;
    y0d = szy/2;
    Osc = 1;
    StimStruct{1} = structpack({'stimtype','Osc','x0d','y0d','da','dsz','dp'});
    stimtype = 'bar';
    Osc = 1;
    StimStruct{2} = structpack({'stimtype','Osc','x0d','y0d','da','ba','dsz','dp'});

    stimtype = 'droplet';
    x0d = szy/2;
    y0d = szy/2; 
    Osc = 0;
    StimStruct{3} = structpack({'stimtype','Osc','x0d','y0d','da','ba','dsz','dp'});
    
    stimtype = 'noise';
    x0d = szy/2;
    y0d = szy/2;
    StimStruct{4} = structpack({'stimtype','Osc','x0d','y0d','da','ba','dsz','dp'});
    x0d = szy/2;
    y0d = szy/2;
    StimStruct{5} = structpack({'stimtype','Osc','x0d','y0d','da','ba','dsz','dp'});

    x0d = szy/2;
    y0d = szy/2;
    StimStruct{6} = structpack({'stimtype','Osc','x0d','y0d','da','ba','dsz','dp'});
    stimtype = 'bar';
    xside = 0;
    StimStruct{7} = structpack({'stimtype','Osc','x0d','y0d','da','ba','dsz','dp','xside'});
    
    
    stimtype = 'filt_noise';
    x0d = szy/2;
    y0d = szy/2;
    StimStruct{8} = structpack({'stimtype','Osc','x0d','y0d','da','ba','dsz','dp'});
end

