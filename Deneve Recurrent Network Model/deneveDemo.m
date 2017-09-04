function [n, out] = deneveDemo(varargin)
% Build a Deneve network consisting of three input layers (retina, eye, and head) and one basis-function layer
% The model is the "feedback" model of Deneve, Latham, & Pouget (2001),
% "Efficient computation and cue integration with noisy population codes".

p = inputParser;
p.addParameter('headWorldOn',false);    %Set to true to simulate an auditory stimulus (in addition to the visual stimulus)
p.addParameter('plotIt',true);          %Plot each iteration.
p.addParameter('nSims',100);
p.addParameter('nIter',10);
p.addParameter('N',20);                 %Number of units per dimension in each layer
p.parse(varargin{:});

%% Build the network

%   FIX THIS -  should we add a prop for each thing here or add one prop with all settings?
%Create layers ("population codes"). All layers are instances of the deneveLayer class.
n = deneveNet('myNet');
addprop(n,'p');
n.p = p.Results;
n.plotIt = n.p.plotIt;

    %Neural layers
addLayer(n,deneveLayer('basis',n.p.N,n.p.N));     %Basis function layer
addLayer(n,deneveLayer('retinal',1,n.p.N));       %Visual layer coding object location on the retina
addLayer(n,deneveLayer('eye',1,n.p.N));           %Eye position layer coding current direction of gaze
addLayer(n,deneveLayer('head',1,n.p.N));          %Object location relative to the head (e.g. auditory)

    %World layers
addLayer(n,deneveLayer('retWorld',1,n.p.N));      %The true state of the world (retinal coordinates)
addLayer(n,deneveLayer('headWorld',1,n.p.N));     %The true state of the world (head coordinates)
addLayer(n,deneveLayer('eyeWorld',1,n.p.N));      %The true state of the eye (head coordinates)

%Specify the network architecture (who talks to whom?)
   %Parameters for von Mises connectivity profiles
retWtPrms       = deneveLayer.defaultVMprms('NET2NET');
eyeWtPrms       = retWtPrms;
headWtPrms      = retWtPrms;
wldWtPrms       = deneveLayer.defaultVMprms('WORLD2NET');
retWtPrms.dim   = 1;         %The basis layer will be tuned for target position on retina along dim 1
eyeWtPrms.dim   = 2;         %The basis layer will be tuned for eye position along dim 2
headWtPrms.dim  = -1;        %The basis layer will be tuned for position rel. head along the diagonal (i.e. h = r + e)

    %Connect the layers
n.basis.setInput({n.retinal,n.eye,n.head},'prms',[retWtPrms,eyeWtPrms,headWtPrms]);
n.retinal.setInput({n.basis,n.retWorld},'prms',[retWtPrms,wldWtPrms]);
n.eye.setInput({n.basis,n.eyeWorld},'prms',[eyeWtPrms,wldWtPrms]);
if n.p.headWorldOn
    n.head.setInput({n.basis,n.headWorld},'prms',[headWtPrms,wldWtPrms]);
else
    n.head.setInput(n.basis,'prms',headWtPrms);
end

%% Run simulations using different target and eye positions.

%Allocate memory to log network state
allocLog(n,n.p.nIter,n.p.nSims*n.p.N*n.p.N); %Allocates for all layers. Call allocLog on deneveLayer objects directly if you don't want all layers to log

%Assign a plot color to each neural layer
n.retinal.plotSetts.lineColor   = [0.8,0,0];
n.eye.plotSetts.lineColor       = [0,0,0.8];
n.head.plotSetts.lineColor      = [0,0.8,0];

%We want the network to listen to the world-layers and add noise only at t==1 and not
%thereafter. So, implement this switch in a custom beforeUpdate() function
%(specified at the bottom of this script).
n.evtFun.beforeUpdate = @beforeUpdate;

%Run the sims.
for i = 1:n.p.N
    disp(num2str(i));
    for j = 1:n.p.N
        for sInd = 1:n.p.nSims
            
            %Get current scenario
            wld.retPos = i;
            wld.eyePos = j;
            wld.headPos = mod((wld.retPos + wld.eyePos - n.p.N/2)-1,n.p.N)+1;
            
            %Use this info in a local plotting function at each iteration of the network
            if n.p.plotIt
                n.evtFun.plot = @(net) myPlot(net,wld);
            end
            
            %Create the input stimuli to each of the input layers (Delta functions)
            [retStim,eyeStim,headStim] = deal(zeros(1,n.p.N));
            retStim(wld.retPos)= 1;
            eyeStim(wld.eyePos)= 1;  
            if n.p.headWorldOn, headStim(wld.headPos) = 1; end
            
            %Assign the stimuli to the network's world layers.
            n.retWorld.setResp(retStim);
            n.eyeWorld.setResp(eyeStim);
            n.headWorld.setResp(headStim);

            %Run the simulation
            n.run();
        end
    end
end
end

function myPlot(n,wld)

subplot(2,1,1); cla
plotState(n.retinal); hold on
plotState(n.eye);
plotState(n.head);
plot([wld.headPos,wld.headPos],ylim,':k','lineWidth',3);

%Plot 2D matrix response of hid
subplot(2,1,2); cla
plotState(n.basis);
title(num2str(n.t));

%Pauses for plotting
if n.t==1
    pause(1);
else
    pause(0.15);
end
end

function beforeUpdate(n)
isTimeZero = n.t==1;

%If time zero, switch off all inputs except those coming from world layers
%Also switch off normalisation/transfer.
%The opposite thereafter.

    %Retina
setEnabled(n.retinal,n.basis.name,~isTimeZero);
setEnabled(n.retinal,n.retWorld.name,isTimeZero);
n.retinal.normalise = ~isTimeZero;
    %Eye
setEnabled(n.eye,n.basis.name,~isTimeZero);
setEnabled(n.eye,n.eyeWorld.name,isTimeZero);
n.eye.normalise = ~isTimeZero;
    %Head
setEnabled(n.head,n.basis.name,~isTimeZero);
if n.p.headWorldOn
    setEnabled(n.head,n.headWorld.name,isTimeZero);
end
n.head.normalise = ~isTimeZero;
    %Basis
setEnabled(n.basis,n.retinal.name,~isTimeZero);
setEnabled(n.basis,n.eye.name,~isTimeZero);
setEnabled(n.basis,n.head.name,~isTimeZero);

%Add noise only if this is the first time point
n.noisy(isTimeZero);
end