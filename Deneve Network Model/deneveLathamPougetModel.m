function n = deneveLathamPougetModel(varargin)

% Returns a deneveNet object that is set up to match the FEEDBACK model of
% Deneve, Latham, & Pouget (2001), "Efficient computation and cuintegration
% with noisy population codes." See deneveDemo.m.

%(one exception: the basis layer here has the same dimensionality as the input
%layers; the original model had N/2 units along each dimension)
p=inputParser;
p.addParameter('N',20);                 %Number of units per dimension in each layer
p.addParameter('addNoise',true);        %Use noisy population responses
p.parse(varargin{:});
p = p.Results;

%% Build the network

%Create layers ("population codes"). All layers are instances of the deneveLayer class.
n = deneveNet('myNet');

    %Neural layers
addLayer(n,deneveLayer('basis',p.N,p.N));       %Basis function layer
addLayer(n,deneveLayer('retinal',1,p.N));       %Visual layer coding object location on the retina
addLayer(n,deneveLayer('eye',1,p.N));           %Eye position layer coding current direction of gaze
addLayer(n,deneveLayer('head',1,p.N));          %Object location relative to the head (e.g. auditory)

    %World layers
addLayer(n,deneveLayer('retWorld',1,p.N));      %The true state of the world (retinal coordinates)
addLayer(n,deneveLayer('headWorld',1,p.N));     %The true state of the world (head coordinates)
addLayer(n,deneveLayer('eyeWorld',1,p.N));      %The true state of the eye (head coordinates)

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
n.head.setInput({n.basis,n.headWorld},'prms',[headWtPrms,wldWtPrms]);

%Assign a plot color to each neural layer
n.retinal.plotSetts.lineColor   = [0.8,0,0];
n.eye.plotSetts.lineColor       = [0,0,0.8];
n.head.plotSetts.lineColor      = [0,0.8,0];

%The order in which layers are updated matters. For example, a feedforward
%network has to be updated in the direction of information flow.
%Here, we need the input layers to take in the state of the world first
%(i.e. pool from the "world" layers), then we need to update the basis layer,
%then update the input layers, then basis etc...
%So, implement these switches/ordering in a custom beforeUpdate() function (specified at the bottom of this script).
%This custom script also ensures that noise is added only at t==0 and not thereafter.
n.evtFun.beforeUpdate = @beforeUpdate;

function beforeUpdate(n)

%This function is called just before each iteration (t) of the network during a simulation

if n.t==1
    %Update only the input layers
    n.basis.locked = true;
    
    %Switch off response normalisation
    n.retinal.normalise = false;
    n.eye.normalise = false;
    n.head.normalise = false;
    
    %Add noise to the responses
    if n.p.addNoise
        n.noisy(true);
    end 
elseif n.t==2
    %Update the basis layer only
    n.basis.locked = false;
    n.retinal.locked = true;
    n.eye.locked = true;
    n.head.locked = true;
    n.noisy(false);
    
elseif n.t==3;
    %Update and normalise all layers therafter
    n.basis.locked = false;
    n.retinal.locked = false;
    n.eye.locked = false;
    n.head.locked = false;
    n.retinal.normalise = true;
    n.eye.normalise = true;
    n.head.normalise = true;
end