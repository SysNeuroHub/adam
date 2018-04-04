function deneveDemo(varargin)
% Build a Deneve network consisting of three input layers (retina, eye, and head) and one basis-function layer
% The model is the "feedback" model of Deneve, Latham, & Pouget (2001),
% "Efficient computation and cue integration with noisy population codes".
%
% See deneveNet.m and deneveLayer.m for more info.
p = inputParser;
p.addParameter('headWorldOn',false);    %Set to true to simulate an auditory stimulus (in addition to the visual stimulus)
p.addParameter('plotIt',true);          %Plot each iteration.
p.addParameter('nSims',100);
p.addParameter('nIter',10);
p.addParameter('addNoise',true);
p.addParameter('N',20);                 %Number of units per dimension in each layer
p.parse(varargin{:});
p = p.Results;

%% Build the network
    %For convenience, here we call a function that returns a ready-made model like that used 
    %in the original paper. Go to that function to see how networks are designed and interconnected.
n = deneveLathamPougetModel('N',p.N,'addNoise',p.addNoise);
addprop(n,'p');
n.p = p;
n.plotIt = p.plotIt;

if ~p.headWorldOn
    %Switch off the input to the head layer (i.e. "function approximation", as per the Deneve paper)
    setEnabled(n.head,n.headWorld.name,false);
else
    %Use the head-layer (e.g. "auditory") stimulus, i.e., "Cue integration"
end

%% Run simulations using different target and eye positions.

%Allocate memory to log network state
allocLog(n,n.p.nIter,n.p.nSims*n.p.N*n.p.N); %Allocates for all layers. Call allocLog on deneveLayer objects directly if you don't want all layers to log

%Create the function the calcluates the head position from retina and eye
headPos = @(r,e) mod((r + e - n.p.N/2)-1,n.p.N)+1;

%Run a simulation for every combination of target position and eye position
for i = 1:n.p.N
    disp(num2str(i));
    for j = 1:n.p.N
        
        %Get current scenario
        wld.retPos(i,j) = i;
        wld.eyePos(i,j) = j;
        wld.headPos(i,j) = headPos(i,j);
        
        %Use this info in a local plotting function at each iteration of the network
        if n.p.plotIt
            n.evtFun.plot = @(net) myPlot(net,wld.headPos(i,j));
        end
        
        %Create the input stimuli to each of the input layers (Delta functions)
        [retStim,eyeStim,headStim] = deal(zeros(1,n.p.N));
        retStim(wld.retPos(i,j))= 1;
        eyeStim(wld.eyePos(i,j))= 1;
        if n.p.headWorldOn, headStim(wld.headPos(i,j)) = 1; end
        
        %Set callback to reset network and again assign these stimuli to the network's world layers.
        n.evtFun.beforeSim = @(n) setStim(n,retStim,eyeStim,headStim);
        
        %Run the simulation
        for sInd = 1:n.p.nSims
            n.run('nIter',n.p.nIter);
            
            %Get a point estimate from each layer and return as error from true position
            [~,err.retinal(i,j,sInd)]    = pointEstimate(n.retinal, wld.retPos(i,j));
            [~,err.eye(i,j,sInd)]        = pointEstimate(n.eye, wld.eyePos(i,j));
            [~,err.head(i,j,sInd)]       = pointEstimate(n.head, wld.headPos(i,j));
        end
    end
end

%% Plot the results
%Calculate mean and standard deviation (in matrices of [Ret,Eye])
mret = x2rad(n.retinal,err.retinal);
meye = x2rad(n.eye,err.eye);
mhed = x2rad(n.head,err.head);

radmean.retinal = circ_mean(mret,[],3);
radmean.eye = circ_mean(meye,[],3);
radmean.head = circ_mean(mhed,[],3);

meanEst.ret = rad2x(n.retinal,radmean.retinal);
meanEst.eye = rad2x(n.eye,radmean.eye);
meanEst.hed = rad2x(n.head,radmean.head);

radstd.ret = circ_std(mret,[],[],3);
radstd.eye = circ_std(meye,[],[],3);
radstd.hed = circ_std(mhed,[],[],3);

stdEst.ret = rad2x(n.retinal,radstd.ret);
stdEst.eye = rad2x(n.eye,radstd.eye);
stdEst.hed = rad2x(n.head,radstd.hed);

%Plot heatmap of means and standard deviations
minVal = min(-1,min([meanEst.ret(:);meanEst.eye(:);meanEst.hed(:)]));
maxVal = max(1,max([meanEst.ret(:);meanEst.eye(:);meanEst.hed(:)]));
clims = [minVal,maxVal];
figure;
subplot(2,3,1);
imagesc(meanEst.ret,clims);
title('Ret Mean Error');
subplot(2,3,2);
imagesc(meanEst.eye,clims);
title('Eye Mean Error');
subplot(2,3,3);
imagesc(meanEst.hed,clims);
title('Head Mean Error');

minVal = min([stdEst.ret(:);stdEst.eye(:);stdEst.hed(:)]);
maxVal = max(1,max([stdEst.ret(:);stdEst.eye(:);stdEst.hed(:)]));
clims = [minVal,maxVal];

if p.nSims > 1
    subplot(2,3,4);
    imagesc(stdEst.ret,clims);
    title('Ret Std Dev');
    xlabel('eye position');
    ylabel('retinal position');
    subplot(2,3,5);
    imagesc(stdEst.eye,clims);
    title('Eye Std Dev');
    subplot(2,3,6);
    imagesc(stdEst.hed,clims);
    title('Head Std Dev');
end

%Plot the tuning curve of the central neuron in the hidden layer
if p.nSims == 1
    figure
    subplot(2,1,1);
    resp = reshape(n.basis.log,p.N,p.N);
    resp=cell2mat(cellfun(@(x) squeeze(x(end,p.N/2,p.N/2)),resp,'uniformoutput',false));
    surf(resp);
    title('Tuning curve of a hidden unit layer');
    subplot(2,1,2);
    toPlot = [round(0.35*p.N), round(0.4*p.N) 0.5*p.N];
    plot(resp(:,toPlot),'linewidth',4);
end

figure;
subplot(3,1,1);
h(1) = histogram(mret(:),10); title('retina');
subplot(3,1,2);
h(2) = histogram(meye(:),10); title('eye');
subplot(3,1,3);
h(3) = histogram(mhed(:),10); title('head');
minWidth = max(arrayfun(@(x) x.BinWidth,h));
set(h,'BinWidth',minWidth);
xlabel('Decoding error');
amUnifyAxesAcrossPlots([],'xy');
keyboard;
end

function setStim(n,retStim,eyeStim,headStim)
%(Re)set the network and assign the current stimulus
n.reset();
n.retWorld.setResp(retStim);
n.eyeWorld.setResp(eyeStim);
n.headWorld.setResp(headStim);
end

function myPlot(n,wldHeadPos)

subplot(2,1,1); cla
plotState(n.retinal); hold on
plotState(n.eye);
plotState(n.head);
plot([wldHeadPos,wldHeadPos],ylim,':k','lineWidth',3);

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
