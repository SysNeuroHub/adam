function [n, out] = blindspot(nSims,plotIt)

%A recurrent network model of V1 that produces an attractive shift of
%position codes into the blindspot. Achieved through asymmetric
%intra-cortical weights at the margin of the blindspot, as might be
%expected if neurons "speak" to neurons on the other side for the purpose
%of filling in. 

nIter = 5;
N = 100;
sigmag = 0.4;

%Create a single layer of neurons that is recurrently connected with itself
n = deneveNet('myNet');

%Network layer.
addLayer(n,deneveLayer('retina',1,N));      %Hidden layer

%World layer
addLayer(n,deneveLayer('retinalWorld',1,N));      %Layer to represent target location in the world

%Create weights that include a blindspot, implemented by making neurons either side effectively closer together in cortex
retPrms = deneveLayer.defaultVMprms('NET2NET');
wldPrms = deneveLayer.defaultVMprms('WORLD2NET');

[a,b] = meshgrid(1:N,1:N);
interNeurDist = b-a;
distToBS = hypot(a-N/2,b-N/2);
acuity = 1 - (0.95.*exp((cos(distToBS.*(2*pi/N))-1)./(sigmag^2)));
interNeurDist = interNeurDist .* acuity;
wfun = @(ind) wldPrms.k.*exp((cos((2*pi/N).*ind)-1)/wldPrms.sigma^2); 
w = wfun(interNeurDist);

%Make the population recurrent
n.retina.setInput({n.retina,n.retinalWorld},'weights',{w,'vonMises'},'prms',[retPrms,wldPrms]);

%Allocate memory to log network state
allocLog(n,nIter,nSims*N);

n.evtFun.beforeUpdate = @beforeUpdate;
n.noisy(false);
n.plotIt = plotIt;

%=========== Run the simulation ==============
%Preallocate matrices
truRet= zeros(N,nSims);
for i = 1:N
    disp(num2str(i));
    for s = 1:nSims
        
        %Reset network layers to zero for each simulation
        n.reset();
        
        %Initialise world layers with random delta functions
        truRet(i,s) = i;        
        r=zeros(1,N);
        r(truRet(i,s))= 1;
        n.retinalWorld.setResp(r);

        %Vectors logging real positions
        plotThisOne = @(n) myPlot(n,truRet(i,s));
        n.evtFun.plot = plotThisOne;
        
        %Run the simulation
        n.run('nIter',nIter);
        
        %Calculate a point estimate in the final state
        estRet(i,s) = pointEstimate(n.retina);  
    end
end

%========== Statistical Data ===========%
%Calculate wrapped error
errRet = err(n.retina,estRet,truRet);

merr = x2rad(n.retina,errRet);
radmean = circ_mean(merr,[],2);
meanVal = rad2x(n.retina,radmean);

%Plot results
x = linspace(-N/2,N/2,N);
plot(x,meanVal,'o-','linewidth',3); hold on;
acuity(N/2,:) = acuity(N/2,:)./(max(acuity(N/2,:))-min(acuity(N/2,:)))-min(acuity(N/2,:));
plot(x,acuity(N/2,:)*diff(ylim)+min(ylim),'linewidth',3);
xlabel('Position rel. to blindspot');
ylabel('Localisation error');

figure
nPlots = 3;
subplot(nPlots,1,1);
toPlot = round(N*[-0.25 -0.075 0 0.075 0.25]);
retWeights = n.retina.input(1).weights;
for i=1:numel(toPlot)
    plot(x,retWeights(:,N/2-toPlot(i)),'linewidth',3); hold on;
end
plot(x, acuity(N/2,:),'linewidth',1);
xlabel('Position rel. to blindspot');
subplot(nPlots,1,2);
imagesc(x,x,retWeights); axis image; title('Connectivity weights');
subplot(nPlots,1,3); 
centerRegion = round(0.25*N):round(0.75*N);
imagesc(x(centerRegion),x(centerRegion),abs(interNeurDist(centerRegion,centerRegion))); axis image;
title('Cortical distance between units');
colorbar;

%Outputs
out.ret = {estRet,truRet,errRet};

function beforeUpdate(n)

%This function is called just before each iteration (t) of the network during a simulation
if n.t==1
    %Update each input layer with the input from the world
    setEnabled(n.retina,'retina',false);
    setEnabled(n.retina,'retinalWorld',true);
    
    %Switch off response normalisation
    n.retina.normalise = false;

elseif n.t==2
    %The stimuli and noise should only be present at t==1, so switch them off now.
    setEnabled(n.retina,'retina',true);
    setEnabled(n.retina,'retinalWorld',false);
    
    %Normalise all layers from herein
    n.retina.normalise = true;
end

function myPlot(n,realRetPos)

plotState(n.retina); ylim([0 40]);
plot([realRetPos, realRetPos],ylim,'k','linewidth',4);
if n.t==1
    pause(0.5)
else
    pause(0.1);
end
cla;
