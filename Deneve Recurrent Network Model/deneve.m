function [n, out] = deneve(varargin)

p = inputParser;
p.addParameter('nSims',100);
p.addParameter('nIter',10);
p.addParameter('headWorldOn',true);
p.addParameter('addNoise',true);
p.addParameter('plotIt',false);
p.addParameter('N',20); %Number of units per dimension
p.addParameter('suppressLayer',[],@(x) isempty(x) | any(strcmpi({'eye','retinal','head'},x)));
p.parse(varargin{:});
p = p.Results;

Kg = 0.8;
sigmag = 0.4;

%Create layers ("population codes"). All layers are instances of the deneveLayer class.
n = deneveNet('myNet');

    %Network layers
addLayer(n,deneveLayer('hidden',p.N,p.N));      %Hidden layer
addLayer(n,deneveLayer('retinal',1,p.N));       %Visual layer of object on retina position
addLayer(n,deneveLayer('eye',1,p.N));           %Eye position signal
addLayer(n,deneveLayer('head',1,p.N));          %Object location relative to the head

    %World layers
addLayer(n,deneveLayer('retWorld',1,p.N));      %Layer to represent world location in retinal coordinates
addLayer(n,deneveLayer('eyeWorld',1,p.N));      %Layer to represent world location in eye centred coordinates
addLayer(n,deneveLayer('headWorld',1,p.N));     %Layer to represent world location in head centred coordinates

%Allocate memory to log network state
allocLog(n,p.nIter,p.nSims*p.N*p.N); %Allocates for all layers. Call allocLog on deneveLayer objects directly if you don't want all layers to log

%Assign a plot color to each layer
[n.retinal.plotSetts.lineColor,n.eye.plotSetts.lineColor,n.head.plotSetts.lineColor] = deal([0.8,0,0],[0,0,0.8],[0,0.8,0]);

%Specify the network architecture (who talks to whom?)
   %Parameters for von Mises connectivity profiles
retWtPrms = deneveLayer.defaultVMprms('NET2NET');
eyeWtPrms = deneveLayer.defaultVMprms('NET2NET');
headWtPrms = deneveLayer.defaultVMprms('NET2NET');
wldWtPrms = deneveLayer.defaultVMprms('WORLD2NET');
retWtPrms.dim = 1;          %Feed the retinal input to dim 1 of the hidden layer
eyeWtPrms.dim = 2;          %Feed the eye input to dim 2 of the hidden layer
headWtPrms.dim = -1;        %Feed the head input along the diagonal

    %Connect the layers
n.hidden.setInput({n.retinal,n.eye,n.head},'prms',[retWtPrms,eyeWtPrms,headWtPrms]);
n.retinal.setInput({n.hidden,n.retWorld},'prms',[retWtPrms,wldWtPrms]);
n.eye.setInput({n.hidden,n.eyeWorld},'prms',[eyeWtPrms,wldWtPrms]);
if p.headWorldOn
    n.head.setInput({n.hidden,n.headWorld},'prms',[headWtPrms,wldWtPrms]);
else
    n.head.setInput(n.hidden,'prms',headWtPrms);
end

%=========== Run the simulations ==============
[truRet,truEye,truHed] = deal(zeros(p.N,p.N,p.nSims));
[estRet,estEye,estHed] = deal(zeros(p.N,p.N,p.nSims,p.nIter));
for i = 1:p.N
    disp(num2str(i));
    for j = 1:p.N
        for s = 1:p.nSims
            %Reset network layers to zero for each simulation
            n.reset;
            
            %Get current scenario
            realRetPos = i;
            realEyePos = j;
            realHedPos = mod((realRetPos + realEyePos - p.N/2)-1,p.N)+1;
            
            %Set the retinal stimulus
            r = zeros(1,p.N);
            r(realRetPos)= 1;
            n.retWorld.setResp(r);
            
            %Set the current eye position
            r = zeros(1,p.N);
            r(realEyePos)= 1;
            n.eyeWorld.setResp(r);
            
            %Switch between delta function and matrix of zeros for head input
            r = zeros(1,p.N);   
            if p.headWorldOn
                r(realHedPos) = 1;
            end
            n.headWorld.setResp(r);
            
            %Run the simulation
            for t=1:p.nIter
                
                %If we are just starting, read in only the world. 
                isTimeZero=t==1;
                setEnabledLayers(n,isTimeZero,p.headWorldOn);
                
                %Calculate the new state of the network
                update(n,isTimeZero);
 
                if isTimeZero
                    %Apply response suppression if requested
                    if ~isempty(p.suppressLayer)                      
                        gainfun(n.(p.suppressLayer),n.(p.suppressLayer).nUnits/2,Kg,sigmag);
                    end
                    
                    %Add noise to the input layers (not the hidden layer)
                    if p.addNoise
                        addNoise(n.retinal);
                        addNoise(n.eye);
                        addNoise(n.head);
                    end
                end
                
                %============== Plot the simulation==============%
                %Plot vector responses of ret, eye, and hed
                if p.plotIt
                    subplot(2,1,1); cla
                    plotState(n.retinal); hold on
                    plotState(n.eye);
                    plotState(n.head);
                    plot([realHedPos,realHedPos],ylim,':k','lineWidth',3);
                    
                    %Plot 2D matrix response of hid
                    subplot(2,1,2); cla
                    plotState(n.hidden);
                    title(num2str(t));
                    
                    %Pauses for plotting
                    if t==1
                        pause(1);
                    else
                        pause(0.15);
                    end
                end

                %Log true positions
                truRet(i,j,s) = realRetPos;
                truEye(i,j,s) = realEyePos;
                truHed(i,j,s) = realHedPos;
                
                %Log point estimates for each layer
                estRet(i,j,s,t) = pointEstimate(n.retinal);
                estEye(i,j,s,t) = pointEstimate(n.eye);
                estHed(i,j,s,t) = pointEstimate(n.head);
                
                %Log the current state of each layer
                logState(n.retinal,'newLog',isTimeZero);
                logState(n.eye,'newLog',isTimeZero);
                logState(n.head,'newLog',isTimeZero);
                logState(n.hidden,'newLog',isTimeZero);
            end
        end
    end
end

%========== Statistical Data ===========%
%Calculate wrapped error
errRet = err(n.retinal,estRet(:,:,:,p.nIter),truRet);
errEye = err(n.eye,estEye(:,:,:,p.nIter),truEye);
errHed = err(n.head,estHed(:,:,:,p.nIter),truHed);

%Outputs
out.ret = {estRet,truRet,errRet};
out.eye = {estEye,truEye,errEye};
out.hed = {estHed,truHed,errHed};

%Plot the tuning curve of the central neuron in the hidden layer
if p.nSims == 1
    subplot(2,1,1);
    resp = reshape(n.hidden.log,p.N,p.N);
    resp=cell2mat(cellfun(@(x) squeeze(x(end,p.N/2,p.N/2)),resp,'uniformoutput',false));
    surf(resp);
    title('Tuning curve of a hidden unit layer');
    subplot(2,1,2);
    toPlot = [round(0.35*p.N), round(0.4*p.N) 0.5*p.N];
    plot(resp(:,toPlot),'linewidth',4);
end
end


function setEnabledLayers(n,isTimeZero,headOn)

%Switch between world and hidden layer inputs for ret, eye and hed
setEnabled(n.retinal,n.hidden.name,~isTimeZero);
setEnabled(n.retinal,n.retWorld.name,isTimeZero);
setEnabled(n.eye,n.hidden.name,~isTimeZero);
setEnabled(n.eye,n.eyeWorld.name,isTimeZero);

if headOn
    setEnabled(n.head,n.hidden.name,~isTimeZero);
    setEnabled(n.head,n.headWorld.name,isTimeZero);
end
end

function update(n,isTimeZero)
normalise = ~isTimeZero;
n.hidden.update(normalise);
n.retinal.update(normalise);
n.eye.update(normalise);
n.head.update(normalise);
end