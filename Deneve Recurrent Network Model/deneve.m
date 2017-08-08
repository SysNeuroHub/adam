function [net, wld, out] = deneve(nSims,headon,ploton,gainon)

%Add another input that is plot or not
%Set constants for functions and loops
nIter = 20;
N = 20;
Kw = 1;
sigmaw = 0.37;
K = 20;
v = 1;
sigma = 0.40;
Kg = 0.8;
sigmag = 0.4;

%Create layers ("population codes"). All layers are instances of the deneveLayer class.
%Network layers
net.hid = deneveLayer('hidden',N,N);    %Hidden layer
net.ret = deneveLayer('retinal',1,N);   %Object location on the retina
net.eye = deneveLayer('eye',1,N);       %Eye position signal
net.hed = deneveLayer('head',1,N);      %Object location relative to the head

%World layers
wld.ret = deneveLayer('retinalWorld',1,N);  %Layer to represent world location in retinal coordinates
wld.eye = deneveLayer('eyeWorld',1,N);      %Layer to represent world location in eye centred coordinates
wld.hed = deneveLayer('headWorld',1,N);     %Layer to represent world location in head centred coordinates

allocLog(net.ret,nIter,(nSims*N*N))
allocLog(net.eye,nIter,(nSims*N*N))
allocLog(net.hed,nIter,(nSims*N*N))
allocLog(net.hid,nIter,(nSims*N*N))


%Set the input weights for each network layer. Weights are symmetric, such that the
%input weight from neuron A to B is the same as the input weight from B to A
%weights are stored in temporary variables to be used in inputs
wfun = @(ind) Kw.*exp((cos((2*pi/N).*ind)-1)/sigmaw^2); %Anonymous function for the bell-shaped input weights
%Meshgrid which replaced for loops
[a,b,c] = meshgrid(1:N,1:N,1:N);
tempretw = wfun(c-b);
tempeyew = wfun(c-a);
temphedw = wfun(c-a-b-10);

temphidw{1} = wfun(b-a);
temphidw{2} = wfun(b-c);
temphidw{3} = wfun(b-a-c-10);

%Calculate input weights for world.  Communication between the eye and ret
%layers with the world are feedforward from the world to the network
%Weights are based on the circular Von Mises function
%Using class weight function
tempWldRet = vmwfun(wld.ret,K,sigma,v);
tempWldEye = vmwfun(wld.eye,K,sigma,v);
tempWldHed = vmwfun(wld.hed,K,sigma,v);

%Specify the network architecture (who talks to whom?)
%keep order of inputs consistent across network inputs
net.hid.setInput({net.ret,net.eye,net.hed},temphidw);
net.ret.setInput({net.hid,wld.ret},{tempretw,tempWldRet});
net.eye.setInput({net.hid,wld.eye},{tempeyew,tempWldEye});

if headon == 1
    net.hed.setInput({net.hid,wld.hed},{temphedw,tempWldHed});
else
    net.hed.setInput({net.hid},{temphedw});
end
    

%=========== Run the simulation ==============
%nSims = 20;
%Preallocate matrices
[truRet,truEye,truHed] = deal(zeros(N,N,nSims));
[estRet,estEye,estHed] = deal(zeros(N,N,nSims,nIter));
       
for i = 1:N
    disp(num2str(i));
    for j = 1:N
        for s = 1:nSims            
            %Initialise world layers with random delta functions
            r = zeros(1,N);
            realRetPos = i;
            r(realRetPos)= 1;
            wld.ret.initialise(r);
            r = zeros(1,N);
            realEyePos = j;
            r(realEyePos)= 1;
            wld.eye.initialise(r);

            %Switch between delta function and matrix of zeros for head input
            r = zeros(1,N);
            realHedPos = mod((realRetPos + realEyePos - 10),N);   %Change to bais head
            if realHedPos == 0         %Set 0 to 20
                realHedPos = N;
            end

            if headon == 1
                r(realHedPos) = 1;
            else
                r(realHedPos) = 0;
            end
            wld.hed.initialise(r);

            %Vectors logging real positions
            truRet(i,j,s) = realRetPos;
            truEye(i,j,s) = realEyePos;
            truHed(i,j,s) = realHedPos;

            %Reset network layers to zero for each simulation
            net.ret.reset(1,N);
            net.eye.reset(1,N);
            net.hed.reset(1,N);
            net.hid.reset(N,N);

            for t=1:nIter
                %Switch between world and hidden layer inputs for ret, eye and hed
                isFirstTime=t==1;
                setEnabled(net.ret,net.hid.name,~isFirstTime);
                setEnabled(net.ret,wld.ret.name,isFirstTime);
                setEnabled(net.eye,net.hid.name,~isFirstTime);
                setEnabled(net.eye,wld.eye.name,isFirstTime);
                
                if headon == 1
                    setEnabled(net.hed,net.hid.name,~isFirstTime);
                    setEnabled(net.hed,wld.hed.name,isFirstTime);
                end

                %Update for each time point to include recurrent inputs
                normalise = ~isFirstTime;
                net.hid.update(normalise);
                net.ret.update(normalise);
                net.eye.update(normalise);
                net.hed.update(normalise);

                %Add noise to the response of ret and eye networks at the first
                %time point, from the world input
                 if isFirstTime
                     if gainon
                         gainfun(net.hed,10,Kg,sigmag);
                     end
                     addNoise(net.ret);
                     addNoise(net.eye);
                     %Only add noise to head input if it is activated
                     if headon == 1
                         addNoise(net.hed);
                     end
                 end

                 est.Ret = pointEstimate(net.ret);
                 est.Eye = pointEstimate(net.eye);
                 est.Hed = pointEstimate(net.hed);
                 
                 logState(net.ret,'newLog',isFirstTime);
                 logState(net.eye,'newLog',isFirstTime);
                 logState(net.hed,'newLog',isFirstTime);
                 logState(net.hid,'newLog',isFirstTime);

             %============== Plot the simulation==============%
                %Plot vector responses of ret, eye, and hed
                if ploton == 1
                    subplot(2,1,1);
                    cla
                    plotState(net.ret,'linestyle','r-o','linewidth',3);
                    plotState(net.eye,'linestyle','b-o','linewidth',3);
                    plotState(net.hed,'linestyle','g-o','linewidth',3);

                    %Plot 2D matrix response of hid
                    subplot(2,1,2);
                    cla
                    plotState(net.hid);
                    title(num2str(t));
                    
                    %Pauses for plotting
                    if t==1
                        pause(2)
                    else
                        pause(0.15);
                    end  %2./t);
                end

                %Matrices for logging estimate positions
                estRet(i,j,s,t) = pointEstimate(net.ret);
                estEye(i,j,s,t) = pointEstimate(net.eye);
                estHed(i,j,s,t) = pointEstimate(net.hed);
            end
        end
    end
end
%========== Statistical Data ===========%
%Calculate wrapped error
errRet = err(net.ret,estRet(:,:,:,nIter),truRet);
errEye = err(net.eye,estEye(:,:,:,nIter),truEye);
errHed = err(net.hed,estHed(:,:,:,nIter),truHed);

%Outputs
out.ret = {estRet,truRet,errRet};
out.eye = {estEye,truEye,errEye};
out.hed = {estHed,truHed,errHed};


%keyboard;
end


%gain function here