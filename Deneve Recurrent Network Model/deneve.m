function deneve
%Set constants for functions and loops
nIter = 50;
N = 20;
Kw = 1;
sigmaw = 0.37;
K = 20;
v = 1;
sigma = 0.40;

%Create layers ("population codes"). All layers are instances of the deneveLayer class.
%Network layers
net.hid = deneveLayer('hidden',N,N);    %Hidden layer
net.ret = deneveLayer('retinal',1,N);   %Object location on the retina
net.eye = deneveLayer('eye',1,N);       %Eye position signal
net.hed = deneveLayer('head',1,N);      %Object location relative to the head

%World layers
wld.ret = deneveLayer('retinalWorld',1,N);  %Layer to represent world location in retinal coordinates
wld.eye = deneveLayer('eyeWorld',1,N);      %Layer to represent world location in eye centred coordinates
%Add wld.hed layer here

%Set the input weights for each network layer. Weights are symmetric, such that the
%input weight from neuron A to B is the same as the input weight from B to A
%weights are stored in temporary variables to be used in inputs
%Function for setting weights could eventually go in class
%Such that weight function is an argument
wfun = @(ind) Kw.*exp((cos((2*pi/N).*ind)-1)/sigmaw^2); %Anonymous function for the bell-shaped input weights
%Preallocate matrices with dimensions [N,N]
tempretw = zeros(N,N,N);
tempeyew = zeros(N,N,N);
temphedw = zeros(N,N,N);
for j=1:N
    for l=1:N
        for m=1:N
            %Pooling weights for each unit in the input/output layers (i.e. an N x N x N matrix)
            tempretw(l,m,j) = wfun(j-l);
            tempeyew(l,m,j) = wfun(j-m);
            temphedw(l,m,j) = wfun(j-l-m);
            
            %Pooling weights for each hidden unit  (i.e. an N x N x N matrix for each input layer)
            temphidw{1}(j,l,m) = tempretw(l,m,j);
            temphidw{2}(j,l,m) = tempeyew(l,m,j);
            temphidw{3}(j,l,m) = temphedw(l,m,j);
        end
    end
end

%Calculate input weights for world.  Communication between the eye and ret
%layers with the world are feedforward from the world to the network
%Weights are based on the circular Von Mises function
vmfun = @(pos,unit) K.*exp((cos((2*pi/N).*pos-(2*pi/N).*unit)-1)/(sigma^2))+v; %Anon function for world to network weights
%Preallocate matrices to size [N,N]
tempWldRet = zeros(N,N);
tempWldEye = zeros(N,N);
for i=1:N
    for j=1:N
       tempWldRet(i,j)= vmfun(i,j);
       tempWldEye(i,j)= vmfun(i,j);
       %Add tempWorldHed here
    end 
end


%Specify the network architecture (who talks to whom?)
%keep order of inputs consistent across network inputs
net.hid.setInput({net.ret,net.eye,net.hed},temphidw);
net.ret.setInput({net.hid,wld.ret},{tempretw,tempWldRet});
net.eye.setInput({net.hid,wld.eye},{tempeyew,tempWldEye});
net.hed.setInput({net.hid},{temphedw});    %Add wld.hed input here

%=========== Run the simulation ==============
nSims = 20;
%Preallocating vectors for true and estimate locations
for s = 1:nSims
    %Initialise world layers with random delta functions
    r = zeros(1,N);
    realRetPos = randi(N);
    r(realRetPos)= 1;
    wld.ret.initialise(r);
    r = zeros(1,N);
    realEyePos = randi(N);
    r(realEyePos)= 1;
    wld.eye.initialise(r);
    %Add wld.hedinitialise and wld.hed.resp here
    
    %Record real positions in vectors
    truRet(s) = realRetPos;
    truEye(s) = realEyePos;
    
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
        %Add hed.input here
        
        %Update for each time point to include recurrent inputs
        net.hid.update(~isFirstTime);
        net.ret.update(~isFirstTime);
        net.eye.update(~isFirstTime);
        net.hed.update(~isFirstTime);
        
        %Add noise to the response of ret and eye networks at the first
        %time point, from the world input
         if isFirstTime
             addNoise(net.ret,N);
             addNoise(net.eye,N);
             %Add addNoise(net.hed,i)
         end
         
         est.Ret = pointEstimate(net.ret);
         est.Eye = pointEstimate(net.eye);
         est.Hed = pointEstimate(net.hed);
        
        %Plot the simulation
        plotState(net,t,est);
    end
    
    %Vectors of peak estimates
    estRetLog(s) = est.Ret;
    estEyeLog(s) = est.Eye;
    estHedLog(s) = est.Hed;
end
keyboard;

%====Function to plot the simulation====

function plotState(net,t,est)
%add a wld input eventually

%Plot ret, eye, and hed responses on the same subplot
subplot(2,1,1);
cla
plot(net.ret.resp,'r-o','linewidth',4);
hold on;
plot(net.eye.resp,'b-o','linewidth',4);
plot(net.hed.resp,'g-o','linewidth',4);

%Add line to indicate peak of ret, eye, and hed networks
%Change to point estimate function (within class) from supp material
retPos = est.Ret;
maxVal = max(net.ret.resp);
plot([retPos retPos],[0 maxVal],'r','linewidth',3);
eyePos = est.Eye;
maxVal = max(net.eye.resp);
plot([eyePos eyePos],[0 maxVal],'b','linewidth',3);
headPos = est.Hed;
maxVal = max(net.hed.resp);
plot([headPos headPos],[0 maxVal],'g','linewidth',3);
%Add line to indicate point satisfying hed=ret+eye rule
plot(mod([retPos+eyePos retPos+eyePos],net.ret.nUnits),[0 maxVal],'k:','linewidth',3);
ylim([0,100]);

%Plot hidden layer on seperate subplot
subplot(2,1,2);
cla
surf(net.hid.resp); zlim([0,10])

%Timing of simulation
if t==1
    pause(2);
else
    pause(0.15);
end%2./t);

