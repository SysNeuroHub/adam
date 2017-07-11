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
vmfun = @(pos,unit) K.*exp((cos((pos-unit).*(2*pi/N))-1)/(sigma^2))+v; %Anon function for world to network weights
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
    
    %Vectors logging real positions
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
             addNoise(net.ret);
             addNoise(net.eye);
             %Add addNoise(net.hed)
         end
         
         est.Ret = pointEstimate(net.ret);
         est.Eye = pointEstimate(net.eye);
         est.Hed = pointEstimate(net.hed);
        
     %============== Plot the simulation==============%
        %Plot vector responses of ret, eye, and hed
        subplot(2,1,1);
        cla
        plotState(net.ret,t,'linestyle','r-o','linewidth',4);
        plotState(net.eye,t,'linestyle','b-o','linewidth',4);
        plotState(net.hed,t,'linestyle','g-o','linewidth',4);
        
        %Plot 2D matrix response of hid
        subplot(2,1,2);
        cla
        plotState(net.hid,t);
    end
    
    %Vectors logging estimate positions
    estRet(s) = pointEstimate(net.ret);
    estEye(s) = pointEstimate(net.eye);
    estHed(s) = pointEstimate(net.hed);
end

%========== Statistical Data ===========%
%Calculate wrapped error
errRet = err(estRet,truRet);
errEye = err(estEye,truEye);

%Calculate mean and standard deviation (in matrices of [Ret,Eye])
meanErr = [mean(errRet),mean(errEye)];
stdErr = [std(errRet),std(errEye)];

%Plot histogram of errors in a seperate figure (no bin size set as yet)
figure;
histogram(errRet);
hold on
histogram(errEye);

keyboard;
end


%=== Currently having an issue putting these in the class===%
function error = err(est,true)
    %Function to determine absolute error
    estRad = x2rad(est);
    truRad = x2rad(true);

    errRad = atan2(sin(estRad-truRad),cos(estRad-truRad));
    error = rad2x(errRad);
end

function rad = x2rad(x)
    %Functon to convert 1-20 scale to radians - change 20 to variable
    rad = x.*(2*pi)./20;
end

function x = rad2x(rad)
    %Functon to convert radians to 1-20 scale - change 20 to variable
    x = rad.*20./(2*pi);
end

