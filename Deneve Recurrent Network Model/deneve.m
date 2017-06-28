function deneve
nIter = 50;
N=20;
Kw = 1;
sigmaw = 0.37;

%Create layers ("population codes"). All layers are instances of the deneveLayer class.
hid = deneveLayer('hidden',N,N);    %Hidden layer
ret = deneveLayer('retinal',1,N);   %Object location on the retina
eye = deneveLayer('eye',1,N);       %Eye position signal
hed = deneveLayer('head',1,N);      %Object location relative to the head

%Specify the network architecture (who talks to whom?)
hid.addInput([ret,eye,hed]);
ret.addInput(hid);
eye.addInput(hid);
hed.addInput(hid);

%Set the input weights for each layer. Weights are symmetric, such that the
%input weight from neuron A to be is the same as the input weight from B to A
wfun = @(ind) Kw.*exp((cos((2*pi/N).*ind)-1)/sigmaw^2); %Anonymous function for the bell-shaped input weights
for j=1:N
    for l=1:N
        for m=1:N
            %Pooling weights for each unit in the input/output layers (i.e. an N x N matrix)
            ret.w{j}{1}(l,m) = wfun(j-l);
            eye.w{j}{1}(l,m) = wfun(j-m);
            hed.w{j}{1}(l,m) = wfun(j-l-m);
            
            %Pooling weights for each hidden unit  (i.e. an 1 x N matrix for each input layer)
            hid.w{l,m}{1}(j) = ret.w{j}{1}(l,m);
            hid.w{l,m}{2}(j) = eye.w{j}{1}(l,m);
            hid.w{l,m}{3}(j) = hed.w{j}{1}(l,m);
        end
    end
end

%=========== Run the simulation ==============
nSims = 20;
for s = 1:nSims
    %Initialise with random numbers
    ret.initialise(poissrnd(50,1,N));
    eye.initialise(poissrnd(50,1,N));
    hed.initialise(zeros(1,N));
    hid.initialise(zeros(N,N));
    eye.r(randi(N))=100;
    ret.r(randi(N))=100;
    for t=1:nIter
        plotState(ret,eye,hed,hid,t);
        
        hid.update();
        ret.update();
        eye.update();
        hed.update();
    end
end
keyboard;

function plotState(ret,eye,hed,hid,t)

subplot(2,1,1);
cla
plot(ret.r,'r-o','linewidth',4);
hold on;
plot(eye.r,'b-o','linewidth',4);
plot(hed.r,'g-o','linewidth',4);

[maxVal,retPos] = max(ret.r);
plot([retPos retPos],[0 maxVal],'r','linewidth',3);
[maxVal,eyePos] = max(eye.r);
plot([eyePos eyePos],[0 maxVal],'b','linewidth',3);
[maxVal,headPos] = max(hed.r);
plot([headPos headPos],[0 maxVal],'g','linewidth',3);
plot(mod([retPos+eyePos retPos+eyePos],ret.nUnits),[0 maxVal],'k:','linewidth',3);
ylim([0,100]);

subplot(2,1,2);
cla
surf(hid.r); zlim([0,10])
if t==1
    pause(3);
else
    pause(0.15);
end%2./t);

