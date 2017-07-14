%Set the number of simulations
nSims = 20;
%Head input switch, 1 for on, 0 for off
head = 0;
%Plot simulation switch, 1 for on, 0 for off
plot = 0;

%Run the deneve function
[net,wld,ret,eye,hed] = deneve(nSims,head,plot);

%Calculate mean and standard deviation (in matrices of [Ret,Eye])
mret = x2rad(net.ret,ret(:,3));
meye = x2rad(net.eye,eye(:,3));
mhed = x2rad(net.hed,hed(:,3));

radMeanErr = circ_mean([mret,meye,mhed]);
meanErr = rad2x(net.ret,radMeanErr);

radStdErr = circ_std([mret,meye,mhed]);
stdErr = rad2x(net.ret,radStdErr);

%Plot histogram of errors in a seperate figure (no bin size set as yet)
figure;
histogram(ret(:,3),'binwidth',0.1);
hold on
histogram(eye(:,3),'binwidth',0.1);
histogram(hed(:,3),'binwidth',0.1);

