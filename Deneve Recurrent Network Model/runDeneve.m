function runDeneve
%Set the number of simulations
nSims = 10;
%Head input switch, 1 for on, 0 for off
head = 0;
%Plot simulation switch, 1 for on, 0 for off
plot = 0;

%Run the deneve function
[net,wld,out] = deneve(nSims,head,plot);

%Calculate mean and standard deviation (in matrices of [Ret,Eye])
mret = x2rad(net.ret,out.ret{3});
meye = x2rad(net.eye,out.eye{3});
mhed = x2rad(net.hed,out.hed{3});

radmean.ret = circ_mean(mret,[],3);
radmean.eye = circ_mean(meye,[],3);
radmean.hed = circ_mean(mhed,[],3);

mean.ret = rad2x(net.ret,radmean.ret);
mean.eye = rad2x(net.eye,radmean.eye);
mean.hed = rad2x(net.hed,radmean.hed);

% radMeanErr = circ_mean([mret,meye,mhed]);
% meanErr = rad2x(net.ret,radMeanErr);

radstd.ret = circ_std(mret,[],[],3);
radstd.eye = circ_std(meye,[],[],3);
radstd.hed = circ_std(mhed,[],[],3);

std.ret = rad2x(net.ret,radstd.ret);
std.eye = rad2x(net.eye,radstd.eye);
std.hed = rad2x(net.hed,radstd.hed);

% radStdErr = circ_std([mret,meye,mhed]);
% stdErr = rad2x(net.ret,radStdErr);


%Plot heatmap of means and standard deviations
figure;
subplot(2,3,1);
imagesc(mean.ret);
title('Ret Mean Error');
subplot(2,3,2);
imagesc(mean.eye);
title('Eye Mean Error');
subplot(2,3,3);
imagesc(mean.hed);
title('Head Mean Error');
subplot(2,3,4);
imagesc(std.ret);
title('Ret Std Dev');
xlabel('eye position');
ylabel('retinal position');
subplot(2,3,5);
imagesc(std.eye);
title('Eye Std Dev');
subplot(2,3,6);
imagesc(std.hed);
title('Head Std Dev');

%Plot histogram of errors in a seperate figure (no bin size set as yet)
figure;
histogram(mean.ret(:),'binwidth',0.1);
hold on
histogram(mean.eye(:),'binwidth',0.1);
histogram(mean.hed(:),'binwidth',0.1);

keyboard

