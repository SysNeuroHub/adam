function runDeneve
%Set the number of simulations
nSims = 5;
%Head input switch, 1 for on, 0 for off
head = 1;
%Plot simulation switch, 1 for on, 0 for off
plot = 0;

%Run the deneve function

[netg,wldg,outg] = deneve(nSims,head,plot,1);
%[netn,wldn,outn] = deneve(nSims,head,plot,0);


%Calculate mean and standard deviation (in matrices of [Ret,Eye])
mret = x2rad(netg.ret,outg.ret{3});
meye = x2rad(netg.eye,outg.eye{3});
mhed = x2rad(netg.hed,outg.hed{3});

radmean.ret = circ_mean(mret,[],3);
radmean.eye = circ_mean(meye,[],3);
radmean.hed = circ_mean(mhed,[],3);

mean.ret = rad2x(netg.ret,radmean.ret);
mean.eye = rad2x(netg.eye,radmean.eye);
mean.hed = rad2x(netg.hed,radmean.hed);

% radMeanErr = circ_mean([mret,meye,mhed]);
% meanErr = rad2x(net.ret,radMeanErr);

radstd.ret = circ_std(mret,[],[],3);
radstd.eye = circ_std(meye,[],[],3);
radstd.hed = circ_std(mhed,[],[],3);

std.ret = rad2x(netg.ret,radstd.ret);
std.eye = rad2x(netg.eye,radstd.eye);
std.hed = rad2x(netg.hed,radstd.hed);

%Plot heatmap of means and standard deviations
clim= minmax(reshape([mean.ret,mean.eye,mean.hed],1,[]));
figure;
subplot(2,3,1);
imagesc(mean.ret,clim);
title('Ret Mean Error');
subplot(2,3,2);
imagesc(mean.eye,clim);
title('Eye Mean Error');
subplot(2,3,3);
imagesc(mean.hed,clim);
title('Head Mean Error');

clim = minmax(reshape([std.ret,std.eye,std.hed],1,[]));
subplot(2,3,4);
imagesc(std.ret,clim);
title('Ret Std Dev');
xlabel('eye position');
ylabel('retinal position');
subplot(2,3,5);
imagesc(std.eye,clim);
title('Eye Std Dev');
subplot(2,3,6);
imagesc(std.hed,clim);
title('Head Std Dev');

keyboard

% %Plot histogram of errors in a seperate figure (no bin size set as yet)
% figure;
% h = histogram(mean.ret(:),30);
% hold on
% histogram(mean.eye(:),h.BinEdges);
% histogram(mean.hed(:),h.BinEdges);
% 
% keyboard


%Plotting by time point
% nret = x2rad(net.ret,norm.ret.err);
% neye = x2rad(net.eye,norm.eye.err);
% nhed = x2rad(net.hed,norm.hed.err);
% 
% gret = x2rad(net.ret,gain.ret.err);
% geye = x2rad(net.eye,gain.eye.err);
% ghed = x2rad(net.hed,gain.hed.err);
% 
% radstd.nret = circ_std(nret);
% radstd.neye = circ_std(neye);
% radstd.nhed = circ_std(nhed);
% radstd.gret = circ_std(gret);
% radstd.geye = circ_std(geye);
% radstd.ghed = circ_std(ghed);
% 
% stdnret = rad2x(net.ret,radstd.nret);
% stdneye = rad2x(net.eye,radstd.neye);
% stdnhed = rad2x(net.hed,radstd.nhed);
% stdgret = rad2x(net.ret,radstd.gret);
% stdgeye = rad2x(net.eye,radstd.geye);
% stdghed = rad2x(net.hed,radstd.ghed);
% ylimit = [0,0.7];
% 
% figure;
% plot(stdnret,'b*-');
% title('Retinal Layer');
% hold on
% plot(stdgret,'r*-');
% xlabel('Time (t)');
% ylabel('Standard Deviation of Error');
% legend('without gain','with gain');
% ylim(ylimit);
% 
% figure;
% plot(stdneye,'b*-');
% title('Eye Layer');
% hold on
% plot(stdgeye,'r*-');
% xlabel('Time (t)');
% ylabel('Standard Deviation of Error');
% legend('without gain','with gain');
% ylim(ylimit);
% 
% figure;
% plot(stdnhed,'b*-');
% title('Head Layer');
% hold on
% plot(stdghed,'r*-');
% xlabel('Time (t)');
% ylabel('Standard Deviation of Error');
% legend('without gain','with gain');
% ylim(ylimit);


% % radStdErr = circ_std([mret,meye,mhed]);
% % stdErr = rad2x(net.ret,radStdErr);
% 
% 

% 
% keyboard

