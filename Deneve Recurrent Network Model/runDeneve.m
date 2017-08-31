function runDeneve

%Run the deneve function
[n,est] = deneve('nSims',20,'headWorldOn',true,'plotIt',false,'suppressLayer','retinal','addNoise',true);

%Calculate mean and standard deviation (in matrices of [Ret,Eye])
mret = x2rad(n.retinal,est.ret{3});
meye = x2rad(n.eye,est.eye{3});
mhed = x2rad(n.head,est.hed{3});

radmean.retinal = circ_mean(mret,[],3);
radmean.eye = circ_mean(meye,[],3);
radmean.head = circ_mean(mhed,[],3);

meanEst.ret = rad2x(n.retinal,radmean.retinal);
meanEst.eye = rad2x(n.eye,radmean.eye);
meanEst.hed = rad2x(n.head,radmean.head);

% radMeanErr = circ_mean([mret,meye,mhed]);
% meanErr = rad2x(net.ret,radMeanErr);

radstd.ret = circ_std(mret,[],[],3);
radstd.eye = circ_std(meye,[],[],3);
radstd.hed = circ_std(mhed,[],[],3);

stdEst.ret = rad2x(n.retinal,radstd.ret);
stdEst.eye = rad2x(n.eye,radstd.eye);
stdEst.hed = rad2x(n.head,radstd.hed);

%Plot heatmap of means and standard deviations
clim= minmax(reshape([meanEst.ret,meanEst.eye,meanEst.hed],1,[]));
figure;
subplot(2,3,1);
imagesc(meanEst.ret,clim);
title('Ret Mean Error');
subplot(2,3,2);
imagesc(meanEst.eye,clim);
title('Eye Mean Error');
subplot(2,3,3);
imagesc(meanEst.hed,clim);
title('Head Mean Error');

clim = minmax(reshape([stdEst.ret,stdEst.eye,stdEst.hed],1,[]));
subplot(2,3,4);
imagesc(stdEst.ret,clim);
title('Ret Std Dev');
xlabel('eye position');
ylabel('retinal position');
subplot(2,3,5);
imagesc(stdEst.eye,clim);
title('Eye Std Dev');
subplot(2,3,6);
imagesc(stdEst.hed,clim);
title('Head Std Dev');

keyboard

% %Plot histogram of errors in a seperate figure (no bin size set as yet)
% figure;
% h = histogram(mean.retinal(:),30);
% hold on
% histogram(mean.eye(:),h.BinEdges);
% histogram(mean.head(:),h.BinEdges);
% 
% keyboard


%Plotting by time point
% nret = x2rad(net.ret,norm.ret.err);
% neye = x2rad(net.eye,norm.eye.err);
% nhed = x2rad(net.hed,norm.hed.err);
% 
% gret = x2rad(net.ret,gain.retinal.err);
% geye = x2rad(net.eye,gain.eye.err);
% ghed = x2rad(net.hed,gain.head.err);
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

