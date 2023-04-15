%% settings

warning off
fitrange=1:200; % range for linear fit
A=5.5e-3*0.2e-3; % probe area
smoothPoints=25; % # of points for smooth(I)
smoothPointsdI=10; % # of points for smooth (dI)
expstart=1; % value in V where to start exponential fit
expoffset=5; % offset in V to left of plasma potential for exponential fitting

%%

U=ken(:,1);
I=ken(:,2);

% smoothed I for derivative
ISmooth=smooth(U,I,smoothPoints);
% get floating potential
params.floatingPot=U(find(ISmooth>0,1,'first'));

% derive characteristic
dI=gradient(ISmooth,U);
dISmooth=smooth(dI,smoothPointsdI);
[~,params.plasmaPotIdx]=max(dISmooth);
params.plasmaPot=U(params.plasmaPotIdx);

% linear fit to isat region
linfit=ezfit(U(fitrange),I(fitrange)*1e5,'lin',mean(dI(10:100))*1e5,I(fitrange(1))*1e5);
linfit.a=linfit.a/1e5;
linfit.c=linfit.c/1e5;

% plot original data
figure(1)
clf
subplot(2,2,1)
plot(U,I)
hold on
plot(U,ISmooth,'r')
plot(linfit)
prettyplot('U [V]','I [A]')



% plot derivative
subplot(2,2,2)
plot(U,dI)
hold on
plot(U,dISmooth,'r')
plot(U(params.plasmaPotIdx),dI(params.plasmaPotIdx),'o')
hold off
prettyplot('U [V]','dI','dI/dU')
legend('char','smoothed','\Phi_p','Location','Northwest')

% calculate I_isat
params.Iisat=linfit(params.plasmaPot);

% update original plot with Plasma potential
subplot(2,2,1)
plot(params.plasmaPot,I(params.plasmaPotIdx),'o')
plot(params.plasmaPot,params.Iisat,'ro')
legend('char','smoothed','linfit','\Phi_p','I_{isat}','Location','Northwest')
hold off

% calculate electron current
Ie=I-linfit(U);
params.Iesat=Ie(params.plasmaPotIdx);
fitrangeExp=find(U>expstart,1):find(U<params(1).plasmaPot-expoffset,1,'last');
expfit=ezfit(U(fitrangeExp),Ie(fitrangeExp)*1e5,'exp',1,1/3);
expfit.a=expfit.a/1e5;
params.Te=1/expfit.b;
[params.ne,params.ni]=density(-params.Iesat,-params.Iisat,A,params.Te);

% plot electron current and exp fit
subplot(2,2,3)
semilogy(U,Ie)
hold on
fitted=expfit(U);
semilogy(U,fitted,'r')
semilogy(U(fitrangeExp),Ie(fitrangeExp),'m+')
semilogy(U(params.plasmaPotIdx),Ie(params.plasmaPotIdx),'o')
axis([-Inf Inf 1e-8 max(Ie)])
prettyplot('U [V]','log(Ie)','Electron current')
legend('Ie','fit points','fit','Location','Northwest')
hold off
%%
subplot(2,2,4)
axis off
text(0.3,0.3,evalc('params'),'FontSize',16)
%%
disp(params)


drawnow
