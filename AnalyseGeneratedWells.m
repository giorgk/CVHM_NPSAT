%% Load generated wells
load('Generated_wells.mat','wellGen');
%% Compute density Map
WellGenDensity = mvksdensity(wellGen(:,1:2),XYbas,'bandwidth',6000);
%%
clf
u = 1;
trisurf(tri, XYbas(:,1), XYbas(:,2), WellGenDensity/max(WellGenDensity), 'edgecolor', 'none')
%hold on
%plot(WelldataXY(:,1), WelldataXY(:,2),'.')
view(0,90);
axis equal
colorbar
alpha(1);
axis off
%% Pumping
histogram(wellGen(:,5)/5.450992969400)
xlim([0 10000])
xlabel('Pumping rate [gpm]')

%% Screen length
histogram((wellGen(:,3)-wellGen(:,4))*3.28084)
xlabel('Screen length [ft]')