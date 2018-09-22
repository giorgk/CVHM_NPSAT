%
%% Find out the total amount water from groundwater recharge
rch = read_Scattered( 'CVHM_RCH.npsat', 2);
% create interpolant
Frch = scatteredInterpolant(rch.p(:,1), rch.p(:,2), rch.v);
%% Read mesh
[msh_nd, msh_el] = readMeshfile('CVHM_mesh.npsat');
msh_el = msh_el + 1;
%
%% Find the an approximation of the amount of recharge each element is 
% going to receive.
TOT_RCH = 0;
for ii = 1:size(msh_el, 1)
    xp = msh_nd(msh_el(ii,:),1);
    yp = msh_nd(msh_el(ii,:),2);
    xc = mean(xp);
    yc = mean(yp);
    A = polyarea(xp, yp);
    
    Vc = Frch(xc, yc);
    TOT_RCH = TOT_RCH + Vc * A;
end
%% Find the total amount from streams
STRM = readStreams('CVHM_streams.npsat');
TOTSTRM = sum(([STRM.Q]').*([STRM.area]'));
%% Total volume to be pumped by wells
% 28988263.07 [m^3/day]
fprintf('%.2f\n',TOT_RCH + TOTSTRM)
%% Read the well database from Rich Pauloo. 
% At this point the well database should not be distributed, therefore I do
% not include it in the CVHM repository, but I read it from my local folder
% ~/Documents/UCDAVIS/CVHM_DATA/WellData
%WellTable = readtable('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/giorgos.csv');
%WellTable = readtable('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/RichWellDataV2/dfspcv.csv');
%% create a shapefile from the table
for ii = 1:size(WellTable,1)
    S(ii,1).Geometry = 'Point';
    S(ii,1).X = WellTable.lon(ii);
    S(ii,1).Y = WellTable.lat(ii);
    temp = WellTable.WCRNumber(ii);
    S(ii,1).WCRNumber = temp{1,1};
    temp = WellTable.type(ii);
    S(ii,1).type = temp{1,1};
    S(ii,1).Top = WellTable.top(ii);
    S(ii,1).Bot = WellTable.bot(ii);
    S(ii,1).depth = WellTable.depth(ii);
    S(ii,1).Q = WellTable.Q(ii);
end
%% OR since Rich created a shapefile why not read it directly
% bot = depth from land surface to the bottom of the perforated interval [ft]
% top = depth from land surface to the top of the perforated interval [ft]
% bot - top is the length of the screened interval
% Q pumping in gpm
% year of construction

S = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/CV_only_PubAg.shp');
%% Write shapefile. The coordinate system is EPSG:4326
shapewrite(S, '/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/CVwelldata');
%% Create a shapefile with the well data only for the CVHM area.
S = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/CVonly_WellData.shp');
%% convert well type to id
% 1 -> 'agriculture'
% 2 -> 'public'
% 3 -> 'domestic'
for ii = 1:length(S)
    if strcmp(S(ii,1).type, 'agriculture')
        S(ii,1).Cat = 1;
    elseif strcmp(S(ii,1).type, 'public')
        S(ii,1).Cat = 2;
    elseif strcmp(S(ii,1).type, 'domestic')
        S(ii,1).Cat = 3;
    else
        S(ii,1).Cat = 4;
    end 
end
%% Make a separate Sag for the agriculture and public wells
Sag = S;
Sag([Sag.Cat]' == 3,:) = [];

%% Add the CVHM Surface elevation 
% (No longer needed since the Bot corresponds to Depth and bot - top is the
% screen length)
[ ELEV, info ] = readArcGisASCIIfile( 'gis_data/PP1766_DIS_Layer1_Top_ASCII.txt' );
Xgrid = info.xllcorner+info.cellsize/2:...
                info.cellsize:...
                info.xllcorner+info.cellsize/2+info.cellsize*(info.ncols-1);
    
Ygrid = info.yllcorner+info.cellsize/2:...
        info.cellsize:...
        info.yllcorner+info.cellsize/2+info.cellsize*(info.nrows-1);

Felev = griddedInterpolant({Xgrid, Ygrid}, ELEV');
%% Interpolate Elevation to well nodes
% Wellelev = Felev(WelldataXY(:,1), WelldataXY(:,2))*0.3048; %convert to mm
%% 
WelldataXY = [[Sag.X]' [Sag.Y]'];
WelldataQ = [Sag.Q]'*5.450992969; %convert to m^3/day
WelldataDepth = [Sag.bot]'*0.3048;
WelldataScreenLen = ([Sag.bot]' - [Sag.top]')*0.3048;
%% Bring the basic shapefile
% Run the first 3 sections of the PrepareGeometryData.m script
% make a list of the barycenters of the bas shapefile
XYbas = nan(length(bas),2);
for ii = 1:length(bas)
    XYbas(ii,:) = [mean(bas(ii,1).X(1:4)) mean(bas(ii,1).Y(1:4))];
end
%% prepare tri for plot
tri = delaunay(XYbas(:,1),XYbas(:,2));
outline = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_NPSAT/gis_data/CVHM_Mesh_outline.shp');
for ii = 1:size(tri,1)
    XYtri(ii,:) = [mean(XYbas(tri(ii,:),1)) mean(XYbas(tri(ii,:),2))];
end
% remove the triangles outside the CVHM outline
[Xout,Yout] = polysplit(outline.X, outline.Y);
in = inpolygon(XYtri(:,1), XYtri(:,2),Xout{1,1}, Yout{1,1});
tri(~in,:) = [];

%% Calculate the KDE for the density of wells
Welldensity = mvksdensity(WelldataXY,XYbas,'bandwidth',3000);

Welldensity(Welldensity==0) = min(Welldensity(Welldensity~=0));
logWellDensity = log10(Welldensity);
logWellDensity_pos = logWellDensity + abs(min(logWellDensity));

maxV = max(logWellDensity_pos);
minV = min(logWellDensity_pos);
Vs   = (logWellDensity_pos - minV) / (maxV - minV);
FwellDens = scatteredInterpolant(XYbas(:,1), XYbas(:,2), Vs);
FwellDens_nrm = scatteredInterpolant(XYbas(:,1), XYbas(:,2), Welldensity/max(Welldensity));
%%  plot density estimation
clf
u = 1;
trisurf(tri, XYbas(:,1), XYbas(:,2), Vs*u + (1-u)*(Welldensity/max(Welldensity)), 'edgecolor', 'none')
%hold on
%plot(WelldataXY(:,1), WelldataXY(:,2),'.')
view(0,90);
axis equal
colorbar
alpha(1);
axis off
%
%% Try distributions to data
% xtest = WelldataQ;
% xtest(xtest<=0,:) = [];
% xtest = log(xtest);
% pd = fitdist(xtest,'Kernel','Kernel','epanechnikov','Bandwidth',0.1); % 
% xv = min(xtest):max(xtest);
% yv = pdf(pd,xv);
% clf
% histogram(xtest,'Normalization','pdf')
% hold on
% plot(xv,yv, 'r','LineWidth',2)
%% Compute dependent spatial statistics between Pumping Depth and Screen length
% find wells with pumping rates
clear WellQ WellD WellS
id_incQ = WelldataQ > 0; % find the records with pumpinh rates
logQ = log10(WelldataQ(id_incQ)); % isolate and log transform those 
XY_Q = WelldataXY(id_incQ,:); % extract the XY of those that have Q
depth = WelldataDepth(id_incQ); % extract the Depth of those that have Q
SL = WelldataScreenLen(id_incQ); % extract the Screen lenght (SL) of those that have Q

% From the wells that have Q find those that have Depth
id_incD = depth > 0;
logD = log10(depth(id_incD)); % extract the depth
XY_QD = XY_Q(id_incD,:); % extract the XY of the records that have Q and D
logQ_D = logQ(id_incD); % extract the Q of the records that have Q and D
SL = SL(id_incD); % extract the SL of the records that have Q and D

% From the wells that have Q and D find those that have SL
id_incS = SL > 0;
logS = log10(SL(id_incS)); % extract the SL
XY_QDS = XY_QD(id_incS,:);  % extract the XY of the records that have Q, D and SL
logD_S = logD(id_incS); % extract the D of the records that have Q, D and SL
logQ_DS = logQ_D(id_incS); % extract the Q of the records that have Q, D and SL
[xg, yg] = meshgrid(0:0.2:5,0:0.2:3);
xgr = reshape(xg,size(xg,1)*size(xg,2),1);
ygr = reshape(yg,size(yg,1)*size(yg,2),1);

%
Thres_init = 5000;
N = size(XYbas,1);
cnt_neg = 0;
for jj = 1:size(XYbas,1)
    [jj N]
    % PUMPING
    dst = sqrt((XYbas(jj,1) - XY_Q(:,1)).^2 + (XYbas(jj,2) - XY_Q(:,2)).^2);
    w = exp(-(1/Thres_init)*dst );
    WellQ.mean(jj,1) = sum(logQ.*(w/sum(w)));
    WellQ.std(jj,1) = std(logQ,w);
    
    % DEPTH
    dst = sqrt((XYbas(jj,1) - XY_QD(:,1)).^2 + (XYbas(jj,2) - XY_QD(:,2)).^2);
    Thres = Thres_init;
    while 1
        w = exp(-(1/Thres)*dst );
        D_Q_fit = fit(logQ_D,logD, 'poly1', 'Weights',w);
        ab = coeffvalues(D_Q_fit);
        if ab(1) < 0
            Thres = Thres + 5000;
        else
            break;
        end
        if Thres > 3000000
            cnt_neg = cnt_neg + 1;
            break;
        end
    end
    ci = predint(D_Q_fit,[0:0.1:4]');
    tempfit=fit([0:0.1:4]',ci(:,1), 'poly1');
    temptab = coeffvalues(tempfit);
    ab95 = confint(D_Q_fit);
    WellD.a(jj,1) = ab(1);
    WellD.b(jj,1) = ab(2);
    WellD.a95(jj,1) = ab(1) - ab95(1,1);
    WellD.b95(jj,1) = ab(2) - ab95(1,2);
    WellD.a_p95(jj,1) = temptab(1) - ab(1);
    WellD.b_p95(jj,1) = temptab(2) - ab(2);
    WellD.Thres(jj,1) = Thres;
    
    % SCREEN LENGTH
    dst = sqrt((XYbas(jj,1) - XY_QDS(:,1)).^2 + (XYbas(jj,2) - XY_QDS(:,2)).^2);
    Thres = Thres_init;
    while 1
        w = exp(-(1/Thres)*dst );
        SL_QD_fit = fit([logQ_DS, logD_S],logS, 'poly11', 'Weights',w);
        abc = coeffvalues(SL_QD_fit);
        if abc(2) < 0 || abc(3) < 0
            Thres = Thres + 5000;
        else
            break;
        end
        if Thres > 3000000
            cnt_neg = cnt_neg + 1;
            break;
        end
    end
    ci = predint(SL_QD_fit,[xgr, ygr]);
    tempfit=fit([xgr ygr],ci(:,1), 'poly11');
    temptab = coeffvalues(tempfit);
    
    
    abc95 = confint(SL_QD_fit);
    WellS.Cx(jj,1) = abc(2);
    WellS.Cy(jj,1) = abc(3);
    WellS.C0(jj,1) = abc(1);
    WellS.Cx95(jj,1) = abc(2) - abc95(1,2);
    WellS.Cy95(jj,1) = abc(3) - abc95(1,3);
    WellS.C095(jj,1) = abc(1) - abc95(1,1);
    WellS.Cx_p95(jj,1) = temptab(2) - abc(2);
    WellS.Cy_p95(jj,1) = temptab(3) - abc(3);
    WellS.C0_p95(jj,1) = temptab(1) - abc(1);
    WellS.Thres(jj,1) = Thres;
end
%% Save the well statics
save('Well_param_stat','WellQ', 'WellD', 'WellS');
%% Create interpolation functions
WellQ.Fmean = scatteredInterpolant(XYbas(:,1), XYbas(:,2), WellQ.mean);
WellQ.Fstd = scatteredInterpolant(XYbas(:,1), XYbas(:,2), WellQ.std);

WellD.Fa = scatteredInterpolant(XYbas(:,1), XYbas(:,2), WellD.a);
WellD.Fb = scatteredInterpolant(XYbas(:,1), XYbas(:,2), WellD.b);
WellD.Fa95 = scatteredInterpolant(XYbas(:,1), XYbas(:,2), WellD.a95);
WellD.Fb95 = scatteredInterpolant(XYbas(:,1), XYbas(:,2), WellD.b95);

WellS.FC0 = scatteredInterpolant(XYbas(:,1), XYbas(:,2), WellS.C0);
WellS.FCx = scatteredInterpolant(XYbas(:,1), XYbas(:,2), WellS.Cx);
WellS.FCy = scatteredInterpolant(XYbas(:,1), XYbas(:,2), WellS.Cy);
WellS.FC095 = scatteredInterpolant(XYbas(:,1), XYbas(:,2), WellS.C095);
WellS.FCx95 = scatteredInterpolant(XYbas(:,1), XYbas(:,2), WellS.Cx95);
WellS.FCy95 = scatteredInterpolant(XYbas(:,1), XYbas(:,2), WellS.Cy95);

%% Plot parameters
figure(1);trisurf(tri, XYbas(:,1), XYbas(:,2), WellD.Fb(XYbas(:,1), XYbas(:,2)), 'edgecolor', 'none')
view(0,90);
axis equal
axis off
colorbar

%% Compute independent spatial statistics for each bas point
% % % % % For the pumping find the weighted mean and weighted standard deviation
% % % % % only nonzero positive Q values will be considered
% % % % id_incQ = WelldataQ > 0;
% % % % logQ = log10(WelldataQ(id_incQ));
% % % % 
% % % % id_incSL = WelldataScreenLen > 0;
% % % % logSL = log10(WelldataScreenLen(id_incSL));% log screen length
% % % % logD = log10(WelldataDepth); % log depth
% % % % Thres = 5000;
% % % % for jj = 1:size(XYbas,1)
% % % %     
% % % %     dst = sqrt((XYbas(jj,1) - WelldataXY(:,1)).^2 + (XYbas(jj,2) - WelldataXY(:,2)).^2);
% % % %     w = exp(-(1/Thres)*dst );
% % % %     
% % % %     wQ = w(id_incQ);
% % % %     wSL = w(id_incSL);
% % % %     
% % % %     WellQmean(jj,1) = sum(logQ.*(wQ/sum(wQ)));
% % % %     WellQstd(jj,1) = std(logQ,wQ);
% % % %     
% % % %     WellSLmean(jj,1) = sum(logSL.*(wSL/sum(wSL)));
% % % %     WellSLstd(jj,1) = std(logSL,wSL);
% % % %     
% % % %     WellDmean(jj,1) = sum(logD.*(w/sum(w)));
% % % %     WellDstd(jj,1) = std(logD,w);
% % % %     
% % % % end
%% Create interpolation functions
% % % % FQmean = scatteredInterpolant(XYbas(:,1), XYbas(:,2), WellQmean);
% % % % FQstd = scatteredInterpolant(XYbas(:,1), XYbas(:,2), WellQstd);
% % % % FSLmean = scatteredInterpolant(XYbas(:,1), XYbas(:,2), WellSLmean);
% % % % FSLstd = scatteredInterpolant(XYbas(:,1), XYbas(:,2), WellSLstd);
% % % % FDmean = scatteredInterpolant(XYbas(:,1), XYbas(:,2), WellDmean);
% % % % FDstd = scatteredInterpolant(XYbas(:,1), XYbas(:,2), WellDstd);
%% plot spatial statistics
tp = 1;
stp = {'Q','SL','D'};
clf
subplot(1,2,1)
eval(['temp = Well' stp{tp} 'mean;']);
trisurf(tri, XYbas(:,1), XYbas(:,2), temp, 'edgecolor', 'none')
view(0,90);
axis equal
axis off
title ([stp{tp} ' mean']);
colorbar
subplot(1,2,2)
eval(['temp = Well' stp{tp} 'std;']);
trisurf(tri, XYbas(:,1), XYbas(:,2), temp, 'edgecolor', 'none')
view(0,90);
axis equal
axis off
title ([stp{tp} ' std']);
colorbar
%}
%% Prepare streams for faster queries
%{
STRM = readStreams('CVHM_streams.npsat');
% add bounding box info
STRM_bbox = nan(length(STRM),4);
for ii = 1:length(STRM)
    STRM_bbox(ii,:) = [min(STRM(ii,1).poly(:,1)) max(STRM(ii,1).poly(:,1)) ...
                       min(STRM(ii,1).poly(:,2)) max(STRM(ii,1).poly(:,2))];
end
%}
%% ============== Main Algorithm ===============
% 10 ft 20ft/day for 100gpm
%% load/process required data
warea = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_NPSAT/gis_data/CVHM_Mesh_inline_buf400.shp');
[wareaX, wareaY] = polysplit(warea.X, warea.Y);
wareaXmin = 99999999999;
wareaXmax = -99999999999;
wareaYmin = 99999999999;
wareaYmax = -99999999999;
for ii = 1:length(wareaX)
    if min(wareaX{ii,1}) < wareaXmin
        wareaXmin = min(wareaX{ii,1});
    end
    if max(wareaX{ii,1}) > wareaXmax
        wareaXmax = max(wareaX{ii,1});
    end
    if min(wareaY{ii,1}) < wareaYmin
        wareaYmin = min(wareaY{ii,1});
    end
    if max(wareaY{ii,1}) > wareaYmax
        wareaYmax = max(wareaY{ii,1});
    end
end
% read the initial estimate of the water table
topelev = read_Scattered('CVHM_top_elev.npsat',2);
Felev = scatteredInterpolant(topelev.p(:,1),topelev.p(:,2),topelev.v);

% abd the bottom of the aquifer
botelev = read_Scattered('CVHM_Bot_elev.npsat',2);
Fbot = scatteredInterpolant(botelev.p(:,1),botelev.p(:,2),botelev.v);

%
%% main RUN
TOT_WELL_Q = 0;
wellGen = nan(100000,5); % X Y T B Q
cnt_well = 1;
clf
plot(wareaX{1,1}, wareaY{1,1},'r')
hold on
udens = 1;
while TOT_WELL_Q < 28988263.07
    % generate a random point inside the CV bounding box until the point is
    % within CV outline
    while true
        xw = wareaXmin + (wareaXmax - wareaXmin)*rand;
        yw = wareaYmin + (wareaYmax - wareaYmin)*rand;
        in = inpolygon(xw, yw, wareaX{1,1}, wareaY{1,1});
        if ~in
            continue;
        end
        % if it is inside the CV area check that it is outside from the
        % islands
        point_in_CV = true;
        for ii = 2:length(wareaX)
            in = inpolygon(xw, yw, wareaX{ii,1}, wareaY{ii,1});
            if in
                point_in_CV = false;
                break;
            end
        end
        if point_in_CV
            break;
        end
    end
    
    % next make sure that it's outside of the river polygons
    id_in = isPoint_inStream([xw yw], STRM, STRM_bbox);
    if ~isempty(id_in)
        continue;
    end
    
    % Make sure that its not too close with the already existing
    % wells
    if cnt_well > 1
        dst = sqrt((xw - wellGen(1:cnt_well-1,1)).^2 + (yw - wellGen(1:cnt_well-1,2)).^2);
        if min(dst) < 200
            continue;
        end
    end
    
    % Finally accept the point with a certaint probability
    
    rw = udens * FwellDens(xw,yw) + (1 - udens) * FwellDens_nrm(xw,yw);
    r= rand;
    if r > rw
        continue;
    end
    
    % if the well has been accepted so far, assign a random pumping rate
    Qw = 10^normrnd(WellQ.Fmean(xw,yw), WellQ.Fstd(xw,yw))/6;
    
    
    elw = Felev(xw,yw);
    botw = Fbot(xw,yw);
    
    % generate a depth and a screen length that fit inside the domain
    count_try = 0;
    add_this_well = false;
    while true
        
        Dmean = WellD.Fa(xw,yw) * log10(Qw) + WellD.Fb(xw,yw);
        D95 = WellD.Fa95(xw,yw) * log10(Qw) + WellD.Fb95(xw,yw);
        Dstd = D95/1.96;
        Dw = 10^normrnd(Dmean, D95);
        
        SLmean = WellS.FCx(xw,yw)*log10(Qw) + WellS.FCy(xw,yw)*log10(Dw) + WellS.FC0(xw,yw);
        SL95 = WellS.FCx95(xw,yw)*log10(Qw) + WellS.FCy95(xw,yw)*log10(Dw) + WellS.FC095(xw,yw);
        SLstd = SL95/1.96;
        SLw = 10^normrnd(SLmean, SLstd);
        
        Bw = elw - Dw;
        Tw = Bw + SLw;
        if Tw - botw > 0 && elw - Tw > 0 && Bw - botw > 0 && elw - Bw > 0
            add_this_well = true;
            break;
        end
        count_try = count_try + 1;
        if count_try > 100
            break;
        end
    end
    
    if add_this_well
        wellGen(cnt_well,:) = [xw yw Tw Bw Qw];
        cnt_well = cnt_well + 1;
        TOT_WELL_Q = TOT_WELL_Q + Qw;
        title(['N:' num2str(cnt_well-1) ', Q:' num2str(28988263.07 - TOT_WELL_Q) '(' num2str(100*TOT_WELL_Q/28988263.07) ')']);
        plot(xw,yw,'.b')
        drawnow
    end
end
wellGen(cnt_well:end,:)=[];
%% print wells into file
load('Generated_wells.mat', 'wellGen25');
wellGen = wellGen25;
% Make sure the pumping is negative
wellGen(:,5) = -wellGen(:,5);
fid = fopen('CVHM_wells.npsat','w');
fprintf(fid, '%d\n', size(wellGen, 1));
fprintf(fid, '%f %f %f %f %f\n', wellGen');
fclose(fid);

%% Analysis PLOTS
CVoutline = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_NPSAT/gis_data/CVHM_Mesh_outline_modif.shp');
%% CVallWells
plot(CVoutline.X,CVoutline.Y,'linewidth',2)
hold on
plot([S.X]', [S.Y]','.')
axis equal
axis off
% print -dpng -r600 CVallWells
%% wells_with_pumping
Q = [Sag.Q]';
Q = Q(~isnan(Q));
XY = [[Sag.X]',[Sag.Y]'];
XY_Q = XY(~isnan([Sag.Q]'),:);
subplot(1,2,1);plot(CVoutline.X,CVoutline.Y)
subplot(1,2,1);hold on
subplot(1,2,1);plot(XY_Q(:,1), XY_Q(:,2),'.')
subplot(1,2,1);title('Wells with Pumping rate')
axis equal
axis off
subplot(1,2,2);histogram(Q)
subplot(1,2,2);xlabel('Well Yield  [gpm]')
%print -dpng -r600 wells_with_pumping
%% Wells_with_D_Q
id = find(~isnan([Sag.Q]') & ~isnan([Sag.depth]'));
D = [Sag.depth]';
Q = [Sag.Q]';
subplot(2,2,[1 3]);hold on
plot(CVoutline.X,CVoutline.Y,'linewidth',2)
subplot(1,2,1);plot(XY(id,1), XY(id,2),'.');
title('Wells with Pumping and Depth');
axis equal
axis off
subplot(2,2,[1 3]);hold off
subplot(2,2,2);histogram(D);
subplot(2,2,2);xlabel('Depth [ft]');
subplot(2,2,2);xlim([0 3000]);
subplot(2,2,4);plot(Q(id),D(id),'.')
subplot(2,2,4);hold on
subplot(2,2,4);xlabel('Pumping rate [gpm]');
subplot(2,2,4);ylabel('Depth [ft]');
D_Q_fit = fit(Q(id),D(id), 'poly1');
ab = coeffvalues(D_Q_fit);
ci = predint(D_Q_fit,[0 :1000:10000]);
ab95 = confint(D_Q_fit);
subplot(2,2,4);plot([0 :1000:10000],[0 :1000:10000].*ab(1) + ab(2),'r','linewidth',2)
subplot(2,2,4);plot([0 :1000:10000]',ci(:,1),'--r')
subplot(2,2,4);plot([0 :1000:10000]',ci(:,2),'--r')
% print -dpng -r600 Wells_with_D_Q
%% Wells_with_D_Q_SL
SL = [Sag.bot]' - [Sag.top]';
id = find(~isnan([Sag.Q]') & ~isnan([Sag.depth]')  & ~isnan(SL));
SL_D_Q_fit = fit([Q(id),D(id)],SL(id), 'poly11');
ab = coeffvalues(SL_D_Q_fit);
[Xg, Yg] = meshgrid(0:2000:10000,0:1000:3000);
cm = SL_D_Q_fit(reshape(Xg,size(Xg,1)*size(Xg,2),1), reshape(Yg,size(Yg,1)*size(Yg,2),1));
ci = predint(SL_D_Q_fit,[reshape(Xg,size(Xg,1)*size(Xg,2),1) reshape(Yg,size(Yg,1)*size(Yg,2),1)]);

%%
clf
subplot(2,2,[1 3]);hold on
plot(CVoutline.X,CVoutline.Y,'linewidth',2)
subplot(1,2,1);plot(XY(id,1), XY(id,2),'.');
title({'Wells with Pumping, Depth', 'and Screenlength'});
axis equal
axis off
subplot(2,2,2);histogram(SL);
subplot(2,2,2);xlabel('Screen length [ft]');
subplot(2,2,2);xlim([0 1000]);
subplot(2,2,4);plot3(Q(id),D(id),SL(id),'.')
subplot(2,2,4);hold on
subplot(2,2,4);xlabel('Pumping rate [gpm]');
subplot(2,2,4);ylabel('Depth [ft]');
subplot(2,2,4);zlabel('Screen length [ft]');
subplot(2,2,4);grid on
subplot(2,2,4);mesh(Xg,Yg,reshape(cm,size(Xg,1),size(Xg,2)),'EdgeColor','r','FaceColor','none');
subplot(2,2,4);mesh(Xg,Yg,reshape(ci(:,1),size(Xg,1),size(Xg,2)),'EdgeColor','b','FaceColor','none','LineStyle','-');
subplot(2,2,4);mesh(Xg,Yg,reshape(ci(:,2),size(Xg,1),size(Xg,2)),'EdgeColor','b','FaceColor','none','LineStyle','-');
%print -dpng -r600 Wells_with_D_Q_SL


