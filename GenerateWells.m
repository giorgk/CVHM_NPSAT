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
WellTable = readtable('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/giorgos.csv');
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
    S(ii,1).Q = WellTable.WellYield(ii);
end
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
WelldataQ = [Sag.Q]';
WelldataDepth = [Sag.Bot]'*0.3048;
WelldataScreenLen = ([Sag.Bot]' - [Sag.Top]')*0.3048;
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
Welldensity = mvksdensity(WelldataXY,XYbas,'bandwidth',5000);

Welldensity(Welldensity==0) = min(Welldensity(Welldensity~=0));
logWellDensity = log10(Welldensity);
logWellDensity_pos = logWellDensity + abs(min(logWellDensity));

maxV = max(logWellDensity_pos);
minV = min(logWellDensity_pos);
Vs   = (logWellDensity_pos - minV) / (maxV - minV);
%%  plot density estimation
clf
trisurf(tri, XYbas(:,1), XYbas(:,2), WellQstd, 'edgecolor', 'none')
hold on
plot(WelldataXY(:,1), WelldataXY(:,2),'.')
view(0,90);
axis equal
alpha(1);
%}
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
%% Compute spatial statistics for each bas point
% For the density simply find out how many wells there are within the threshold

% For the pumping find the weighted mean and weighted standard deviation
% only nonzero positive Q values will be considered
id_incQ = WelldataQ > 0;
logQ = log10(WelldataQ(id_incQ));

id_incSL = WelldataScreenLen > 0;
logSL = log10(WelldataScreenLen(id_incSL));% log screen length
logD = log10(WelldataDepth); % log depth
Thres = 5000;
for jj = 1:size(XYbas,1)
    dst = sqrt((XYbas(jj,1) - WelldataXY(:,1)).^2 + (XYbas(jj,2) - WelldataXY(:,2)).^2);
    w = exp(-(1/Thres)*dst );
    
    wQ = w(id_incQ);
    wSL = w(id_incSL);
    
    WellQmean(jj,1) = 10^sum(logQ.*(wQ/sum(wQ)));
    WellQstd(jj,1) = 10^std(logQ,wQ);
    
    WellSLmean(jj,1) = 10^sum(logSL.*(wSL/sum(wSL)));
    WellSLstd(jj,1) = 10^std(logSL,wSL);
    
    WellDmean(jj,1) = 10^sum(logD.*(w/sum(w)));
    WellDstd(jj,1) = 10^std(logD,w);
    
end
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

%% Prepare streams for faster queries
STRM = readStreams('CVHM_streams.npsat');
% add bounding box info
STRM_bbox = nan(length(STRM),4);
for ii = 1:length(STRM)
    STRM_bbox(ii,:) = [min(STRM(ii,1).poly(:,1)) max(STRM(ii,1).poly(:,1)) ...
                       min(STRM(ii,1).poly(:,2)) max(STRM(ii,1).poly(:,2))];
end

%% ============== Main Algorithm ===============
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

TOT_WELL_Q = 0;
%% main RUN
while TOT_WELL_Q < 28988263.07
    % generate a random point inside the CV bounding box until the point is
    % within CV outline
    while true
        xw = wareaXmin + (wareaXmax - wareaXmin)*rand;
        yw = wareaYmin + (wareaYmax - wareaYmin)*rand;
        in = inpolgon(xw, yw, wareaX{1,1}, wareaY{1,1});
        if ~in
            continue;
        end
        % if it is inside the CV area check that it is outside from the
        % islands
        point_in_CV = true;
        for ii = 2:length(wareaX)
            in = inpolgon(xw, yw, wareaX{ii,1}, wareaY{ii,1});
            if in
                point_in_CV = false;
                break;
            end
        end
        if point_in_CV
            break;
        end
    end
    
end

