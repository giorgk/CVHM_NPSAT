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
%% Add the CVHM Surface elevation 
[ ELEV, info ] = readArcGisASCIIfile( 'gis_data/PP1766_DIS_Layer1_Top_ASCII.txt' );
Xgrid = info.xllcorner+info.cellsize/2:...
                info.cellsize:...
                info.xllcorner+info.cellsize/2+info.cellsize*(info.ncols-1);
    
Ygrid = info.yllcorner+info.cellsize/2:...
        info.cellsize:...
        info.yllcorner+info.cellsize/2+info.cellsize*(info.nrows-1);

Felev = griddedInterpolant({Xgrid, Ygrid}, ELEV');
%% Interpolate Elevation to well nodes
Wellelev = Felev(WelldataXY(:,1), WelldataXY(:,2))*0.3048; %convert to mm
%% 
WelldataXY = [[S.X]' [S.Y]'];
WelldataQ = [S.Q]';
WelldataTop = [S.Top]'*0.3048;
WelldataDepth = Wellelev - WelldataTop;
WelldataScreenLen = (WelldataTop - [S.Bot]')*0.3048;
%% Bring the basic shapefile
% Run the first 3 sections of the PrepareGeometryData.m script
% make a list of the barycenters of the bas shapefile
XYbas = nan(length(bas),2);
for ii = 1:length(bas)
    XYbas(ii,:) = [mean(bas(ii,1).X(1:4)) mean(bas(ii,1).Y(1:4))];
end
%}
%% Calculate the KDE for the density of wells
f = mvksdensity(WelldataXY,XYbas,'bandwidth',50);
%% Compute spatial statistics for each bas point
% For the density simply find out how many wells there are within the threshold

% For the pumping find the weighted mean and weighted standard deviation

Thres = 30000;
for jj = 1:size(XYbas,1)
    dst = sqrt((WelldataXY(:,1) - XYbas(jj,1)).^2 + (WelldataXY(:,2) - XYbas(jj,2)).^2);
    id = find(dst < Thres);
    Welldensity(jj,1) = length(id);
    
    w = Thres - dst(id,1);
    w = w/sum(w);
    
    WellQmean(jj,1) = sum(WelldataQ(id,1).*w);
    WellQstd(jj,1) = std(WelldataQ(id,1),w);
    
    WellSLmean(jj,1) = sum(WelldataScreenLen(id,1).*w);
    WellSLstd(jj,1) = std(WelldataScreenLen(id,1),w);
    
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
%% plot
trisurf(tri, XYbas(:,1),XYbas(:,2),Welldensity/max(Welldensity),'edgecolor','none')
%% Prepare streams for faster queries
STRM = readStreams('CVHM_streams.npsat');
% add bounding box info
STRM_bbox = nan(length(STRM),4);
for ii = 1:length(STRM)
    STRM_bbox(ii,:) = [min(STRM(ii,1).poly(:,1)) max(STRM(ii,1).poly(:,1)) ...
                       min(STRM(ii,1).poly(:,2)) max(STRM(ii,1).poly(:,2))];
end
    
%% ============== Main Algorithm ===============

