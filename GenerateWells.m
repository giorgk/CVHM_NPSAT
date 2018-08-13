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
WelldataXY = [[S.X]' [S.Y]'];
%% Bring the basic shapefile
% Run the first 3 sections of the PrepareGeometryData.m script
% make a list of the barycenters of the bas shapefile
XYbas = nan(length(bas),2);
for ii = 1:length(bas)
    XYbas(ii,:) = [mean(bas(ii,1).X(1:4)) mean(bas(ii,1).Y(1:4))];
end
%% Calculate a density indicator for each bas point
Thres = 20000;
for jj = 1:size(XYbas,1)
    Dens(jj,1) = sum(sqrt((WelldataXY(:,1) - XYbas(jj,1)).^2 + (WelldataXY(:,2) - XYbas(jj,2)).^2) < Thres);
end