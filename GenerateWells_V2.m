%% WELL GENERATION ALGORITHM V2 
%
% In this version we are going to take into account the pumping 
% distribution of CVHM and combined information from the well dataset and
% CVHM model
%% Load Well data set
% bot = depth from land surface to the bottom of the perforated interval [ft]
% top = depth from land surface to the top of the perforated interval [ft]
% bot - top is the length of the screened interval
% Q pumping in gpm
% year of construction
% S = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/CV_only_PubAg.shp');
S = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/CVpbAG_mth.shp');
%% Process well data set
WellTable = readtable('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/giorgos_with_method_ll.csv');
%% Create shapefile from the table
for ii = 1:size(WellTable,1)
    S(ii,1).Geometry = 'Point';
    S(ii,1).X = WellTable.lon(ii,1);
    S(ii,1).Y = WellTable.lat(ii,1);
    S(ii,1).WCRNumber = cell2mat(WellTable.WCRNumber(ii,1));
    S(ii,1).top = str2double(cell2mat(WellTable.top(ii,1)));
    S(ii,1).bot = str2double(cell2mat(WellTable.bot(ii,1)));
    S(ii,1).depth = str2double(cell2mat(WellTable.depth(ii,1)));
    S(ii,1).year = str2double(cell2mat(WellTable.year(ii,1)));
    S(ii,1).Q = str2double(cell2mat(WellTable.Q(ii,1)));
    S(ii,1).type = str2double(cell2mat(WellTable.type(ii,1)));
    if strcmp(cell2mat(WellTable.type(ii,1)), 'agriculture')
        S(ii,1).type = 0;
    elseif strcmp(cell2mat(WellTable.type(ii,1)), 'public')
        S(ii,1).type = 1;
    else
        S(ii,1).type = -9;
    end
    
    if strcmp(cell2mat(WellTable.MethodofDeterminationLL(ii,1)), 'Derived from Address') || ...
           strcmp(cell2mat(WellTable.MethodofDeterminationLL(ii,1)), 'GPS')
        S(ii,1).method = 0;
    else
        S(ii,1).method = 1;
    end
    
end
shapewrite(S, '/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/CVpbAG_mth');
%% Well Age analysis
% discard the wells without year of construction.
% discard prehistoric and other old wells prior to 1900
% discard future wells
id_dlt = isnan([S.year]');
S(id_dlt,:) = [];
id_dlt = [S.year]' < 1900;
S(id_dlt,:) = [];
id_dlt = [S.year]' > 2018;
S(id_dlt,:) = [];
%% plot histogram of age
histogram([S.year]')
xlabel('Year of construction');
%% Find the number of wells for a given threshold
age_thres = 25;
for age_thres = 1:80
    Nwells(age_thres,1) = length(find([S.year]' > 2018 - age_thres));
end
%%
plot(1:80, Nwells/1000,'linewidth',2)
xlabel('Well Age [year]')
ylabel('# of wells [10^3]')
grid on
grid minor
%% spacial distribution for a given age threshold
age_thres = 30;
id_dlt = [S.year]' < 2018 - age_thres;
Sact = S;
Sact(id_dlt,:) = [];
plot([Sact.X]',[Sact.Y]','.')
axis equal
axis off
%% save to local shapefile.
shapewrite(Sact, '/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/CVwells_30y')
%% CVHM WELL PUMPING
% starting with the Bas active shapefile we will add all the well info on
% each cell.
Bas = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_NPSAT/gis_data/BAS_active.shp');
%%
for ii = [1 2 3 45 6 7 8 9 10]
    Bas = rmfield(Bas,['lay' num2str(ii) '_act']);
end
for ii = 1:10
    Bas = rmfield(Bas,['strt_hd_' num2str(ii,'%02.0f')]);
end
Bas = rmfield(Bas, 'grid_active');

CVHM_wells = Bas;
%% load Well Budgets
cbcf_path = '/home/giorgk/Documents/UCDAVIS/CVHM_DATA/ClaudiaRun/';
m=4;y=1961;
t2=0;
for ii=1:510
    if m==13
        y
        m = 1;
        y = y+1;
    end
    cnst=floor((ii-1)/10);
    if t2<cnst*10+1
        clear OUT
        load([cbcf_path 'cbcf_' num2str(cnst*10+1) '_' num2str(10*cnst+10) '.mat']);
        t2=cnst*10+1;
    end
    Welldates(ii,1)=y;
    Welldates(ii,2)=m;
    FARMWELLS{ii,1} = OUT{ii,1}{10,2};
    MNWELLS{ii,1} = OUT{ii,1}{9,2};
    m = m + 1;
end
%% OR
load('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/Well_Budget.mat');
%% 
yr = 1977; % year start
mnt = 10; % month start
ye = 2003; % year end
me = 2; % month end
Ndays = [];
while 1
    Ndays = [Ndays; eomday(yr,mnt)];
    mnt = mnt + 1;
    if mnt > 12
        mnt = 1;
        yr = yr + 1;
    end
    if yr == ye && mnt > me
        break;
    end
end
%
%% Average the stresses for the simulation period
% get istart and iend from GWaterBudget_calc
iid = istart:iend;
MNW = zeros(441, 98);
FRW = zeros(441, 98);
for ii = 1:length(Ndays)
    for jj = 1:size(FARMWELLS{iid(ii),1},1)
        FRW(FARMWELLS{iid(ii),1}(jj,2), FARMWELLS{iid(ii),1}(jj,3)) = ...
            FRW(FARMWELLS{iid(ii),1}(jj,2), FARMWELLS{iid(ii),1}(jj,3)) + ...
            FARMWELLS{iid(ii),1}(jj,4)*Ndays(ii,1);
    end
    for jj = 1:size(MNWELLS{iid(ii),1},1)
        MNW(MNWELLS{iid(ii),1}(jj,2), MNWELLS{iid(ii),1}(jj,3)) = ...
            MNW(MNWELLS{iid(ii),1}(jj,2), MNWELLS{iid(ii),1}(jj,3)) + ...
            MNWELLS{iid(ii),1}(jj,4)*Ndays(ii,1);
    end
end
FRW = FRW./sum(Ndays);
MNW = MNW./sum(Ndays);
%% add the well rates to CVHM_wells shapefiles
for ii = 1:size(CVHM_wells,1)
    CVHM_wells(ii,1).farm = FRW(CVHM_wells(ii,1).ROW, CVHM_wells(ii,1).COLUMN_);
    CVHM_wells(ii,1).mnw = MNW(CVHM_wells(ii,1).ROW, CVHM_wells(ii,1).COLUMN_);
end
shapewrite(CVHM_wells, 'gis_data/CVHM_wellBud');
%
%% Find the volume of water that is used for Urban and Ag for each farm
CVHMwells = shaperead('gis_data/CVHM_wellBud');
Farms = shaperead('gis_data/FMP.shp');
Farms_poly = shaperead('gis_data/FARMS_poly.shp');
for ii = 1:length(Farms_poly)
    Farms_poly(ii,1).UrbanPump = 0;
    Farms_poly(ii,1).AgPump = 0;
end
ROWfrm = [Farms.ROW]';
COLfrm = [Farms.COLUMN_]';
%
for ii = 1:size(CVHMwells,1)
    r = CVHMwells(ii,1).ROW;
    c = CVHMwells(ii,1).COLUMN_;
    ir = find(ROWfrm == r & COLfrm == c);
    ifarm = find([Farms_poly.dwr_sbrgns]' == Farms(ir,1).dwr_sbrgns);
    if Farms(ir,1).LU_2000 == 2
        Farms_poly(ifarm,1).UrbanPump = Farms_poly(ifarm,1).UrbanPump + ...
            CVHMwells(ii,1).farm + CVHMwells(ii,1).mnw;
    else
        Farms_poly(ifarm,1).AgPump = Farms_poly(ifarm,1).AgPump + ...
            CVHMwells(ii,1).farm + CVHMwells(ii,1).mnw;
    end
end
shapewrite(Farms_poly, 'gis_data/FARMS_polyPump.shp');
%
%% Generate Pumping wells
Farm_basin = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_NPSAT/gis_data/FARMS_poly.shp');
Farms_poly = shaperead('gis_data/FARMS_polyPump.shp');
CVHMwells = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/CVwells_30y_proj.shp');
cvhmDomain = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_NPSAT/gis_data/BAS_active_outline.shp');
%
%%
clear WELLS4CVHM
for ii = 1:length(Farms_poly)
    if ~isempty(find([0] == Farms_poly(ii,1).dwr_sbrgns,1))
        continue;
    end
        
    display(['Farm: ' num2str(Farms_poly(ii,1).dwr_sbrgns)])
    in_wells = inpolygon([CVHMwells.X]',[CVHMwells.Y]', ...
                         Farms_poly(ii,1).X, Farms_poly(ii,1).Y) & ... 
               inpolygon([CVHMwells.X]',[CVHMwells.Y]', ...
                         cvhmDomain.X, cvhmDomain.Y);
                     
    tempW = CVHMwells(in_wells,1);
     
    ipub = find([tempW.type]'==1);
    iag = find([tempW.type]'==0);
    tempWpb = tempW(ipub,1);
    tempWag = tempW(iag,1);
    tempWag_mod = tempWag;
    tempWpb_mod = tempWpb;
    if length(tempWpb_mod) < 15
        % use all public records that belong to basin
        ibas = Farm_basin(ii,1).Basins;
        in_basin = false(length(CVHMwells),1);
        for kk = 1:length(Farm_basin)
            if Farm_basin(kk,1).Basins ~= ibas
                continue;
            end
            in_basin(inpolygon([CVHMwells.X]',[CVHMwells.Y]', ...
                         Farm_basin(kk,1).X, Farm_basin(kk,1).Y)) = true;
        end
        tempWbas = CVHMwells(in_basin,1);
        ipub = find([tempWbas.type]'==1);
        tempWpb = tempWbas(ipub,1);
    end

    % ----------AG wells --------
    %%% Pumping ECDF
    idQ = ~isnan([tempWag.Q]') & [tempWag.Q]' > 0; % isolate the records that have pumping
    display(sum(idQ))
    logQ = log10([tempWag(idQ,1).Q]');
    [fQ, xQ]=ecdf(logQ); % create an ECDF based on the data
    
    %%% Depth ECDF
    idD = ~isnan([tempWag.depth]') & [tempWag.depth]' > 0;
    display(sum(idD))
    logD = log10([tempWag(idD,1).depth]');
    [fD, xD]=ecdf(logD);
    
    %%% Screen length ECDF
    tempSL = [tempWag.bot]' - [tempWag.top]';
    idS = ~isnan(tempSL) & tempSL < 2000  & tempSL > 0;
    display(sum(idS))
    logS = log10([tempWag(idS,1).bot]' - [tempWag(idS,1).top]');
    [fS, xS]=ecdf(logS);
    
    clear Stats
    Stats.fQ = fQ; Stats.lowQ = 0.15;
    Stats.xQ = xQ; Stats.uppQ = 0.85;
    Stats.fD = fD; Stats.lowD = 0;
    Stats.xD = xD; Stats.uppD = 1;
    Stats.fS = fS; Stats.lowS = 0;
    Stats.xS = xS; Stats.uppS = 1;
    
    for jj = 1:length(tempWag_mod)
        [ tempWag_mod(jj,1).Q, ...
          tempWag_mod(jj,1).depth, ...
          tempWag_mod(jj,1).SL ] = AssignQDS(Stats, ...
                                             tempWag_mod(jj,1).Q, ...
                                             tempWag_mod(jj,1).depth, ...
                                             tempWag_mod(jj,1).bot - tempWag_mod(jj,1).top);
    end

    
    % --------------- Urban wells-------------
    %%% Pumping ECDF
    idQ = ~isnan([tempWpb.Q]') & [tempWpb.Q]' > 0; % isolate the records that have pumping
    display(sum(idQ))
    if sum(idQ) < 15
        % use all public records that belong to basin
        ibas = Farm_basin(ii,1).Basins;
        in_basin = false(length(CVHMwells),1);
        for kk = 1:length(Farm_basin)
            if Farm_basin(kk,1).Basins ~= ibas
                continue;
            end
            in_basin(inpolygon([CVHMwells.X]',[CVHMwells.Y]', ...
                         Farm_basin(kk,1).X, Farm_basin(kk,1).Y)) = true;
        end
        tempWbas = CVHMwells(in_basin,1);
        ipub = find([tempWbas.type]'==1);
        tempWpbQ = tempWbas(ipub,1);
        idQ = ~isnan([tempWpbQ.Q]') & [tempWpbQ.Q]' > 0;
        display(sum(idQ))
        logQ = log10([tempWpbQ(idQ,1).Q]');
        [fQ, xQ]=ecdf(logQ);
    else
        logQ = log10([tempWpb(idQ,1).Q]');
        [fQ, xQ]=ecdf(logQ);
    end
    
    
    
    %%% Depth ECDF
    idD = ~isnan([tempWpb.depth]') & [tempWpb.depth]' > 0;
    display(sum(idD))
    logD = log10([tempWpb(idD,1).depth]');
    [fD, xD]=ecdf(logD);
    
    %%% Screen length ECDF
    tempSL = [tempWpb.bot]' - [tempWpb.top]';
    idS = ~isnan(tempSL) & tempSL < 2000  & tempSL > 0;
    display(sum(idS))
    logS = log10([tempWpb(idS,1).bot]' - [tempWpb(idS,1).top]');
    [fS, xS]=ecdf(logS);
    
    clear Stats
    Stats.fQ = fQ; Stats.lowQ = 0.15;
    Stats.xQ = xQ; Stats.uppQ = 0.85;
    Stats.fD = fD; Stats.lowD = 0;
    Stats.xD = xD; Stats.uppD = 1;
    Stats.fS = fS; Stats.lowS = 0;
    Stats.xS = xS; Stats.uppS = 1;
    
    for jj = 1:length(tempWpb_mod)
        [ tempWpb_mod(jj,1).Q, ...
          tempWpb_mod(jj,1).depth, ...
          tempWpb_mod(jj,1).SL ] = AssignQDS(Stats, ...
                                             tempWpb_mod(jj,1).Q, ...
                                             tempWpb_mod(jj,1).depth, ...
                                             tempWpb_mod(jj,1).bot - tempWpb_mod(jj,1).top);
    end
    
    
    
    
    WELLS4CVHM(ii,1).farm_id = Farms_poly(ii,1).dwr_sbrgns;
    WELLS4CVHM(ii,1).UrbanPump = Farms_poly(ii,1).UrbanPump;
    WELLS4CVHM(ii,1).AgPump = Farms_poly(ii,1).AgPump;
    WELLS4CVHM(ii,1).Pubwells = tempWpb_mod;
    WELLS4CVHM(ii,1).Agwells = tempWag_mod;
    WELLS4CVHM(ii,1).Npb = length(tempWpb_mod);
    WELLS4CVHM(ii,1).Nag = length(tempWag_mod);
end
WELLS4CVHM(17,:) = [];
%
%% Scale the pumping rates to match the desired pumping
Total_pump_vol = 28988263.07; %m^3/day
Min_pumping = 400;
init_vol = abs(sum([WELLS4CVHM.UrbanPump]')) + abs(sum([WELLS4CVHM.AgPump]'));
for ii = 1:length(WELLS4CVHM)
    
    urb_perc = abs(WELLS4CVHM(ii,1).UrbanPump) / init_vol;
    ag_perc = abs(WELLS4CVHM(ii,1).AgPump) / init_vol;
    
    urb_act = Total_pump_vol*urb_perc;
    ag_act = Total_pump_vol*ag_perc;
    
    %%% Public wells
    Q_init = [WELLS4CVHM(ii,1).Pubwells.Q]';
    Q_perc = Q_init/sum(Q_init);
    Q_temp = urb_act*Q_perc;
    c = min(Q_temp);
    if c < Min_pumping
        cnt = 1;
        while 1
            Q_init = [WELLS4CVHM(ii,1).Pubwells.Q]';
            [c, d] = sort(Q_init);
            Q_init(d(1:cnt,1),:) = [];
            Q_perc = Q_init/sum(Q_init);
            Q_temp = urb_act*Q_perc;
            c = min(Q_temp);
            if c > Min_pumping
                break;
            else
                cnt = cnt + 1;
            end
        end
        WELLS4CVHM(ii,1).Pubwells(d(1:cnt,1),:) = [];
    end
   
    Q_init = [WELLS4CVHM(ii,1).Pubwells.Q]';
    Q_perc = Q_init/sum(Q_init);
    for jj = 1:length(Q_perc)
        WELLS4CVHM(ii,1).Pubwells(jj,1).Q_act = urb_act*Q_perc(jj);
    end
    
    %%% Ag wells
    Q_init = [WELLS4CVHM(ii,1).Agwells.Q]';
    Q_perc = Q_init/sum(Q_init);
    Q_temp = ag_act*Q_perc;
    c = min(Q_temp);
    if c < Min_pumping
        cnt = 1;
        while 1
            Q_init = [WELLS4CVHM(ii,1).Agwells.Q]';
            [c, d] = sort(Q_init);
            Q_init(d(1:cnt,1),:) = [];
            Q_perc = Q_init/sum(Q_init);
            Q_temp = ag_act*Q_perc;
            c = min(Q_temp);
            if c > Min_pumping
                break;
            else
                cnt = cnt + 1;
            end
        end
        WELLS4CVHM(ii,1).Agwells(d(1:cnt,1),:) = [];
    end
    
    Q_init = [WELLS4CVHM(ii,1).Agwells.Q]';
    Q_perc = Q_init/sum(Q_init);
    for jj = 1:length(Q_perc)
        WELLS4CVHM(ii,1).Agwells(jj,1).Q_act = ag_act*Q_perc(jj);
    end
end
%% Check total volume
Q_ag = 0;
Q_pb = 0;
for ii = 1:length(WELLS4CVHM)
    Q_ag = Q_ag + sum([WELLS4CVHM(ii,1).Agwells.Q_act]);
    Q_pb = Q_pb + sum([WELLS4CVHM(ii,1).Pubwells.Q_act]);
end
%% Put it all in one variable table
WW = nan(20000,7);
cnt = 1;
for ii = 1:size(WELLS4CVHM,1)
    for jj = 1:length(WELLS4CVHM(ii,1).Agwells)
        WW(cnt,:) = [WELLS4CVHM(ii,1).Agwells(jj,1).X ...
                     WELLS4CVHM(ii,1).Agwells(jj,1).Y ...
                     WELLS4CVHM(ii,1).Agwells(jj,1).depth ...
                     WELLS4CVHM(ii,1).Agwells(jj,1).SL ...
                     WELLS4CVHM(ii,1).Agwells(jj,1).Q_act ...
                     ii jj];
        cnt = cnt + 1;
    end
    
    for jj = 1:length(WELLS4CVHM(ii,1).Pubwells)
        WW(cnt,:) = [WELLS4CVHM(ii,1).Pubwells(jj,1).X ...
                     WELLS4CVHM(ii,1).Pubwells(jj,1).Y ...
                     WELLS4CVHM(ii,1).Pubwells(jj,1).depth ...
                     WELLS4CVHM(ii,1).Pubwells(jj,1).SL ...
                     WELLS4CVHM(ii,1).Pubwells(jj,1).Q_act ...
                     ii jj];
        cnt = cnt + 1;
    end
end
WW(cnt:end,:) = [];
WW_bck = WW; % make a backup table
%% Randomize well locations
WW = WW_bck;
for ii = 1:size(WW,1)
    ii
    c = sqrt((WW(:,1) - WW(ii,1)).^2 + (WW(:,2) - WW(ii,2)).^2);
    id = find(c < 400);
    if length(id) > 1
        xm = WW(ii,1);
        ym = WW(ii,2);
        for jj = 1:length(id)
            cnt = 0;
            thres = 400;
            ifarm = WELLS4CVHM(WW(id(jj),6),1).farm_id;% This is the id of the farm
            ifarm = find([Farms_poly.dwr_sbrgns]' == ifarm); % This is the index of the farm in Farm_poly
            
            while 1
                while 1 % make sure the generated well is inside the farm
                    xr = xm - 600 + (xm + 600 - (xm - 600))*rand;
                    yr = ym - 600 + (ym + 600 - (ym - 600))*rand;
                    in = inpolygon(xr, yr, Farms_poly(ifarm,1).X, Farms_poly(ifarm,1).Y) & ...
                         inpolygon(xr, yr, cvhmDomain.X, cvhmDomain.Y);
                    if in
                        break;
                    end
                end
                if jj == 1
                    WW(id(jj),1) = xr;
                    WW(id(jj),2) = yr;
                    break;
                else
                    dst = sqrt((WW(id(1:jj-1),1) - xr).^2 + (WW(id(1:jj-1),2)- yr).^2);
                    if min(dst) > thres
                        WW(id(jj),1) = xr;
                        WW(id(jj),2) = yr;
                        break;
                    end
                end
                cnt = cnt + 1;
                if cnt > 10
                    thres = thres - 10;
                    thres
                    cnt = 0;
                end
            end
        end
    end
end
%% check minimum distances
for ii = 1:size(WW,1)
    tii = setdiff([1:size(WW,1)]',ii);
    mindst(ii,1) = min(sqrt((WW(tii,1) - WW(ii,1)).^2 + (WW(tii,2) - WW(ii,2)).^2));
end
%% relocate wells that are closer than 400 m
[c d] = sort(mindst);
for ii = 1:length(c)
    ii
    if c(ii) > 400
        break;
    end
    tii = setdiff([1:size(WW,1)]',d(ii));
    xc = WW(d(ii),1);
    yc = WW(d(ii),2);
    ifarm = WELLS4CVHM(WW(d(ii),6),1).farm_id;% This is the id of the farm
    ifarm = find([Farms_poly.dwr_sbrgns]' == ifarm); % This is the index of the farm in Farm_poly
    radius = 10;
    cnt = 0;
    while 1
        % pick a random direction and a random radius
        while 1
            r_rad = radius*rand;
            tr = 2*pi*rand;
            xr = xc + r_rad*cos(tr);
            yr = yc + r_rad*sin(tr);
            in = inpolygon(xr, yr, Farms_poly(ifarm,1).X, Farms_poly(ifarm,1).Y) & ...
                 inpolygon(xr, yr, cvhmDomain.X, cvhmDomain.Y);
            if in
                break;
            end
        end
        
        dst = min(sqrt((WW(tii,1) - xr).^2 + (WW(tii,2) - yr).^2));
        if dst > 400
            WW(d(ii),1) = xr;
            WW(d(ii),2) = yr;
            break;
        end
        cnt = cnt + 1;
        if cnt > 100
            radius = radius + 10;
            cnt = 0;
        end
    end
end
%% Correct those that their Screen Lengths are greater or equal to depth
id = find(WW(:,3) <= WW(:,4));
WW(id,3) = WW(id,4) + 30;
%% Assign top and bottom according to the artificial elevation and CVHM bottom
topelev = read_Scattered('CVHM_top_elev.npsat',2);
Ftop = scatteredInterpolant(topelev.p(:,1),topelev.p(:,2),topelev.v);

% abd the bottom of the aquifer
botelev = read_Scattered('CVHM_Bot_elev.npsat',2);
Fbot = scatteredInterpolant(botelev.p(:,1),botelev.p(:,2),botelev.v);
%}
%%
cnvrs = 0.3048;
for ii = 1:size(WW,1)
    Wt = Ftop(WW(ii,1), WW(ii,2));
    Bs = Fbot(WW(ii,1), WW(ii,2));
    aquif_height = Wt - Bs;
    if aquif_height < WW(ii,3)*cnvrs
        if aquif_height > 30
            botw = Bs + 10;
        else
            wtf = true;
        end
    else
        botw = Wt - WW(ii,3)*cnvrs;
    end
    
    if botw + WW(ii,4)*cnvrs > Wt
        if aquif_height > 30
            topw = Wt - 10;
        else
            wtf = true;
        end
    else
        topw = botw+WW(ii,4)*cnvrs;
    end
    WW(ii,8) = topw;
    WW(ii,9) = botw;
    WW(ii,10) = Wt;
    WW(ii,11) = Bs;
end

%% Correct the wells with 
id = find(WW(:,9) < WW(:,11))
%% Print to file

fid = fopen('CVHM_wells.npsat','w');
fprintf(fid, '%d\n', size(WW, 1));
fprintf(fid, '%f %f %f %f -%f\n', [WW(:,1) WW(:,2) WW(:,8) WW(:,9) WW(:,5)]');
fclose(fid);
