%% WELL GENERATION ALGORITHM V2 
%{
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
%}
%%
clear WELLS4CVHM
for ii = 1:length(Farms_poly)
    if ~isempty(find([0] == Farms_poly(ii,1).dwr_sbrgns,1))
        continue;
    end
        
    display(['Farm: ' num2str(Farms_poly(ii,1).dwr_sbrgns)])
    in_wells = inpolygon([CVHMwells.X]',[CVHMwells.Y]', ...
                         Farms_poly(ii,1).X, Farms_poly(ii,1).Y);
                     
    tempW = CVHMwells(in_wells,1);
     
    ipub = find([tempW.type]'==1);
    iag = find([tempW.type]'==0);
    tempWpb = tempW(ipub,1);
    tempWag = tempW(iag,1);
    tempWag_mod = tempWag;
    tempWpb_mod = tempWpb;
                     
    if Farms_poly(ii,1).dwr_sbrgns == 19
        % For the 19th farm use the information from 20 and 21
        in_wells = false(length(CVHMwells),1);
        for kk = 1:length(Farms_poly)
            if ~isempty(find([19 20 21] == Farms_poly(kk,1).dwr_sbrgns,1))
                in_wells(inpolygon([CVHMwells.X]',[CVHMwells.Y]', ...
                         Farms_poly(kk,1).X, Farms_poly(kk,1).Y)) = true;
            end
        end
        tempW = CVHMwells(in_wells,1);
     
        ipub = find([tempW.Itype]'==1);
        iag = find([tempW.Itype]'==0);
        tempWpb = tempW(ipub,1);
        tempWag = tempW(iag,1);
    end
    
   

    % ----------AG wells --------
    %%% PUMPING
    idQ = ~isnan([tempWag.Q]') & [tempWag.Q]' > 0; % isolate the records that have pumping
    display(sum(idQ))
    [f, x]=ecdf([tempWag(idQ,1).Q]'); % create an ECDF based on the data
    % generate random pumping based on the ECDF 
    for jj = 1:length(tempWag_mod)
        if isnan(tempWag_mod(jj,1).Q)
             tempWag_mod(jj,1).Q = interp1q(f,x,rand);
        end
    end
    
    %%% DEPTH
    % isolate the records that have pumping and depth
    idQD = ~isnan([tempWag.Q]') & [tempWag.Q]' > 0 & ...
           ~isnan([tempWag.depth]') & [tempWag.depth]' > 0;
    display(sum(idQD))
    QD_fit = fit(log10([tempWag(idQD,1).Q]'),log10([tempWag(idQD,1).depth]'), 'poly1');
    % generate random depth if depth is missing
    for jj = 1:length(tempWag_mod)
        if isnan(tempWag_mod(jj,1).depth)
            Dm = QD_fit(log10(tempWag_mod(jj,1).Q));
            Ds = predint(QD_fit, log10(tempWag_mod(jj,1).Q),0.68,'observation','off');
            tempWag_mod(jj,1).depth = 10^normrnd(Dm, Dm-Ds(1));
        end
    end
    
    %%% SCREEN LENGTH
    % isolate the records that have pumping and depth and screen length
    idQDS = ~isnan([tempWag.Q]') & [tempWag.Q]' > 0 & ...
            ~isnan([tempWag.depth]') & [tempWag.depth]' > 0 & ...
            ~isnan([tempWag.top]') & ~isnan([tempWag.bot]') & ...
            [tempWag.bot]' - [tempWag.top]' > 0;
    display(sum(idQDS))    
    QDS_fit=fit([log10([tempWag(idQDS,1).Q]'), log10([tempWag(idQDS,1).depth]')], ...
        log10([tempWag(idQDS,1).bot]' - [tempWag(idQDS,1).top]'), 'poly11');
    for jj = 1:length(tempWag_mod)
        if isnan([tempWag_mod(jj,1).top]) || isnan([tempWag_mod(jj,1).bot])
            Sm = QDS_fit(log10([tempWag_mod(jj,1).Q]'), log10([tempWag_mod(jj,1).depth]'));
            Ss = predint(QDS_fit, [log10([tempWag_mod(jj,1).Q]'), log10([tempWag_mod(jj,1).depth]')],0.68,'observation','off');
            tempWag_mod(jj,1).SL = 10^normrnd(Sm, Sm-Ss(1));
        else
            tempWag_mod(jj,1).SL = tempWag_mod(jj,1).bot - tempWag_mod(jj,1).top;
        end
    end
    
    % --------------- Urban wells-------------
    %%% PUMPING
    idQ = ~isnan([tempWpb.Q]') & [tempWpb.Q]' > 0; % isolate the records that have pumping
    display(sum(idQ))
    if sum(idQ) < 10
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
        ipub = find([tempWbas.Itype]'==1);
        tempWpb = tempWbas(ipub,1);
        idQ = ~isnan([tempWpb.Q]') & [tempWpb.Q]' > 0;
        display(sum(idQ))
    end
    
    [f, x]=ecdf([tempWpb(idQ,1).Q]'); % create an ECDF based on the data
    % generate random pumping based on the ECDF 
    for jj = 1:length(tempWpb_mod)
        if isnan(tempWpb_mod(jj,1).Q)
             tempWpb_mod(jj,1).Q = interp1q(f,x,rand);
        end
    end
    
    %%% DEPTH
    % isolate the records that have pumping and depth
    idQD = ~isnan([tempWpb.Q]') & [tempWpb.Q]' > 0 & ...
           ~isnan([tempWpb.depth]') & [tempWpb.depth]' > 0;
    display(sum(idQD))
    QD_fit = fit(log10([tempWpb(idQD,1).Q]'),log10([tempWpb(idQD,1).depth]'), 'poly1');
    % generate random depth if depth is missing
    for jj = 1:length(tempWpb_mod)
        if isnan(tempWpb_mod(jj,1).depth)
            Dm = QD_fit(log10(tempWpb_mod(jj,1).Q));
            Ds = predint(QD_fit, log10(tempWpb_mod(jj,1).Q),0.68,'observation','off');
            tempWpb_mod(jj,1).depth = 10^normrnd(Dm, Dm-Ds(1));
        end
    end
    
    %%% SCREEN LENGTH
    % isolate the records that have pumping and depth and screen length
    idQDS = ~isnan([tempWpb.Q]') & [tempWpb.Q]' > 0 & ...
            ~isnan([tempWpb.depth]') & [tempWpb.depth]' > 0 & ...
            ~isnan([tempWpb.top]') & ~isnan([tempWpb.bot]') & ...
            [tempWpb.bot]' - [tempWpb.top]' > 0;
    display(sum(idQDS))    
    QDS_fit=fit([log10([tempWpb(idQDS,1).Q]'), log10([tempWpb(idQDS,1).depth]')], ...
        log10([tempWpb(idQDS,1).bot]' - [tempWpb(idQDS,1).top]'), 'poly11');
    for jj = 1:length(tempWpb_mod)
        if isnan([tempWpb_mod(jj,1).top]) || isnan([tempWpb_mod(jj,1).bot])
            Sm = QDS_fit(log10([tempWpb_mod(jj,1).Q]'), log10([tempWpb_mod(jj,1).depth]'));
            Ss = predint(QDS_fit, [log10([tempWpb_mod(jj,1).Q]'), log10([tempWpb_mod(jj,1).depth]')],0.68,'observation','off');
            tempWpb_mod(jj,1).SL = 10^normrnd(Sm, Sm-Ss(1));
        else
            tempWpb_mod(jj,1).SL = tempWpb_mod(jj,1).bot - tempWpb_mod(jj,1).top;
        end
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
%}
%% Scale the pumping rates to match the desired pumping
Total_pump_vol = 28988263.07; %m^3/day
init_vol = abs(sum([WELLS4CVHM.UrbanPump]')) + abs(sum([WELLS4CVHM.AgPump]'));
for ii = 1:length(WELLS4CVHM)
    
    urb_perc = abs(WELLS4CVHM(ii,1).UrbanPump) / init_vol;
    ag_perc = abs(WELLS4CVHM(ii,1).AgPump) / init_vol;
    
    urb_act = Total_pump_vol*urb_perc;
    ag_act = Total_pump_vol*ag_perc;
    
    Q_init = [WELLS4CVHM(ii,1).Pubwells.Q]';
    Q_perc = Q_init/sum(Q_init);
    for jj = 1:length(Q_perc)
        WELLS4CVHM(ii,1).Pubwells(jj,1).Q_act = urb_act*Q_perc(jj);
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