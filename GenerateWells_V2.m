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
S = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/CV_only_PubAg.shp');
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
Farms_poly = shaperead('gis_data/FARMS_polyPump.shp');
CVHMwells = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/CVwells_30y.shp');
for ii = 1:length(CVHMwells)
    if strcmp(CVHMwells(ii,1).type,'public')
        CVHMwells(ii,1).Itype = 1;
    else
        CVHMwells(ii,1).Itype = 0;
    end
end
%}
%%
for ii = 1:length(Farms_poly)
    in_wells = inpolygon([CVHMwells.X]',[CVHMwells.Y]', ...
                         Farms_poly(ii,1).X, Farms_poly(ii,1).Y);
                     
    tempW = CVHMwells(in_wells,1);
     
    ipub = find([tempW.Itype]'==1);
    iag = find([tempW.Itype]'==0);
    tempWpb = tempW(ipub,1);
    tempWag = tempW(iag,1);
    tempWag_mod = tempWag;
    tempWpb_mod = tempWpb;

    % ----------AG wells --------
    %%% PUMPING
    idQ = ~isnan([tempWag.Q]); % isolate the records that have pumping
    [f, x]=ecdf([tempWag(idQ,1).Q]'); % create an ECDF based on the data
    % generate random pumping based on the ECDF 
    for jj = 1:length(tempWag)
        if isnan(tempWag(jj,1).Q)
             tempWag_mod(jj,1).Q = interp1q(f,x,rand);
        end
    end
    
    %%% DEPTH
    % isolate the records that have pumping and depth
    idQD = ~isnan([tempWag.Q]) & ~isnan([tempWag.depth]);
    QD_fit = fit(log10([tempWag(idQD,1).Q]'),log10([tempWag(idQD,1).depth]'), 'poly1');
    % generate random depth if depth is missing
    for jj = 1:length(tempWag)
        if isnan(tempWag(jj,1).depth)
            Dm = QD_fit(log10(tempWag_mod(jj,1).Q));
            Ds = predint(QD_fit, log10(tempWag_mod(jj,1).Q),0.68,'observation','off');
            tempWag_mod(jj,1).depth = 10^normrnd(Dm, Dm-Ds(1));
        end
    end
    
    %%% SCREEN LENGTH
    % isolate the records that have pumping and depth and screen length
    idQDS = ~isnan([tempWag.Q]) & ~isnan([tempWag.depth]) & ...
            ~isnan([tempWag.top]) & ~isnan([tempWag.bot]);
        
    QDS_fit=fit([log10([tempWag(idQDS,1).Q]'), log10([tempWag(idQDS,1).depth]')], ...
        log10([tempWag(idQDS,1).bot]' - [tempWag(idQDS,1).top]'), 'poly11');
    for jj = 1:length(tempWag)
        if isnan([tempWag(jj,1).top]) || isnan([tempWag(jj,1).bot])
            Sm = QDS_fit(log10([tempWag_mod(jj,1).Q]'), log10([tempWag_mod(jj,1).depth]'));
            Ss = predint(QDS_fit, [log10([tempWag_mod(jj,1).Q]'), log10([tempWag_mod(jj,1).depth]')],0.68,'observation','off');
            tempWag_mod(jj,1).SL = 10^normrnd(Sm, Sm-Ss(1));
        else
            tempWag_mod(jj,1).SL = tempWag(jj,1).bot - tempWag(jj,1).top;
        end
    end
    
    % --------------- Urban wells-------------
    %%% PUMPING
    idQ = ~isnan([tempWpb.Q]); % isolate the records that have pumping
    [f, x]=ecdf([tempWpb(idQ,1).Q]'); % create an ECDF based on the data
    % generate random pumping based on the ECDF 
    for jj = 1:length(tempWpb)
        if isnan(tempWpb(jj,1).Q)
             tempWpb_mod(jj,1).Q = interp1q(f,x,rand);
        end
    end
    
    %%% DEPTH
    % isolate the records that have pumping and depth
    idQD = ~isnan([tempWpb.Q]) & ~isnan([tempWpb.depth]);
    QD_fit = fit(log10([tempWpb(idQD,1).Q]'),log10([tempWpb(idQD,1).depth]'), 'poly1');
    % generate random depth if depth is missing
    for jj = 1:length(tempWpb)
        if isnan(tempWpb(jj,1).depth)
            Dm = QD_fit(log10(tempWpb(jj,1).Q));
            Ds = predint(QD_fit, log10(tempWpb_mod(jj,1).Q),0.68,'observation','off');
            tempWpb_mod(jj,1).depth = 10^normrnd(Dm, Dm-Ds(1));
        end
    end
    
    %%% SCREEN LENGTH
    % isolate the records that have pumping and depth and screen length
    idQDS = ~isnan([tempWpb.Q]) & ~isnan([tempWpb.depth]) & ...
            ~isnan([tempWpb.top]) & ~isnan([tempWpb.bot]);
        
    QDS_fit=fit([log10([tempWpb(idQDS,1).Q]'), log10([tempWpb(idQDS,1).depth]')], ...
        log10([tempWpb(idQDS,1).bot]' - [tempWpb(idQDS,1).top]'), 'poly11');
    for jj = 1:length(tempWpb)
        if isnan([tempWpb(jj,1).top]) || isnan([tempWpb(jj,1).bot])
            Sm = QDS_fit(log10([tempWpb_mod(jj,1).Q]'), log10([tempWpb_mod(jj,1).depth]'));
            Ss = predint(QDS_fit, [log10([tempWpb_mod(jj,1).Q]'), log10([tempWpb_mod(jj,1).depth]')],0.68,'observation','off');
            tempWpb_mod(jj,1).SL = 10^normrnd(Sm, Sm-Ss(1));
        else
            tempWpb_mod(jj,1).SL = tempWpb(jj,1).bot - tempWpb(jj,1).top;
        end
    end
    
end
