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
%}
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
    
