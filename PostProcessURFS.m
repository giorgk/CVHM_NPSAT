%% First put them all in one variable
%
urf_path = 'output/URFS/';
%%
tempURFS = [];
load([urf_path 'WellURF_0_0_50.mat']);
tempURFS = [tempURFS; WellURF];
load([urf_path 'WellURF_0_50_99.mat'])
tempURFS = [tempURFS; WellURF];
load([urf_path 'WellURF_0_100_127.mat'])
tempURFS = [tempURFS; WellURF];
load([urf_path 'WellURF_1_0_127.mat'])
tempURFS = [tempURFS; WellURF];
load([urf_path 'WellURF_2_5_127.mat'])
tempURFS = [tempURFS; WellURF];
clear WellURF
%%
clear URFS
Nurfs = length(tempURFS);
URFS(Nurfs,1).Eid = nan;
URFS(Nurfs,1).Sid = nan;
URFS(Nurfs,1).X = nan;
URFS(Nurfs,1).Y = nan;
URFS(Nurfs,1).Vland = nan;
URFS(Nurfs,1).Vcds = nan;
URFS(Nurfs,1).riv = 0;
%% raw URFS
for ii = 1:Nurfs
    URFS(ii,1).X = tempURFS(ii,1).p_lnd(1);
    URFS(ii,1).Y = tempURFS(ii,1).p_lnd(2);
    URFS(ii,1).Eid = tempURFS(ii,1).Eid;
    URFS(ii,1).Sid = tempURFS(ii,1).Sid;
    URFS(ii,1).Vland = tempURFS(ii,1).v_lnd(1);
    URFS(ii,1).Vcds = tempURFS(ii,1).v_cds(1);
    URFS(ii,1).riv = 0;
    URFS(ii,1).urf = tempURFS(ii,1).URF;
end
clear tempURFS
%% with fitting
ft = fittype( '(1/(x*b*sqrt(2*pi)))*exp((-(log(x)-a)^2)/(2*b^2))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.499112701571756 0.336174051321482];
X = 1:200;
for ii = 1:Nurfs
    if rem(ii,1000) == 0
        display(ii)
    end
    URFS(ii,1).X = tempURFS(ii,1).p_lnd(1);
    URFS(ii,1).Y = tempURFS(ii,1).p_lnd(2);
    URFS(ii,1).Eid = tempURFS(ii,1).Eid;
    URFS(ii,1).Sid = tempURFS(ii,1).Sid;
    URFS(ii,1).Vland = tempURFS(ii,1).v_lnd(1);
    URFS(ii,1).Vcds = tempURFS(ii,1).v_cds(1);
    URFS(ii,1).riv = 0;
    xx = tempURFS(ii,1).URF;
    
    % Simplify URF
%     xx(xx<0) = 0;
%     if max(xx) < 1e-10
%         URFS(ii,1).urf = [];
%     else
%         tol = max(xx)*0.025;
%         xs = unique(simplifyURF([1:200 ;xx]', max(xx)*0.025),'rows');
%         URFS(ii,1).urf.x = xs(:,1);
%         URFS(ii,1).urf.y = xs(:,2);
%     end
    
    if max(xx) < 1e-5
        URFS(ii,1).urf.mu = -100;
        URFS(ii,1).urf.std = 1;
    else
        % or fit a distribution
        [xData, yData] = prepareCurveData( X, xx );
        [fitresult, gof] = fit( xData, yData, ft, opts );
        ab = coeffvalues(fitresult);
        URFS(ii,1).urf.mu = ab(1);
        URFS(ii,1).urf.std = ab(2);
    end
    
end
clear tempURFS
%% Find out which streamlines end up on rivers
Xp = [URFS.X]';
Yp = [URFS.Y]';
streams = readStreams('CVHM_streams.npsat');
%}
%% 
for ii = 1:length(streams)
    if rem(ii,100) == 0
        display(ii)
    end
    id = find(inpolygon(Xp, Yp, streams(ii,1).poly(:,1),streams(ii,1).poly(:,2)));
    
    for jj = 1:size(streams(ii,1).poly,1)-1
        L = [streams(ii,1).poly(jj,:) streams(ii,1).poly(jj+1,:)];
        dst = Dist_Point_LineSegment(Xp, Yp, L);
        idd = find(dst < 30);
    end
    id = unique([id;idd]);
    
    for jj = 1:length(id)
        URFS(id(jj),1).riv = 1;
        URFS(id(jj),1).urf = [];
        %URFS(ii,1).urf.mu = -100;
        %URFS(ii,1).urf.std = 1;
    end
end
%% Create shapefile with the end points of streamlines
S(length(URFS),1).Geometry = 'Point';
for ii = 1:length(URFS)
    S(ii,1).Geometry = 'Point';
    S(ii,1).X = URFS(ii,1).X;
    S(ii,1).Y = URFS(ii,1).Y;
    S(ii,1).Eid = double(URFS(ii,1).Eid);
    S(ii,1).Sid = double(URFS(ii,1).Sid);
    
    S(ii,1).Vland = URFS(ii,1).Vland;
    S(ii,1).Vcds = URFS(ii,1).Vcds;
    S(ii,1).riv = URFS(ii,1).riv;
    S(ii,1).mu = URFS(ii,1).urf.mu;
    S(ii,1).std = URFS(ii,1).urf.std;
end
%%
% shapewrite(S, [urf_path 'CVHM_strmln_Endpnt']);
shapewrite(S, [urf_path 'CVHM_strmlnFitted_Endpnt']);
%% Update the coordinates in CVHM_ALLURFS to 3310
S = shaperead([urf_path 'CVHM_strmln_EndpntTA']);
load([urf_path 'CVHM_ALLURFS']);
%%
for ii = 1:length(S)
    URFS(ii,1).X = S(ii,1).X;
    URFS(ii,1).Y = S(ii,1).Y;
end
%% Create a shapefile with the generated wells
Wells = readWells('CVHM_wells3.npsat');
S(size(Wells,1),1).Geometry = 'Point';
for ii = 1:size(Wells,1)
    S(ii,1).Geometry = 'Point';
    S(ii,1).X = Wells(ii,1);
    S(ii,1).Y = Wells(ii,2);
    S(ii,1).Eid = ii-1;
    S(ii,1).top = Wells(ii,3);
    S(ii,1).bot = Wells(ii,4);
    S(ii,1).Q = Wells(ii,5);
end
shapewrite(S,[urf_path 'CVHMwells']);
%%
Wells = shaperead([urf_path 'CVHMwellsTA']);
%% Fit lognormal to urfs
load('/home/giorgk/Documents/UCDAVIS/Mantis/Local/CVHM/CVHM_ALLURFS_TA.mat', 'URFS')
%%
for ii = 1:size(URFS)
    
    
end


