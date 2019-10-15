%% Paths
sim_results_dir = '/home/giorgk/Documents/UCDAVIS/CVHM_NPSAT/Sim_1995m7_1999m12/Results/';
file_prefix = 'WELLURFS_';
nproc = 64;
%% Fit URFs or load fitted
URFparam = fitURF([sim_results_dir file_prefix], nproc);
save([sim_results_dir 'fittedParamURF.mat'],'URFparam');
%%
load([sim_results_dir 'fittedParamURF.mat'],'URFparam');
%% Write wells and exit points as shapefile in order to convert them to the
% TA coordinate system.
clear S
S(2000000,1).Geometry = 'Point';
S(2000000,1).X = 0;
S(2000000,1).Y = 0;
S(2000000,1).Eid = 0;
S(2000000,1).Sid = 0;
cnt = 1;
for ii = 0:nproc-1
    ii
    clear w
    w = load([sim_results_dir file_prefix num2str(ii,'%04d') '.mat']);
    for k = 1:length(w.wellurfs)
        for j = 1:length(w.wellurfs(k,1).DATA)
            S(cnt,1).Geometry = 'Point';
            S(cnt,1).X = w.wellurfs(k,1).DATA(j,1).p_lnd(1);
            S(cnt,1).Y = w.wellurfs(k,1).DATA(j,1).p_lnd(2);
            S(cnt,1).Eid = double(w.wellurfs(k,1).DATA(j,1).Eid);
            S(cnt,1).Sid = double(w.wellurfs(k,1).DATA(j,1).Sid);
            cnt = cnt + 1;
        end
    end
end
S(cnt:end,:) = [];
%% write shapefile
shapewrite(S,[sim_results_dir 'endpnts']);
%% export the shapefile as 3310
% and read it back here
S = shaperead([sim_results_dir 'endpntsTA']);
%% calculate the pixel ij that the points correspond to
pnt = [[S.X]' [S.Y]'];
[II, JJ] = calc_IJ_Mantis( pnt );
%% Find out which streamlines end up on rivers
streams = readStreams([sim_results_dir '../cvhm_1995m7_1999m12_Streams.npsat']);
Scvhm = shaperead([sim_results_dir 'endpnts']);
pnt_cvhm = [[Scvhm.X]' [Scvhm.Y]'];
inRivers = false(size(pnt,1),1);
for ii = 1:length(streams)
    if rem(ii,100) == 0
        display(ii)
    end
    id = find(inpolygon(pnt_cvhm(:,1), pnt_cvhm(:,2), streams(ii,1).poly(:,1),streams(ii,1).poly(:,2)));
    inRivers(id,1) = true;
end
%% print to file for server
DATA = nan(length(S),7); % Eid Sid I J mu std wgh
for ii = 1:length(S)
    if inRivers(ii)
        mu = 0;
        std = 0;
    else
        mu = URFparam.mu(S(ii,1).Eid+1, S(ii,1).Sid+1);
        std = URFparam.std(S(ii,1).Eid+1, S(ii,1).Sid+1);
    end
    DATA(ii,:) = [S(ii,1).Eid+1 S(ii,1).Sid+1 II(ii, 1) JJ(ii,1) ...
        mu std......
        URFparam.wgh(S(ii,1).Eid+1, S(ii,1).Sid+1)];
end
[ii, jj] = find(isnan(DATA));
DATA(ii,:) = [];
fid = fopen([sim_results_dir 'Sim_1995m7_1999m12_URFS.dat'],'w');
fprintf(fid, '%d CVHM_95_99\n', size(DATA,1));
fprintf(fid, '%d %d %d %d %f %f %f\n', DATA');
fclose(fid);
%% Create a well shapefile to convert it to TA
% Run this snippet inside the Simulation folder
wells = readWells('cvhm_1995m7_1999m12_Wells_iter_1.npsat');
for ii = 1:size(wells,1)
   S(ii,1).Geometry = 'Point';
   S(ii,1).X = wells(ii,1);
   S(ii,1).Y = wells(ii,2);
   S(ii,1).Eid = ii;
end
shapewrite(S, 'wellShapefile');
%% Read it 
S = shaperead('wellShapefileTA');
fid = fopen('cvhm_1995m7_1999m12_Wells.dat','w');
fprintf(fid, '%d CVHM_95_99\n', length(S));
fprintf(fid, '%d %f %f\n',[[S.Eid]' [S.X]' [S.Y]']');
fclose(fid);
