%% 
% This scripts updates the well tops and bottom based on a different initial 
% estimation. Usually you have run the model once and have obtained a water
% table surface. The water table surface is printed a point cloud file
% which we will use it to update the initial elevation and the well screens
% so that the are under the water table
%
%% make init_surf from parallel run
nproc = 64;
Elev_new = [];
for ii = 0:nproc-1
    fid = fopen(['output/Ref7/cvhm_top_000_' num2str(ii,'%04d') '.xyz'],'r');
    Np = fscanf(fid, '%d',1);
    temp = fscanf(fid, '%f',Np*4);
    fclose(fid);
    temp = reshape(temp, 4, Np)';
    Elev_new = [Elev_new;temp];
end
%% Write the assembled new elevation into one file
% This will have duplicalted points. Therefore create a scatter interpolant
% using the code snippet bellow to take care of this
fid = fopen('output/init_surf_ref7.xyz','w');
fprintf(fid, '%d\n', length(FnewElev.Values));
fprintf(fid, '%f %f %f\n', [FnewElev.Points FnewElev.Values]');
fclose(fid);

%% Read new water table as point cloud
fid = fopen('output/init_surf_0.xyz','r');
Np = fscanf(fid, '%d',1);
temp = fscanf(fid, '%f',Np*4);
fclose(fid);
Elev_new = reshape(temp, 4, Np)';
%% Read initial water table
topelev = read_Scattered('CVHM_top_elev.npsat',2);
%% plot the points
plot3(Elev_new(:,1), Elev_new(:,2), Elev_new(:,4),'.r')
hold on
%plot3(Elev_new(:,1), Elev_new(:,2), Elev_new(:,3),'.k')
plot3(topelev.p(:,1), topelev.p(:,2), topelev.v(:,1),'.b')
%% try removing the extreme points
Elev_new_mod = Elev_new;
%Elev_new_mod(Elev_new_mod(:,4) > 225,:) = [];
%dst = sqrt((Elev_new_mod(:,1) - Elev_new_mod(18498,1)).^2 + (Elev_new_mod(:,2) - Elev_new_mod(18498,2)).^2);
%Elev_new_mod(dst < 4700,:) = [];

plot3(Elev_new_mod(:,1), Elev_new_mod(:,2), Elev_new_mod(:,4),'.r')
hold on
%plot3(Elev_new(:,1), Elev_new(:,2), Elev_new(:,3),'.k')
plot3(topelev.p(:,1), topelev.p(:,2), topelev.v(:,1),'.b')
%% Generate an interpolation function out of the new modified elevation
FnewElev = scatteredInterpolant(Elev_new_mod(:,1), Elev_new_mod(:,2), Elev_new_mod(:,4));
FnewElev.Method = 'nearest';
FnewElev.ExtrapolationMethod = 'nearest';
%% Print a new top for the model
newp = FnewElev(topelev.p(:,1), topelev.p(:,2));
writeScatteredData('CVHM_top_elev5.npsat', ...
                   struct('PDIM',2,'TYPE','HOR','MODE','SIMPLE'), ...
                   [topelev.p(:,1) topelev.p(:,2) newp]);
%}
%% Update the wells using the new water table
%% read data
% wells
%
wells = readWells('CVHM_wells.npsat');
% top 
topelev = read_Scattered('CVHM_top_elev3.npsat',2);
Ftop = scatteredInterpolant(topelev.p(:,1), topelev.p(:,2), topelev.v(:,1));
% Old top this will help identifying the depth 
Oldtopelev = read_Scattered('CVHM_top_elev.npsat',2);
Ftop_old = scatteredInterpolant(Oldtopelev.p(:,1), Oldtopelev.p(:,2), Oldtopelev.v(:,1));
% Bottom
botelev = read_Scattered('CVHM_Bot_elev.npsat',2);
Fbot = scatteredInterpolant(botelev.p(:,1), botelev.p(:,2), botelev.v(:,1));
%}
%% update well top and bottom
% update those that their top is above the water table
for ii = 1:size(wells,1)
    xw = wells(ii,1);
    yw = wells(ii,2);
    tw = wells(ii,3);
    bw = wells(ii,4);
    
    % dont allow very short wells
    sl = tw - bw;
    if sl < 6.096
        sl = 6.096;
    end
    
    wt_n = Ftop(xw,yw); % new water table
    wt_o = Ftop_old(xw,yw); % old water table
    bs = Fbot(xw, yw); % base of aquifer
    
    % if the screen length is longer than the aquifer depth, 
    % make the screen length shorter 
    aq_depth = wt_n - bs;
    if sl + 20 > aq_depth
        sl = aq_depth - 20;
    end
    
    % This is the offset in the water table between 2 consecutive runs
    % negative means the water table dropped so we have to move the well
    % lower
    dw = wt_n - wt_o;
    % if the top of the screen was found higher than the new water table
    % then lower by an extra offset
    extra_dw = 0;
    if tw > wt_n
        extra_dw = -30;
    end
    
    % the new well top is the old plus the offset
    tw_n = tw + dw + extra_dw;
    bw_n = tw_n - sl;
    
    % make sure its not too close to water table
    if wt_n - tw_n < 10
        tw_n = wt_n - 10;
        bw_n = tw_n - sl;
    end
    
    % if the bottom of the well was lowered bellow the base move the we;;
    % screen a bit higher
    if bw_n < bs
        dup = 10 + bs - bw_n;
        if tw_n + dup + 10 < wt_n
            tw_n = tw_n + dup;
            bw_n = bw_n + dup;
        else
            wtf = true
        end
    end
    wells(ii,3) = tw_n;
    wells(ii,4) = bw_n;
    %end
end
%% Print the updated wells
fid = fopen('CVHM_wells3.npsat','w');
fprintf(fid, '%d\n', size(wells, 1));
fprintf(fid, '%f %f %f %f %f\n', wells');
fclose(fid);