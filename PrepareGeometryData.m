%% ============ MESH FILE based on CVHM grid ==================
% From the complete CVHM database
% https://ca.water.usgs.gov/projects/central-valley/cvhm-database.html
% extract and load the BAS file 
%
bas = shaperead('gis_data/BAS6');
%% Isolate the active grid.
for ii = 1:size(bas,1)
    active = 0;
    for jj = [1 2 3 45 6 7 8 9 10]
        active = active + bas(ii,1).(['lay' num2str(jj) '_act']);
    end
    if active > 0
        bas(ii,1).grid_active = 1;
    else
        bas(ii,1).grid_active = 0;
    end
end
%% delete non active cells
id = find([bas.grid_active]'== 0);
bas(id,:) = [];
%%
shapewrite(bas,'gis_data/BAS_active');
%
%% Make a list of unique nodes and elements
msh_nd = [];
msh_el = [];
for ii = 1:size(bas,1)
    ii
    el = [];
    for jj = 1:4
       px = bas(ii,1).X(jj);
       py = bas(ii,1).Y(jj);
       % test if this node already exists
       if isempty(msh_nd)
           msh_nd = [px py];
           el = 1;
       else
           dst = sqrt((px - msh_nd(:,1)).^2 + (py - msh_nd(:,2)).^2);
           id = find(dst < 0.1);
           if isempty(id)
               msh_nd = [msh_nd;px py];
               el = [el size(msh_nd,1)];
           else
               if length(id) == 1
                   el = [el id];
               else
                   wtf = true;
               end
           end
       end
    end
    msh_el = [msh_el; el];
end
%% write mesh file
writeMeshfile('CVHM_mesh.npsat', msh_nd, msh_el);
%% ============ BOTTOM FILE based on CVHM grid ==================
% load data
[ BOT, info ] = readArcGisASCIIfile( 'gis_data/PP1766_DIS_Layer10_Bot_ASCII.txt' );
% remove data outside grid
BOT(BOT == BOT(1,1)) = nan;
%%
Xgrid = info.xllcorner+info.cellsize/2:...
        info.cellsize:...
        info.xllcorner+info.cellsize/2+info.cellsize*(info.ncols-1);
    
Ygrid = info.yllcorner+info.cellsize/2:...
        info.cellsize:...
        info.yllcorner+info.cellsize/2+info.cellsize*(info.nrows-1);
%%
% find cells with elevation information
[II, JJ] = find(BOT ~= 0 & ~isnan(BOT));
%% expand the nodes by a mile
buff = shaperead('gis_data/CVHM_Mesh_outline_buffer');
%% make unique list of buffer outline points
buff_pnt = [buff(1,1).X(1) buff(1,1).Y(1)];
for ii = 1:length(buff(1,1).X) - 1
    dst = sqrt((buff_pnt(:,1) - buff(1,1).X(ii)).^2 + (buff_pnt(:,2) - buff(1,1).Y(ii)).^2);
    if min(dst) > 0.1
        buff_pnt = [buff_pnt; buff(1,1).X(ii) buff(1,1).Y(ii)];
    end
end
%% Interpolate from the BOT raster the buffer points
Fbot = scatteredInterpolant(Xgrid(JJ)', Ygrid(II)', ...
    BOT(sub2ind(size(BOT),II,JJ)));
Fbot.Method = 'nearest';
Fbot.ExtrapolationMethod = 'nearest';
buff_pnt_Z = Fbot(buff_pnt(:,1), buff_pnt(:,2));
%% write Bottom elevation file
BOTALL = [Xgrid(JJ)', Ygrid(II)', ...
    BOT(sub2ind(size(BOT),II,JJ));...
    buff_pnt buff_pnt_Z];
fid = fopen('CVHM_Bot_elev.npsat','w');
fprintf(fid, 'SCATTERED\n');
fprintf(fid, 'HOR\n');
fprintf(fid, 'SIMPLE\n');
fprintf(fid, '%d %d\n', [size(BOTALL,1) 1]);
fprintf(fid, '%f %f %f\n', BOTALL');
fclose(fid);



