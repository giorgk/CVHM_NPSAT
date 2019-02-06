urf_path = '/media/giorgk/RAVEN/CVHM_URFS/';
%%
WellURF = [];
jj = 0;
for ii = 0:49
    temp = readURFs([urf_path 'CVHMwells_' num2str(jj,'%04d') '_' num2str(ii,'%04d') '_streamlines.urfs']);
    WellURF = [WellURF;temp];
end
%%
%WellURF = [];
jj = 0;
for ii = 50:99
    temp = readURFs([urf_path 'CVHMwells_' num2str(jj,'%04d') '_' num2str(ii,'%04d') '_streamlines.urfs']);
    WellURF = [WellURF;temp];
end
%%
%WellURF = [];
jj = 0;
for ii = 100:127
    temp = readURFs([urf_path 'CVHMwells_' num2str(jj,'%04d') '_' num2str(ii,'%04d') '_streamlines.urfs']);
    WellURF = [WellURF;temp];
end
%%
%WellURF = [];
jj = 0;
for jj = 1:2
    for ii = 0:127
        temp = readURFs([urf_path 'CVHMwells_' num2str(jj,'%04d') '_' num2str(ii,'%04d') '_streamlines.urfs']);
        WellURF = [WellURF;temp];
    end
end
%% 
% stuck on '/media/giorgk/New Volume/CVHM_URFS/CVHMwells_0002_0005_streamlines.urfs'
%WellURF = [];
for jj = 2:2
    for ii = 5:127
        temp = readURFs([urf_path 'CVHMwells_' num2str(jj,'%04d') '_' num2str(ii,'%04d') '_streamlines.urfs']);
        WellURF = [WellURF;temp];
    end
end
%% Write age info on shapefile
clear S
S(size(ageURF, 1),1).Geometry = [];
for ii = 1:size(ageURF, 1)
    S(ii,1).Geometry = 'Point';
    S(ii,1).X = ageURF(ii,1);
    S(ii,1).Y = ageURF(ii,2);
    S(ii,1).el = ageURF(ii,3);
    S(ii,1).age = ageURF(ii,4)/365;
end
shapewrite(S, 'strmln_Age')