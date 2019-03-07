% pkg load mapping
% pkg load io
% load Bas shapefile with the active cells only
bas = shaperead('BAS_active');
% load the average stresses
load('AvStresses')
% First remove the extra fields
for ii = [1 2 3 45 6 7 8 9 10]
    bas = rmfield(bas, ['lay' num2str(ii) '_act']);
end

for ii = 1:10
    bas = rmfield(bas, ['strt_hd_' num2str(ii,'%02d')]);
end
R = [bas.ROW]';
C = [bas.COLUMN_]';

for ii = 1:length(bas)
    bas(ii,1).rch_77_99 = AVrch_77_99(bas(ii,1).ROW, bas(ii,1).COLUMN_);
    bas(ii,1).rch_91_03 = AVrch_91_03(bas(ii,1).ROW, bas(ii,1).COLUMN_);
    bas(ii,1).well_77_99 = AVwell_77_99(bas(ii,1).ROW, bas(ii,1).COLUMN_);
    bas(ii,1).well_91_03 = AVwell_91_03(bas(ii,1).ROW, bas(ii,1).COLUMN_);
    
    id = find(STRMS_77_99(:,2) == bas(ii,1).ROW & STRMS_77_99(:,3) == bas(ii,1).COLUMN_);
    if isempty(id)
       bas(ii,1).str_77_99 = 0;
    else
       bas(ii,1).str_77_99 = STRMS_77_99(id,4);
    end
    
    id = find(STRMS_91_03(:,2) == bas(ii,1).ROW & STRMS_91_03(:,3) == bas(ii,1).COLUMN_);
    if isempty(id)
       bas(ii,1).str_91_03 = 0;
    else
       bas(ii,1).str_91_03 = STRMS_91_03(id,4);
    end
end
