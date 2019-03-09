pkg load mapping
pkg load io
% load Bas shapefile with the active cells only
bas = shaperead('BAS_active');
% load the average stresses
load('../AvStresses')
% First remove the extra fields
for ii = [1 2 3 45 6 7 8 9 10]
    bas = rmfield(bas, ['lay' num2str(ii) '_act']);
end

for ii = 1:10
    bas = rmfield(bas, ['strt_hd_' num2str(ii,'%02d')]);
end


for ii = 1:length(bas)
    bas(ii).rch_77_99 = AVrch_77_99(bas(ii).ROW, bas(ii).COLUMN_);
    bas(ii).rch_91_03 = AVrch_91_03(bas(ii).ROW, bas(ii).COLUMN_);
    bas(ii).well_77_99 = AVwell_77_99(bas(ii).ROW, bas(ii).COLUMN_);
    bas(ii).well_91_03 = AVwell_91_03(bas(ii).ROW, bas(ii).COLUMN_);
    
    id = find(STRMS_77_99(:,2) == bas(ii).ROW & STRMS_77_99(:,3) == bas(ii).COLUMN_);
    if isempty(id)
       bas(ii).str_77_99 = 0;
    else
       bas(ii).str_77_99 = sum(STRMS_77_99(id,4));
    end
    
    id = find(STRMS_91_03(:,2) == bas(ii).ROW & STRMS_91_03(:,3) == bas(ii).COLUMN_);
    if isempty(id)
       bas(ii).str_91_03 = 0;
    else
       bas(ii).str_91_03 = sum(STRMS_91_03(id,4));
    end
end
