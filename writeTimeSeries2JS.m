function writeTimeSeries2JS(outfile, varname, D, M, Y, data, append)
% writeTimeSeries2JS(outfile, varname, D, M, Y, data)
% writes time series data into a javascript format file
% outfile is the name of the js file
% varname is the name of the variable
% D, M, Y, are the day month and year
% data are the data
% append: if true the is going to add this variable at the end of the file

if append
    fid = fopen([outfile '.js'],'a');
    fprintf(fid,'\n');
else
    fid = fopen([outfile '.js'],'w');
end
fprintf(fid, 'var %s = [\n', varname);
for ii = 1:length(D)
    fprintf(fid, '{ x: new Date(%d,%d,%d), y: %0.5f },\n', [Y(ii) M(ii) D(ii) data(ii)]);
end
fprintf(fid, '];');

fclose(fid);