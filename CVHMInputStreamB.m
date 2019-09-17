function CVHMInputStreamB(opt)


% Read stream polygons from houdini.
% Run this part after the animation in houdini file "CVHM_HOU.hipnc" has
% been played. The animation creates one file per river path.
% Here we read the files and create the stream input file
cnt = 1;
clear STRM_POLY
STRM_POLY(1,1).Q = [];
STRM_POLY(1,1).Geo = [];
for ii = 1:53 % check the number of generated files
    fid = fopen([opt.simFolder filesep 'stream_temp_files/temp_' opt.timestring '_out_' num2str(ii) '.txt'], 'r');
    Npoly = textscan(fid,'%d',1);
    for jj = 1:Npoly{1,1}
        C = textscan(fid,'%d %f',1);
        STRM_POLY(cnt,1).Q = C{1,2};
        poly = zeros(C{1,1},2);
        for k = 1:C{1,1}
            C = textscan(fid,'%f %f',1);
            poly(k,:) = [C{1,1} C{1,2}] ;
        end
        poly = poly*100000;
        D = pdist(poly);
        
        if min(D) < 0.5
            while 1
                D = squareform(D);
                break_in = false;
                for r = 1:size(D,1)
                    for c = 1:size(D,2)
                        if r == c
                            continue;
                        else
                            if D(r,c) == 0
                                poly(c,:) = [];
                                break_in = true;
                            end
                            if break_in
                                break;
                            end
                        end
                        if break_in
                            break;
                        end
                    end
                    if break_in
                        break;
                    end
                end
                D = pdist(poly);
                if min(D) > 0.5
                    break;
                end
            end
        end
        STRM_POLY(cnt,1).Geo = poly;
        cnt = cnt + 1;
    end
    fclose(fid);
end
STRM_POLY(find([STRM_POLY.Q]' == 0),:) = [];
% remove the segments without area
dlt_id = [];
for ii = 1:size(STRM_POLY,1)
    Np = size(STRM_POLY(ii,1).Geo,1);
    if Np <=2
        dlt_id = [dlt_id; ii];
    end
end
STRM_POLY(dlt_id,:) = [];

% Write the stream file (finally!)
fid = fopen([opt.simFolder filesep opt.prefix '_' opt.timestring  '_Streams.npsat'],'w');
fprintf(fid, '%d\n', size(STRM_POLY,1));
for ii = 1:size(STRM_POLY,1)
    Np = size(STRM_POLY(ii,1).Geo,1);
    if Np <=2
        continue;
    end
    
    if Np == 3
        A = polyarea(STRM_POLY(ii,1).Geo(:, 1), STRM_POLY(ii,1).Geo(:, 2));
        fprintf(fid, '%d %f\n', [size(STRM_POLY(ii,1).Geo,1) STRM_POLY(ii,1).Q/A]);
        fprintf(fid, '%f %f\n', STRM_POLY(ii,1).Geo');
    elseif Np == 4
        A = polyarea(STRM_POLY(ii,1).Geo([1 2 4 3],1), STRM_POLY(ii,1).Geo([1 2 4 3],2));
        fprintf(fid, '%d %f\n', [size(STRM_POLY(ii,1).Geo,1) STRM_POLY(ii,1).Q/A]);
        fprintf(fid, '%f %f\n', STRM_POLY(ii,1).Geo([1 2 4 3],:)');
    else
        check_that = true;
    end
end
fclose(fid);