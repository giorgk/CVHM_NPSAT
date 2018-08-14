function id_in = isPoint_inStream(p, STRM, STRMbbox)
% id_in = isPoint_inStream(p, STRM, STRMbbox) tests if the point p is
% inside the polygons defined in the STRM structure. The third input is a
% matrix with length equal to STRM and has 4 columns [xmin xmax ymin ymax].
% This defines the bounding boxes of STRM

id = find(p(1,1) > STRMbbox(:,1) & p(1,1) < STRMbbox(:,2) & ...
          p(1,2) > STRMbbox(:,3) & p(1,2) < STRMbbox(:,4));

id_in = [];
for ii = 1:length(id)
    in = inpolygon(p(1,1), p(1,2), STRM(id(ii),1).poly(:,1), STRM(id(ii),1).poly(:,2));
    id_in = [id_in; ii];
end

