function [ Q_out, D_out, S_out ] = AssignQDS_v2(F_QD, F_DS, Q_in, D_in, S_in)
if ~isnan(Q_in) && ~isnan(D_in) && ~isnan(S_in)
    Q_out = Q_in;
    D_out = D_in;
    S_out = S_in;
    return;
end

cnt = 0;
while 1
    if isnan(Q_in)
        Qr = min(F_QD.X) + (max(F_QD.X) - min(F_QD.X))*rand;
    else
        Qr = Q_in;
    end
    
    % minimum depth is set to 50 ft
    if isnan(D_in)
        dmin = max([min(F_QD.Y) min(F_DS.X) 1.699]);
        dmax = min([max(F_QD.Y) max(F_DS.X)]);
        Dr = dmin  + (dmax - dmin)*rand;
    else
        Dr = D_in;
    end
    
    % minimum screen length 50 ft
    if isnan(S_in)
        slmin = max([min(F_DS.Y) 1.699]);
        slmax = min([max(F_DS.Y) 3.1]);
        Sr = slmin + (slmax - slmin)*rand;
    else
        Sr = S_in;
    end

    r1 = rand;
    r2 = rand;
    %R1 = interp2(F_QD.X, F_QD.Y, F_QD.V, Qr, Dr);
    %R2 = interp2(F_DS.X, F_DS.Y, F_DS.V, Dr, Sr);
    R1 = F_QD.F(Qr, Dr);
    R2 = F_DS.F(Dr, Sr);
    %R1 = F_QD(Qr, Dr);
    %R2 = F_DS(Dr, Sr);

    if (r1 < R1 && r2 < R2)
        Q_out = Qr;
        D_out = Dr;
        S_out = Sr;
        break;
    end
    cnt = cnt + 1;
    if cnt > 1000
        Q_in = nan;
        D_in = nan;
        S_in = nan;
        cnt = 0;
    end
end
