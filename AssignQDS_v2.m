function [ Q_out, D_out, S_out ] = AssignQDS_v2(F_QD, F_DS, Q_in, D_in, S_in)
if ~isnan(Q_in) && ~isnan(D_in) && ~isnan(S_in)
    Q_out = Q_in;
    D_out = D_in;
    S_out = S_in;
    return;
end

while 1
    if isnan(Q_in)
        Qr = min(F_QD.X) + (max(F_QD.X) - min(F_QD.X))*rand;
    else
        Qr = Q_in;
    end

    if isnan(D_in)
        Dr = min(F_QD.Y) + (max(F_QD.Y) - min(F_QD.Y))*rand;
    else
        Dr = D_in;
    end

    if isnan(S_in)
        Sr = min(F_DS.Y) + (max(F_DS.Y) - min(F_DS.Y))*rand;
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
end
