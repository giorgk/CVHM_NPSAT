function [ Q_out, D_out, SL_out ] = AssignQDS(S, Q_in, D_in, SL_in )

ww = 0.15;

    if ~isnan(Q_in)
        if Q_in < 10^interp1q(S.fQ, S.xQ, S.lowQ)
            r = S.lowQ + ww*rand;
            Q_in = 10^interp1q(S.fQ, S.xQ, r);
        end
        if Q_in > 10^interp1q(S.fQ, S.xQ, S.uppQ)
            r = S.uppQ - ww + ww*rand;
            Q_in = 10^interp1q(S.fQ, S.xQ, r);
        end
    end
    if ~isnan(D_in)
        if D_in < 10^interp1q(S.fD, S.xD, S.lowD)
            r = S.lowD + ww*rand;
            D_in = 10^interp1q(S.fD, S.xD, r);
        end
        if D_in > 10^interp1q(S.fD, S.xD, S.uppD)
            r = S.uppD - ww + ww*rand;
            D_in = 10^interp1q(S.fD, S.xD, r);
        end
    end
    if ~isnan(SL_in)
        if SL_in < 10^interp1q(S.fS, S.xS, S.lowS)
            r = S.lowS + ww*rand;
            SL_in = 10^interp1q(S.fS, S.xS, r);
        end
        if SL_in > 10^interp1q(S.fS, S.xS, S.uppS)
            r = S.uppS - ww + ww*rand;
            SL_in = 10^interp1q(S.fS, S.xS, r);
        end
    end

S.w = 0.2;

    if isnan(Q_in)
        if isnan(D_in)
            if isnan(SL_in) % all unknown
                Q_out = DSnan(S);
                D_out = Q_Snan(S, Q_out);
                SL_out = QDknown(S, Q_out, D_out);
            else % SL known
                Q_out = SL_Dnan(S,SL_in);
                D_out = QSLknown(S, Q_out, SL_in);
                SL_out = SL_in;
            end
        else
            if isnan(SL_in) % D known only
                Q_out = D_SLnan(S, D_in);
                D_out = D_in;
                SL_out = QDknown(S, Q_out, D_out);
            else % D and SL known
                Q_out = DSLknown(S,D_in, SL_in);
                D_out = D_in;
                SL_out = SL_in;
            end
        end
    else
        if isnan(D_in)
            if isnan(SL_in) % Q known only
                Q_out = Q_in;
                D_out = Q_Snan(S, Q_out);
            else % Q and SL known
                Q_out = Q_in;
                D_out = QSLknown(S, Q_out, SL_in);
                SL_out = SL_in;
            end
        else
            if isnan(SL_in) % Q and D known
                Q_out = Q_in;
                D_out = D_in;
                SL_out = QDknown(S, Q_out, D_out);
                
            else
                Q_out = Q_in;
                D_out = D_in;
                SL_out = SL_in;
            end
        end
    end
end

function rr = randQ(S)
    rr = Stats.lowQ + (S.uppQ - S.lowQ)*rand;
end

function rr = randD(S)
    rr = S.lowD + (S.uppD - S.lowD)*rand;
end

function rr = randS(S)
    rr = S.lowS + (S.uppS - S.lowS)*rand;
end

function Q = DSnan(S)
    r = randQ(S);
    Q = 10^interp1q(S.fQ, S.xQ,r);
end

function D = Q_Snan(S, Q)
    % find the ecdf of the given Q
    rQ = interp1q(S.xQ, S.fQ, log10(Q));
    % generate a random number around that
    w = S.w;
    rD = rQ-w + 2*w*rand;
    D = 10^interp1q(S.fD, S.xD, rD);
end

function SL = QDknown(S, Q, D)
    rQ = interp1q(S.xQ, S.fQ, log10(Q));
    rD = interp1q(S.xD, S.fD, log10(D));
    w = S.w;
    rSL = (rQ + rD)/2 - w + 2*w*rand;
    SL = 10^interp1q(S.fS, S.xS, rSL);
end

function Q = SL_Dnan(S,SL)
    rSL = interp1q(S.xS, S.fS, log10(SL));
    w = S.w;
    rQ = rSL - w + 2*w*rand;
    Q = 10^interp1q(S.fQ, S.xQ, rQ);
end

function D = QSLknown(S, Q, SL)
    rQ = interp1q(S.xQ, S.fQ, log10(Q));
    rSL = interp1q(S.xS, S.fS, log10(SL));
    w = S.w;
    rD = (rQ + rSL)/2 - w + 2*w*rand;
    D = 10^interp1q(S.fD, S.xD, rD);
end

function Q = D_SLnan(S, D)
    rD = interp1q(S.xD, S.fD, log10(D));
    w = S.w;
    rQ = rD-w + 2*w*rand;
    Q = 10^interp1q(S.fQ, S.xQ, rQ);
end

function Q = DSLknown(S,D, SL)
    rD = interp1q(S.xD, S.fD, log10(D));
    rSL = interp1q(S.xS, S.fS, log10(SL));
    w = S.w;
    rQ = (rD + rSL)/2 - w + 2*w*rand;
    Q = 10^interp1q(S.fQ, S.xQ, rQ);
end