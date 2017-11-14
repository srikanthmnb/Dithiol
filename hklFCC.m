function M_hklSmOrder = hklFCC( Smax )

% hklFC: Calculate the allowed Bragg reflection (hkl) in FCC and the
% corresponding multiplicities
% -------------------------------------------------------------------
% Input: Smax, maximum h^2+k^2+l^2 = S;
% Output: a matrix of 5 columns (h, k, l, S, multiplicity) 
%         with increasing order of S
% Example: hklFCC(16)
%         h     k     l     S     Mu
%         1     1     1     3     8
%         2     0     0     4     6
%         2     2     0     8    12
%         3     1     1    11    24
%         2     2     2    12     8
%         4     0     0    16     6
%         ...
% -------------------------------------------------------------------
% H. Zhang @ Sept 13, 2016


    % Calculation allowed Bragg reflection in FCC
    
    % if all h k l are even
    Me_hkl = [];
    h = 0;
    k = 0;
    l = 0;
    Ae_S = [];
    S = h^2+k^2+l^2;
    Ae_Mu = [];
    while h^2 <= Smax

        h = h+2;
        for k=0:2:h
            for l=0:2:k
                S = h^2+k^2+l^2;
                if S<= Smax
                    Me_hkl = [Me_hkl; h, k, l];
                    Ae_S = [Ae_S,S];
                    Mu =  MultiFCC(h,k,l); %Use function MultiFCC to get Mu
                    Ae_Mu = [Ae_Mu, Mu];
                end
            end
        end
    end

    % if all h k l are odd
    h = 1;
    k = 1;
    l = 1;
    Mo_hkl = [h,k,l];

    S = h^2+k^2+l^2;
    Ao_S = [S];

    Mu =  MultiFCC(h,k,l);
    Ao_Mu = [Mu];

    while h^2 <= Smax

        h = h+2;
        for k=1:2:h
            for l=1:2:k
                S = h^2+k^2+l^2;
                if S<= Smax
                    Mo_hkl = [Mo_hkl; h, k, l];
                    Ao_S = [Ao_S,S];
                    Mu =  MultiFCC(h,k,l);
                    Ao_Mu = [Ao_Mu, Mu];

                end
            end
        end
    end

    M_hklSm = [Me_hkl,Ae_S',Ae_Mu';Mo_hkl,Ao_S',Ao_Mu'];
    M_hklSmOrder = sortrows(M_hklSm,4);

end

