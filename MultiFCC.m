function MU = MultiFCC(H, K, L)

%   MultiFCC returns multiplicities of(H K L) relections in FCC structure
%   H, K, L must be integers AND they cannot all ZERO
%   If H, K, L < 0, get absolute value


A = [abs(H),abs(K),abs(L)];

Aorder = sort(A,'descend');

h = Aorder(1);
k = Aorder(2);
l = Aorder(3);


    if l > 0
        if k == l % hhh or hkk
            if h == k % M = 8 for hhh
                MU = 8; 
            else % M = 24 for hkk
                MU = 24; 
            end
        else   % hhl or hkl
            if h == k % M = 24 for hhl
                MU = 24; 
            else % M = 48 for hkl
                MU = 48; 
            end
        end

    else % l = 0

        if k == l % M = 6 for h00
            if h > 0
                MU = 6; 
            else
                MU = 0; % return 0 for 000
            end
        else
            if h == k % M = 12 for hh0
                MU = 12;
            else  % M = 12 for hk0
                MU = 24;
            end
        end

    end

end

