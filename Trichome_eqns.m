function dydt = Trichome_eqns(t, y, ctr, D, k)
% Trichome patterning model, including TTG1, GL1, GL3, TRY, CPC, AC1
% (complex of TTG1 and GL3), AC2 (complex of GL1 and GL3).
% Input arguments:
% ctr: cell grid indices 
% D: spatial coupling matrix
% k: parameter vector

dydt = zeros(size(y));

TTG1 = y(ctr,:);
GL1 = y(ctr+1,:);
GL3 = y(ctr+2,:);
TRY = y(ctr+3,:);
CPC = y(ctr+4,:);
AC1 = y(ctr+5,:); % TTG1 - GL3
AC2 = y(ctr+6,:); % GL1  - GL3

% TTG1
dydt(ctr,:)   = k(1) - TTG1.*(k(2) + k(3).*GL3) + k(2)*k(4).*(D*TTG1);

% GL1
dydt(ctr+1,:) = k(5) + k(6).*AC2 - GL1.*(k(7) + k(8).*GL3);

% GL3
dydt(ctr+2,:) = k(9) + (k(10)*k(11).*AC1.^2)./(k(11)+AC1.^2) ...
                + (k(12)*k(13).*AC2.^2)./(k(13)+AC2.^2) ...
                - GL3.*(k(14) + k(3).*TTG1 + k(8).*GL1 + ...
                k(15).*TRY + k(16).*CPC);         
            
% TRY
dydt(ctr+3,:) = k(17).*AC1.^2 - TRY.*(k(18) + GL3.*k(15)) ...
                + k(18)*k(19).*(D*TRY);      

% CPC
dydt(ctr+4,:) = k(20).*AC2.^2 - CPC.*(k(21) + k(16).*GL3) + ...
                k(21)*k(22)*(D*CPC);             

% AC1: TTG1-GL3
dydt(ctr+5,:) = k(3).*GL3.*TTG1 - k(23).*AC1;           

% AC2: GL1-GL3
dydt(ctr+6,:) = k(8).*GL3.*GL1 - k(24).*AC2;            
end