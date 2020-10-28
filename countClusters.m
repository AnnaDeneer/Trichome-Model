function clusters = countClusters(ssAC, th, xmax, ymax, D)
% Algorithm to find and count clusters using steady state AC concentrations
% (ssAC = ssAC1+ssAC2). Trichomes are identified by ssAC > th (th:
% threshold).
% Example:
% [t, y] = simModel(k, 0);
% ssAC1 = y(end, ctr+5);
% ssAC2 = y(end, ctr+6);
% ssAC = ssAC1+ssAC2; 
% th = 0.5*max(ssAC);
% c = countClusters(ssAC, th);

T = find(ssAC>=th); % get trichome indices

[i,j] = find(D); % get both x,y coordinates from sparse D matrix
inCluster = zeros(1,ymax*xmax);
for idx = 1 : numel(T)
        %Current trichome neighbor indices:
        nidx = i(j==T(idx) & j~=i);
        if any(ismember(T, nidx)) %check which of neighbours are trichomes         
            inCluster(T(idx)) = 1;
        end
end
clusters = sum(inCluster);