function distList = nn(ssAC, th)
% Returns the nearest neighbour distribution for trichome cells. ssAC is
% the sum of the concentrations of active complex 1 and 2. th is the
% threshold above which a cell is considered to be a trichome.
   
    T = find(ssAC>=th);
    
    % If only one or two trichome peak(s) found, no nn distances need to be
    % calculated
    if numel(T) <= 2
        distList = NaN;
        return
    end
    
    y = 20 - mod(T-1,20);   % y-coords
    x = 1 + floor((T-1)/(20)) - y/2; % xcoords
    [~, d] = knnsearch([x(:) y(:)], [x(:) y(:)], 'k', 2);
    distList = d(:,2)';
end


