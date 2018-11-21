function connectedComponents = findConnectedComponents(adj)
% Given an adjacency matrix adj, returns the sets of vertices corresponding
% to connected components
%
% For adj = [0 0 1; 0 0 0; 1 0 0], it returns {[1 3] [2]}
    n = size(adj, 1);
    remaining = 1:n;
    cc = {};
    while ~isempty(remaining)
        comp = [];
        test = remaining(1);
        remaining = remaining(2:end);
        while ~isempty(test)
            t = test(1);
            test = test(2:end);
            comp = [comp t];
            newV = remaining(find(adj(t, remaining)));
            test = [test newV];
            remaining = setdiff(remaining, newV);
        end
        cc = horzcat(cc, {comp});
    end
    connectedComponents = cc;
end
