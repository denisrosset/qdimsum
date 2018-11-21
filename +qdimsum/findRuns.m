function runs = findRuns(vector, settings)
% Given a vector of ascending values, returns ranges of indices
% [start(i) ... end(i)] such that vec(end(i))-vec(start(i)) <= settings.blockDiagEigTol
    runs = {};
    start = 1;
    for i = 1:length(vector) + 1
        if (i == length(vector) + 1) || vector(i) - vector(start) > settings.blockDiagEigTol
            runs = horzcat(runs, {start:i-1});
            start = i;
        end
    end
    if ~isequal(settings.blockDiagEigHist, [])
        for i = 1:length(runs)
            run = runs{i};
            diff = vector(run(2:end)) - vector(run(1:end-1));
            settings.blockDiagEigHist.register(diff);
        end
    end
end
