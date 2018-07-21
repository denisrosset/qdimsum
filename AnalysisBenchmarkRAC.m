data3 = load('BenchmarkRAC3.mat');
sol3 = 1/2*(1+1/sqrt(3));
diff3 = data3.objs - sol3;
tbasis3 = cellfun(@(t) t.sampling + t.rank, data3.timings);
tblockdiag3 = cellfun(@(t) t.groupDecomposition + t.blockDiagonalization, data3.timings);
tsolver3 = cellfun(@(t) t.solver, data3.timings);
names = {'{\ttfamily none}' '{\ttfamily reynolds}' '{\ttfamily isotypic}' ...
         '{\ttfamily irreps} (N)' '{\ttfamily blocks} (N)' ...
         '{\ttfamily irreps} (A)' '{\ttfamily blocks} (A)'};

data7 = load('BenchmarkRAC7.mat');
sol7 = 1/2*(1+1/sqrt(7));
diff7 = data7.objs - sol7;
tbasis7 = cellfun(@(t) t.sampling + t.rank, data7.timings);
tblockdiag7 = cellfun(@(t) t.groupDecomposition + t.blockDiagonalization, data7.timings);
tsolver7 = cellfun(@(t) t.solver, data7.timings);
blockSizes3 = [70 70 28 7 7 7 7];
blockSizes7 = [0 750 180 7 7 7 7];
i = 1;
str = [names{i} ' & ' num2str(blockSizes(i))];
str = [str ' & ' num2str(mean(tblockdiag(i,:)))];
str = [str ' & ' num2str(mean(tbasis(i,:)))];
str = [str ' & ' num2str(mean(tsolver(i,:)))];
str = [str ' & ' num2str(mean(abs(diff(i,:))))];
str = [str ' & & & & & \\ \hline'];
disp(str);
for i = 2:7
    str = [names{i} ' & ' num2str(blockSizes3(i))];
    str = [str ' & ' sprintf('$%.2e}$', mean(tblockdiag3(:,i)))];
    str = [str ' & ' sprintf('$%.2e}$', mean(tbasis3(:,i)))];
    str = [str ' & ' sprintf('$%.2e}$', mean(tsolver3(:,i)))];
    str = [str ' & ' sprintf('$%.2e}$', mean(abs(diff3(:,i))))];
    disp(str);
    str = [' & ' num2str(blockSizes7(i))];
    str = [str ' & ' sprintf('$%.2e}$', mean(tblockdiag7(:,i-1)))];
    str = [str ' & ' sprintf('$%.2e}$', mean(tbasis7(:,i-1)))];
    str = [str ' & ' sprintf('$%.2e}$', mean(tsolver7(:,i-1)))];
    str = [str ' & ' sprintf('$%.2e}$', mean(abs(diff7(:,i-1))))];
    str = [str ' \\ \hline'];
    disp(str)
end
