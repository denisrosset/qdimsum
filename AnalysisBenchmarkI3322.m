data3 = load('BenchmarkI3322.mat');
sol3 = 5;
diff3 = data3.objs - sol3;
tbasis3 = cellfun(@(t) t.sampling + t.rank, data3.timings);
tblockdiag3 = cellfun(@(t) t.groupDecomposition + t.blockDiagonalization, data3.timings);
tsolver3 = cellfun(@(t) t.solver, data3.timings);
names = {'{\ttfamily none}' '{\ttfamily reynolds}' '{\ttfamily isotypic}' ...
         '{\ttfamily irreps}' '{\ttfamily blocks}'};
data7 = load('BenchmarkI3322big.mat');
sol7 = 5;
diff7 = data7.objs - sol7;
tbasis7 = cellfun(@(t) t.sampling + t.rank, data7.timings);
tblockdiag7 = cellfun(@(t) t.groupDecomposition + t.blockDiagonalization, data7.timings);
tsolver7 = cellfun(@(t) t.solver, data7.timings);
blockSizes3 = [52 52 26 13 13];
blockSizes7 = [0 244 122 61 61];
i = 1;
str = [names{i} ' & ' num2str(blockSizes3(i))];
str = [str ' & ' num2str(mean(tblockdiag3(i,:)))];
str = [str ' & ' num2str(mean(tbasis3(i,:)))];
str = [str ' & ' num2str(mean(tsolver3(i,:)))];
str = [str ' & ' num2str(mean(abs(diff3(i,:))))];
str = [str ' & & & & & \\ \hline'];
disp(str);
for i = 2:5
    str = [names{i} ' & ' num2str(blockSizes7(i))];
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
