
table = readtable('F:/test_MEM_GC_encod.xlsx')
% formula = 'GC~freq0+freq1+(direction|corr_incorrect)+(subject_number|corr_incorrect)+(channel_number|corr_incorrect)';
% formula = 'GC~freq0+freq1+direction+(channel_number|direction)+(subject_number|direction)+(1|subject_number)';
formula = 'GC~freq0+freq1+direction+channel_number+subject_number+(1|subject_number)';

stat = fitlme(table,formula)



F = fitted(stat);
R = response(stat);
figure();
plot(R,F,'rx')
xlabel('Response')
ylabel('Fitted')

figure();
plotResiduals(stat,'fitted')