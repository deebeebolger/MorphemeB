
% Finds the indices of the electrodes of interest.
[y,x] = ismember({AllChans.chaninfo(1:GPcfg.Channum).labels},GPcfg.eoi);
[y1,x1] = sort(x(y));
I = find(x);
Inds = I(x1);

% Now can plot the pvalues over time for those electrodes.
% Helps to decide which time windows present significant differences.
t = cell(length(Inds),1);
%figure;
for counter = 1:length(Inds)
%     subplot(3,3,counter)
%     plot(Time, pvalue_corr{Inds(counter),1},'.','MarkerSize',8);
%     hold on
%     h2=hline(0.05,'r','');
    t{counter,1} = Time([pvalue_corr{Inds(counter),1} <= 0.05]);
end


