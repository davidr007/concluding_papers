function eps = tot_fracvar(a,mean_medium);
meana = mean(a(:));
eps = sqrt(mean(a(:).^2)/mean_medium^2);


