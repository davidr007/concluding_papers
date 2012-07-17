

Event1 = dlmread('Event1.asc', '\t', 8, 0);
Event2 = dlmread('Event2.asc', '\t', 8, 0);
Event3 = dlmread('Event3.asc', '\t', 8, 0);

figure 
h(1) = plot(117.181, -30.460, 'o')
hold on
h(2) = plot(117.066, -30.537, '^')
h(3) = plot(117.057, -30.498, '+')
h(4) = plot(117.056, -30.494, 's')
legend(h, {'Receiver BK7A','Event1','Event2', 'Event3'})
plot([117.1,117.1521],[-30.53,-30.53],'b-','LineWidth', 2)
text(117.12, -30.527, '5 km')


figure
subplot(3,1,1)
title('Waveforms not normalised')
ind1 = find(Event1(:,1)==0);
plot(Event1(ind1:end,1), Event1(ind1:end,2),'b')
hold on
ind2 = find(Event2(:,1)==0);
plot(Event2(ind2:end,1)+0.16, Event2(ind2:end,2),'r')

subplot(3,1,2)
maxxcorr = seis_xcorr(Event1(ind1:end,2), Event2(ind2:end,2), Event1(ind1:end,1), 2)
plot(maxxcorr(:,1),maxxcorr(:,2))

figure
subplot(3,1,1)
title('Waveforms normalised by maximum value')
ind1 = find(Event1(:,1)==0);
plot(Event1(ind1:end,1), Event1(ind1:end,2)/max(Event1(:,2)),'b')
hold on
ind2 = find(Event2(:,1)==0);
plot(Event2(ind2:end,1)+0.16, Event2(ind2:end,2)/max(Event2(:,2)),'r')

subplot(3,1,2)
maxxcorr = seis_xcorr(Event1(ind1:end,2)/max(Event1(:,2)), Event2(ind2:end,2)/max(Event2(:,2)), Event1(ind1:end,1), 1)
plot(maxxcorr(:,1),maxxcorr(:,2))


figure
plot(Event1(:,1), Event1(:,2),'b')
hold on
plot(Event2(:,1), Event2(:,2),'r')

figure
ind = find(Event1(:,1)==0);
plot(Event1(ind:end,1), Event1(ind:end,2),'b')
hold on
ind = find(Event2(:,1)==0);
plot(Event2(ind:end,1)+0.16, Event2(ind:end,2),'r')