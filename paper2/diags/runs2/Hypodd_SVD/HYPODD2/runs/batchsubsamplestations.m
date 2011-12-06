% This is a batch script to subsample the stations..





station_counts = 1:10;

for i = 1:length(station_counts)
    disp(['Doing directory ', num2str(i), ' of ', num2str(length(station_counts))])
    subsample_hypoDDdata_to_stations(['stat',num2str(station_counts(i))],station_counts(i))
end