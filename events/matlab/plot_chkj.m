function a = plot_chkj(outfname, stfname,time)
% reads the CHKJ96PE2 output and plots
% outname = output name chkj
% stfname = name of stf function
% time = time steps


fid = fopen(outfname)
a = fscanf(fid,'%i %e %e %e', [4 Inf]);
a = a';
fclose(fid)
stf = load(stfname);
[n m] = size(stf);

subplot(4,1,1)
plot(0:time:(n-2)*time, stf(2:end,6))

subplot(4,1,2)
plot(a(:,1)*time,a(:,2))

subplot(4,1,3)
plot(a(:,1)*time, a(:,3))

subplot(4,1,4)
plot(a(:,1)*time, a(:,4))
