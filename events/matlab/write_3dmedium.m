function write_3dmedium(fname,vs,dd)
% a is a 3d matrix of shear modulus
[nz,nx,ny] = size(vs);


count = 1
fid = fopen(fname,'wb'); 
for nzn = 1:nz
    for nyn = 1: nx
        for nxn = 1: ny
            disp(['Now doing nxn= ', num2str(nxn), ...
                ' nyn= ', num2str(nyn), ...
                ' nzn= ', num2str(nzn), ': count = ', num2str(count)]);
            fwrite(fid,nxn, 'integer*4');
            fwrite(fid,nyn, 'integer*4');
            fwrite(fid,nzn, 'integer*4');
            fwrite(fid,vs(nzn,nxn,nyn),'real*4');
            fwrite(fid,vs(nzn,nxn,nyn)/1.65,'real*4');
            fwrite(fid,dd,'real*4');
            count = count+1;
        end 
    end
end
fclose(fid)
