
function [bscan] = read3DleeD(path)

    file1 = fopen(path,'r');                      % open file to read
    s = fread(file1,3,'int32', 'ieee-le');     % read 3-dimentional data size
    data = fread(file1,s(1)*s(2)*s(3),'double','ieee-le');       % read phase data in double precision
    fclose('all');                                % close file

    bscan = reshape(data,[s(1) s(2) s(3)]);

end