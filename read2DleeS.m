
function [bscan] = read2DleeS(path)

    file1 = fopen(path,'r');                      % open file to read
    size = fread(file1,2,'int32', 'ieee-le');     % read 3-dimentional data size
    data = fread(file1,size(1)*size(2),'single','ieee-le');       % read phase data in double precision
    fclose('all');                                % close file

    bscan = reshape(data,[size(1) size(2)]);

end