function [bscan] = read3DleeS(path)

    file1 = fopen(path,'r');                      % open file to read
    S = fread(file1,3,'int32', 'ieee-le');        % read 3-dimentional data size
    % 4 bytes x 3 numbers = 12 bytes offset to data	
    data = fread(file1,S(1)*S(3)*S(2),'single','ieee-le');
    % read phase data in single (4 bytes) precision
    fclose('all');                                % close file

    bscan = reshape(data,[S(1) S(2) S(3)]);

end