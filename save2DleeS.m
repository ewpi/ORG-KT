
function [] = save2DleeS(path,dataToSave)

    fileID = fopen(path,'w');                            % open file to write
    fwrite(fileID,size(dataToSave),'int32','ieee-le');   % save data size 
    fwrite(fileID,dataToSave,'single','ieee-le');        % save data in double precision
    fclose(fileID);                                      % close file

end

