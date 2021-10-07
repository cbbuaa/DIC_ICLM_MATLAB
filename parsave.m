function parsave(dataFile,disp,strain,ZNCC,iterNum,Params)
% the reason why I put the save function here is to facilitate the usage of
% the parallel computaion: parfor.
% Author: Bin Chen;
% E-mail: binchen@kth.se
% Update: 2021-06-04

save(dataFile,'disp','strain','ZNCC','iterNum','Params');

end