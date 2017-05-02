flushinput(B);
dataY = [2 0 0 1 0 0];

fwrite(B, dataY, 'float32');
pause(1);
while(1)
    waitforbuttonpress;
    fwrite(B, dataY, 'float32');
    A =  recieveMatrixBt(B)
end

% ak = recieveMatrixBt(B);
% pk = recieveMatrixBt(B);
% wk = recieveMatrixBt(B);
% q = recieveMatrixBt(B);
% Pap = recieveMatrixBt(B);
% T1 = recieveMatrixBt(B);
% T2 = recieveMatrixBt(B);
% T3 = recieveMatrixBt(B);
% res = ak*pk*ak' + wk*q*wk';


