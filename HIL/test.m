flushinput(B);
dataY = [1 2 3 4 5 6.5];
dataX = 0.1*[1 2 3 4 5 6 7 8];
fwrite(B, dataY, 'float32');
fwrite(B, dataX, 'float32');
while(1)
    recieveMatrixBt(B)
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


