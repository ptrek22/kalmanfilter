flushinput(B);
dataY = [1 2 3 4 5 6.5];
dataX = 0.1*[1 2 3 4 5 6 7 8];
fwrite(B, dataY, 'float32');
while(1)
    fwrite(B, dataX, 'float32');
    arec = fread(B, 8, 'float32')
end