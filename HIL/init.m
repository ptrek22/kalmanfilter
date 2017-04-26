B = Bluetooth('HC-05', 1);
B.terminator = 10; %/n
fopen(B);
flushinput(B);
%fread(B, 8, 'float32')
%fwrite(B, 5.12, 'float32')