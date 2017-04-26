function [ A] = recieveMatrixBt(B)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
size = fread(B,2,'uint8')';
A = reshape(fread(B, size(1)*size(2), 'float32'),size);
end

