clc;
% close all;
clear all

% Processing raw data
Fd = 257;
gain = 240;

signal = 21;
strS = ['I' num2str(signal,'%.2d')];
load(       ['D:\Dropbox\Signals\incartdb\' strS '\' strS 'm.mat']);
fid = fopen(['D:\Dropbox\Signals\incartdb\' strS '\annotations.txt']);

val = val/gain;
[Ch,Ts] = size(val);

fgetl(fid);
i = 0;
while ~feof(fid)
   i = i+1;
   line = fgetl(fid);
   mark(i) = str2double(line(15:21));
   annot(1,i) = line(27);
end
fclose(fid);

mark = mark(2:end-1);
annot = annot(2:end-1);

save(['D:\Dropbox\Signals\incartdb\' strS '\' strS 'proc.mat'],...
   'Fd','gain','Ch','Ts','mark','annot','val')

