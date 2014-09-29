clc;close all;clear all

% Processing raw data
ch_num = 64;
fid = fopen('d:\Dropbox\Signals\EEG Motor Movement-Imagery Dataset\channels.txt');
channels = cell(1,ch_num);
for i = 1:ch_num
   channels{i} = fgetl(fid);
end
fclose(fid);

Fd = 160;
for patient = 1:1
   strS = ['S' num2str(patient,'%.3d')];
   
   for run = 3:3
   
      strR = ['R' num2str(run,'%.2d')];
      load(       ['d:\Dropbox\Signals\EEG Motor Movement-Imagery Dataset\' strS '\' strR '\' strS strR '_edfm.mat']);
      fid = fopen(['d:\Dropbox\Signals\EEG Motor Movement-Imagery Dataset\' strS '\' strR '\annotations.txt']);
      
      raw = cell(1,ch_num);
      for ch = 1:ch_num
         raw{ch} = val(ch,:);
      end
      Ts = length(raw{1});
      
      fgetl(fid);fgetl(fid);fgetl(fid);
      i = 0;
      while ~feof(fid)
         i = i+1;
         line = fgetl(fid);
         mark(i) = str2double(line(30:35)) + 1;  % originally started from 0
         annot(i) = str2double(line(58));
      end
      mark(i+1) = Ts;
      fclose(fid);
      
      mov0 = cell(1,ch_num);
      mov1 = cell(1,ch_num);
      mov2 = cell(1,ch_num);
      for ch = 1:ch_num
         t0=0; t1=0; t2=0;
         for ann = 1:length(annot)-1
            if annot(ann) == 0
               t0 = t0 + 1;
               mov0{ch}(t0,:) = raw{ch}(mark(ann):mark(ann+1));
            end
            if annot(ann) == 1
               t1 = t1 + 1;
               mov1{ch}(t1,:) = raw{ch}(mark(ann):mark(ann+1));
            end
            if annot(ann) == 2
               t2 = t2 + 1;
               mov2{ch}(t2,:) = raw{ch}(mark(ann):mark(ann+1));
            end
         end
      end


      save(['d:\Dropbox\Signals\EEG Motor Movement-Imagery Dataset\Processed\' strS '\' strR '\' strS strR],...
         'ch_num','channels','Fd','raw','Ts','mov0','mov1','mov2')
   end
end
