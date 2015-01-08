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
for patient = 2
   strS = ['S' num2str(patient,'%.3d')];
   
   for run = 3
   
      strR = ['R' num2str(run,'%.2d')];
      load(       ['d:\Dropbox\Signals\EEG Motor Movement-Imagery Dataset\Raw\' strS '\' strR '\' strS strR '_edfm.mat']);
      fid = fopen(['d:\Dropbox\Signals\EEG Motor Movement-Imagery Dataset\Raw\' strS '\' strR '\annotations.txt']);
      
      raw = cell(1,ch_num);
      for ch = 1:ch_num
         raw{ch} = val(ch,:);
         raw{ch}(end+1) = 0;
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
      fclose(fid);
      
      len_m(1) = mark(2);
      len_m(2) = mark(3)-mark(2);
      if annot(end) > 0
         mark(end+1) = mark(end)+len_m(2);  % Last movement is longer than others
      else
         mark(end+1) = mark(end)+len_m(1);
      end
      
      mov0 = cell(1,ch_num);
      mov1 = cell(1,ch_num);
      mov2 = cell(1,ch_num);
      for ch = 1:ch_num
         m0=0; m1=0; m2=0;
         for ann = 1:length(annot)
            if annot(ann) == 0
               m0 = m0 + 1;
               mov0{ch}(m0,:) = raw{ch}(mark(ann):mark(ann+1));
            end
            if annot(ann) == 1
               m1 = m1 + 1;
               mov1{ch}(m1,:) = raw{ch}(mark(ann):mark(ann+1));
            end
            if annot(ann) == 2
               m2 = m2 + 1;
               mov2{ch}(m2,:) = raw{ch}(mark(ann):mark(ann+1));
            end
         end
      end


      save(['d:\Dropbox\Signals\EEG Motor Movement-Imagery Dataset\Processed\' strS '\' strR '\' strS strR],...
         'ch_num','channels','Fd','raw','Ts','mov0','mov1','mov2','annot')
   end
end
