function savePort(Port, filename, appenddate)
%savePort   Saving signal portrait 'Port' to file 'filename'
%   and, maybe, append current date and time to it.

if nargin < 3, appenddate = 0; end
if appenddate == 1
   filename = [filename '_' datestr(now,'yyyymmdd-HHMMSS')];
end

if nargin < 2, filename = datestr(now,'yyyymmdd-HHMMSS'); end

fid = fopen([filename '.txt'],'w');

for ncell = 1:length(Port)
   fprintf(fid,'cell %d\n',ncell);
   
   for nper = 1:length(Port{ncell}.period)
      fprintf(fid,'period %d\n',nper);
      
      for nwin = 1:length(Port{ncell}.period{nper}.window)
         fprintf(fid,'window %d\n',nwin);
         
         arrsize = length(Port{ncell}.period{nper}.window{nwin}.win_basis);
         fprintf(fid,'size %d\n',arrsize);
         
         fprintf(fid,'win_basis\n');
         for nrow = 1:arrsize
            for ncol = 1:arrsize
               fprintf(fid,'%f ',Port{ncell}.period{nper}.window{nwin}.win_basis(nrow,ncol));
            end
            fprintf(fid,'\n');
         end
         
         fprintf(fid,'svproj\n');
         for nrow = 1:arrsize
            for ncol = 1:arrsize
               fprintf(fid,'%f ',Port{ncell}.period{nper}.window{nwin}.svproj(nrow,ncol));
            end
            fprintf(fid,'\n');
         end
      end
   end
end

fclose(fid);