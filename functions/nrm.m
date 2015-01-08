function [normalized] = nrm(input)

   normalized = input ./ sqrt( input * input' );
   
end