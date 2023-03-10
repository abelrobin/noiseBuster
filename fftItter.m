function [dataOut] = fftItter(dataIn, lowNotch, highNotch, numbIter)

%first want to see if the data loaded in is a peg file or just a csv
%its going to be better to take in a structure always I think
%also need to get the data from the user
    %frequencey
%need some logical tests to determine if the data is 2 dimensional or just
%contains a single dimension (could accept FID or VUV data).

%next want to put the data into FFTDenoise

%next want to loop over the FFTDenoise function as many times as requested
%by numbIter

%finally want to output the data, I dont see any point of outputting plots
%for the time being.

end