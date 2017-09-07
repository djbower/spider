% Quick and dirty plot for m file output

% Using Octave

% This is truly ugly, intended just to work with the output formats that are there. Once we are done with our diagnostic work, we should have a much better way of examining the (binary) output.

close all; clear; clc; more off

outputFiles = glob('output/*.m'); % Note these will be in lexicographical order (thus 10 comes before 2..)

for i = 1:length(outputFiles)
  source(outputFiles{i});
  w = who();
  disp(w{1})
  data = eval(w{1});                            % hideous - get first variable
  clear(w{1})
  subplot(1,2,1)
  plot(cumsum(data),'-x','Color',0.5*rand(1,3)); hold on; % sum up
  subplot(1,2,2)
  plot(data(2:end),'-x','Color',0.5*rand(1,3)); hold on;  % ignore first point
end
