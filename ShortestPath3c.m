% Problem3c.m Solution Script
% 23Feb 2017
% *** Remove all headers and white space from input text***
% lengths of Shortest paths from all vertexes to vertex m

clear variables
close all
clc

%% GET DATA FROM FILE
%open file
fid = fopen('Project3Problem3-1.txt');

%read lines while data
rline = fgets(fid);
rowidx = 0;

while ischar(rline)
    % inc count
    rowidx = rowidx + 1;
      
    % splits the string at the specified delimiter
    C = strsplit(rline, ' ');
    
    % convert nodes to indexes using ascii codes (a = 1, b = 2, etc.)
    edgeStart(rowidx) = double(C{1}) - double('a') + 1;
    edgeEnd(rowidx) = double(C{2}) - double('a') + 1;
    edgeWeight(rowidx) = str2num(C{3});
    %disp(double(C{1}) - double('a') + 1);
    
    % go for the next line
    rline = fgetl(fid);
end

fclose(fid);

%% PROCESS DATA - FIND SHORTEST PATHS
% Shortest paths from a to all
numberOfNodes = max([edgeStart, edgeEnd]);

% Build A and B matrices from end to start
% Size A is num of inequal by num of nodes - numel is num of elements
A = zeros(numel(edgeWeight), numberOfNodes);

% reverse the direction of the edges by swapping 1's:
% previously: edgeStart(j) = -1) and edgeEnd(j)= 1
for j = 1:numel(edgeWeight)
    A(j, edgeStart(j)) = 1;
    A(j, edgeEnd(j)) = -1;
end
b = edgeWeight';

% Add constraints < 0
% identity matrix 
A = [A; -eye(numberOfNodes)];
% set zeros
b = [b; zeros(numberOfNodes, 1)];

% single equality constraint - distance to m = 0; 
Aeq = zeros(1, numberOfNodes);
startNode = double('m') - double('a') + 1
Aeq(1, startNode) = 1;
beq = 0;

% Minimize Constraint to Max negative sum of distances
f = -ones(numberOfNodes, 1);

% This is where the magic happens ... call the linprog function to utilize
% the simplex method
 [x, fval, exitflag] = linprog(f, A, b, Aeq, beq);
 
 % Open the output file
 fid = fopen('Problem3C_Solution.txt', 'w');
 
 fprintf(fid, 'Problem 3 C solution:\n\n');
 
 % for all numbers in x, print the results
 for j = 1:numel(x)
     fprintf(fid, 'Distance from m to %c = %2.0f \n', char('a' + j - 1), x(j));
 end
 
fclose(fid);