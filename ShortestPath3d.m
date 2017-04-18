% Problem3d.m Solution Script
% 23Feb 2017
% *** Remove all headers and white space from input text***
% lengths of Shortest paths from all vertices to vertex i, and
% lengths of shortest paths from i to all vertices

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

%% PROCESS DATA - PART I - FIND SHORTEST PATHS FROM VERTICES TO i
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

% single equality constraint - distance to i = 0; 
% Account for unreachable nodes l & m

Aeq = zeros(3, numberOfNodes);
startNode = double('i') - double('a') + 1;
nodeL = double('l') - double('a') + 1;
nodeM = double('m') - double('a') + 1;

Aeq(1, startNode) = 1;
Aeq(2, nodeL) = 1;
Aeq(3, nodeM) = 1;
beq = [0; 99999;99999];

% Minimize Constraint to Max negative sum of distances
f = -ones(numberOfNodes, 1);

% This is where the magic happens ... call the linprog function to utilize
% the simplex method
 [x, fval, exitflag] = linprog(f, A, b, Aeq, beq);
   
 %% PROCESS DATA - PART II - FIND SHORTEST PATHS FROM I TO ALL OTHER VERTICES
 
 distanceToNodei = x;
 
 % numberOfNodes -> highest numbered node
 numberOfNodes = max([edgeStart, edgeEnd]);
 
 % Build a and b matrices from edgeStart -> edgeEnd
 % A -> number of inequalities by num of nodes
 A = zeros(numel(edgeWeight), numberOfNodes);
 
 for j = 1:numel(edgeWeight)
    A(j, edgeStart(j)) = -1;
    A(j, edgeEnd(j)) = 1;
 end

b = edgeWeight';

% Add constraints < 0
% identity matrix 
A = [A; -eye(numberOfNodes)];
% set zeros
b = [b; zeros(numberOfNodes, 1)];

% single equality constraint
Aeq = zeros(1, numberOfNodes);
startNode = double('i') - double('a') + 1;
Aeq(1, startNode) = 1;
beq = 0;

% Minimize Constraint to Max negative sum of distances
f = -ones(numberOfNodes, 1);

% call linprog
[x, fval, exitflag] = linprog(f, A, b, Aeq, beq);
 
distanceFromNodei = x;
 
 %% COMBINE RESULTS FROM PART I AND PART II
 for i = 1:numberOfNodes
     for j = 1:numberOfNodes
         distFromTo(i,j) = distanceToNodei(i) + distanceFromNodei(j);
         if distFromTo(i,j) > 999
             distFromTo(i,j) = NaN;
         end
     end
 end
 
 % Open the output file and print results
 fid = fopen('Problem3D_Solution.txt', 'w');
 fprintf(fid, 'Problem 3.D Solution:\n\n');
 for i = 1:numberOfNodes
     for j = i:numberOfNodes
        fprintf(fid, 'Shortest Distance from %c to %c passing through i = %2.0f \n',char('a'+i-1), char('a'+j-1), distFromTo(i,j));
     end
 end
 
 fclose(fid);
