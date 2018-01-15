z=xlsread("C:\Users\amithmayya\Downloads\Dynamic Optimization\Project\Flight Price Data.xlsx");
z(1:1+size(z,1):end) = Inf;
cities=ones(size(z,1),2);
[a,b]=Testing_TSP(cities,z,79);
