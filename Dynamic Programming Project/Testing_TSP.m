%=======Traveling Salesman Problem Dynamic Programming Function ===========
%= This function solves the Traveling Salesman Problem (TSP) using Dynamic=
%= programming (DP).                                                      =
%= The function is based on the paper by Held and Karp from 1962. The DP  =
%= is guaranteed to provide the accurate (optimal) result to the TSP, but =
%= the time complexity of this algorithm is O(2^n n^2), which limits the  =
%= the use of this algorithm to 15 cities or less.                        =
%= NOTE: For reasonable runtime, please do not try to calculate a tour of =
%        more than 13 cities. DP is not for large sets of cities.         =
%= Inputs:                                                                =
%       Cities: an n-by-2 (or n-by-3) matrix of city locations. If this   =
%               input is not given, a set of 10 cities will be randomly   =
%               generated.                                                =
%       Dmatrix: This is an n-by-n matrix of the distances between the    =
%               cities. If not given, this matrix will be calculated at   =
%               initialization.                                           =
%= Outputs:                                                               =
%       OptimalTour: a vector with indices of the cities, according to the=
%               the optimal order of cities to visit. It is a closed tour =
%               that will return to city 1, which is assumed to be the    =
%               salesman's origin (depot).                                =
%=      mincost: the total length of the optimal route.                   =
%= Developed by: Elad Kivelevitch                                         =
%= Version: 1.0                                                           =
%= Date: 15 May 2011                                                      =
%==========================================================================

function [OptimalTour_1,mincost_1]=Testing_TSP(cities,Dmatrix,Budget)
% Initialization Process
if nargin==0
    cities=random('unif',0,100,10,2);
    BudgetConstraints=1000;
end
[NumOfCities,dummy]=size(cities);
Primes=primes(NumOfCities*10);
BudgetConstraints=Budget;

if nargin<2 %No Dmatrix used    
    D=diag(inf*ones(1,NumOfCities)); %Assign an inifinite cost to traveling from a city to itself
    for i=1:NumOfCities %Assign the distances between pairs of cities
        for j=i+1:NumOfCities
            D(i,j)=norm(cities(i,:)-cities(j,:));
            D(j,i)=D(i,j);
        end
    end
else
    D=Dmatrix;
    
end
display(D);
NumOfDataSets=1;
for i=2:NumOfCities
    NumOfDataSets=NumOfDataSets+nchoosek(NumOfCities,i);
end
Data(NumOfDataSets).S=[];
Data(NumOfDataSets).l=0;
Data(NumOfDataSets).cost=inf;
Data(NumOfDataSets).pre=[];
Data(NumOfDataSets).m=[];
LookUpTable(NumOfDataSets)=0;
%Define a data structure that holds the following pieces of data we need
%for later. This data structure uses the same notation used in the paper 
%by Held and Karp (1962):
% S - the set of cities in the tour.
% l - the last city visited in the set S. 
% cost - the cost of a tour, which includes all city in S and ends at l.
%In addition, the following data items are used in the dataset for reducing
%runtime:
% Pre - the index of predecessor dataset, i.e. the one with Spre=S-{l}
% m - the city in S-{l} that yielded the lowest cost C(Spre,m)+D(m,l).
% This index will facilitate the generation of the optimal tour without
% further calculations.
Data(1).S=[1];
Data(1).l=1;
Data(1).cost=0;
Data(1).Pre=[];
Data(1).m=[];

%Data Structure to hold data for each stage
Data_1(NumOfCities).a=[];
Data_1(NumOfCities).b=0;

Data_1(1).a=[];
Data_1(1).b=0;

%The minimum cost to travel to one city from starting point (Stage 2)
%Temp_10 placeholder for minimum cost and Temp_11 placeholder for City to
%visit
Temp_10=BudgetConstraints;
Temp_11=[];
for n=2:NumOfCities
    if (D(1,n)<=Temp_10)
        Temp_10=D(1,n);
        Temp_11=[1 n];
    else
        continue;
    end
end
%Input the minimum cost and city visited into the data structure created
%for each stage
Data_1(2).a=Temp_11;
Data_1(2).b=Temp_10;

for s=2:NumOfCities
    Data(s).S=[Data(1).S,s];
    Data(s).l=s;
    Data(s).cost=D(s,1);
    Data(s).Pre=1;
    Data(s).m=1;
    LUT=calcLUT(Data(s).S,s,Primes);
    LookUpTable(s)=LUT;
end
IndexStartPrevStep=2; %index into Data that marks the beginning of the previous step
IndexLastStep=NumOfCities; %index into Data that marks the end of the previous step
CurrentData=IndexLastStep; %index into Data that marks the current dataset
%This is the main Dynamic Programming loop
for s=3:NumOfCities
    OptimalTour=[];
    %generate possible sets with s-1 cities out of the possible N-1 cities
    %(not including city 1)
    TempSets=nchoosek(2:NumOfCities,s-1);
    NumOfSets=size(TempSets);
    for j=1:NumOfSets(1) %iterate through all the sets
        for k=1:NumOfSets(2) %iterate through all the elements in each set
            SminuskSet=[1,TempSets(j,1:k-1),TempSets(j,k+1:NumOfSets(2))]; %this is the set S-{k}               
            candidatecost(2:length(SminuskSet))=inf;
            indices=[];
            for mm=2:length(SminuskSet) %choose a city in S-{k} that will be last
                LUV=calcLUT(SminuskSet,SminuskSet(mm),Primes);
                index=find(LUV==LookUpTable(IndexStartPrevStep:IndexLastStep));  
                index=index+IndexStartPrevStep-1;
                if index==0
                    candidatecost(mm)=inf;                    
                else
                    candidatecost(mm)=Data(index).cost+D(SminuskSet(mm),TempSets(j,k));
                    indices(mm)=index;
                end                
            end
            [mincost,indexcost]=min(candidatecost(2:end));
            CurrentData=CurrentData+1;
            Data(CurrentData).S=[1,TempSets(j,:)];
            Data(CurrentData).l=TempSets(j,k);
            Data(CurrentData).cost=mincost;
            Data(CurrentData).Pre=indices(indexcost+1);
            Data(CurrentData).m=SminuskSet(indexcost+1);
            LookUpTable(CurrentData)=calcLUT(Data(CurrentData).S,TempSets(j,k),Primes);
        end
    end
    IndexStartPrevStep=IndexLastStep+1;
    IndexLastStep=CurrentData;  
nn=0;
candidatecost_1=[];
%Now add the distance back from the last city to city 1
for i=IndexStartPrevStep:IndexLastStep
    nn=nn+1;
    candidatecost_1(nn)=Data(i).cost;
end 
%Find the one that minimizes the total distance
[mincost,indexcost_1]=min(candidatecost_1);
Temp_1=Data(IndexStartPrevStep+indexcost_1-1);
%Generate the optimal tour by traversing back from the last city to its
%predecessors

while ~isempty(Temp_1.Pre)
    OptimalTour=[OptimalTour,Temp_1.l];
    Temp_1=Data(Temp_1.Pre);
end
OptimalTour=[OptimalTour,1];
Data_1(s).a=OptimalTour;
Data_1(s).b=mincost;


end

%Loop to identify which stage yields the minimum cost with the available
%budget and the travel plan for this minimum cost.
Temp_placeholder=[];
Temp_mincost=0;
for m=3:NumOfCities
    if (BudgetConstraints<Data_1(m).b) && (BudgetConstraints>=Data_1(m-1).b)
    Temp_placeholder=Data_1(m-1).a;
    Temp_mincost=Data_1(m-1).b;
    break;
    else 
     Temp_placeholder=Data_1(m).a;
     Temp_mincost=Data_1(m).b;      
    end
end
OptimalTour_1=Temp_placeholder;
mincost_1=Temp_mincost;
     



function LUT=calcLUT(vec,last,Primes)
j=length(vec);
LUT=Primes(last);
for i=2:j
    LUT=LUT*Primes(vec(i));
end
