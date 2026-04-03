function [xp]=perturbFEM(x)
% Perturb a structural configuration by moving material from a
% cell with density==1 to a cell with density==0.
% The perturbation is such that the new cell must share at
% least one edge with an existing cell. This avoids checker-
% board and islands, but not protrusions
% dWo, 2/22/2004
global Nrows Ncolumns


        
matcells=find(x==1);

if mod(Nrows,2)
    indF=ceil(Nrows/2)*Ncolumns;
else
    indF=[Nrows/2  (Nrows/2)+1]*Ncolumns;
end

% reshape x into Nrows x Ncolumns matrix
xm=reshape(x,Ncolumns,Nrows); xm=xm'; xm=flipud(xm);
emptycells=find(xm==0);
XM=[1 zeros(1,Ncolumns) 0; ...
    ones(Nrows,1) xm zeros(Nrows,1); ...
    1 zeros(1,Ncolumns) 0]; % padd boundary with ones and zeros

% find a cell to remove
% indremove=[];
% first search for any islands
% for ind=1:length(matcells);
%    removecell=matcells(ind);
%    rindx= ceil(removecell/Nrows);
%    rindy=removecell-((rindx-1)*Nrows);
%    rindy=rindy+1; rindx=rindx+1; %account for boundary in XM
%    neighbors=sum([ XM(rindy,rindx-1) XM(rindy,rindx+1) XM(rindy-1,rindx) XM(rindy+1,rindx) ]);
%    if neighbors==0
        % have found an island
%        indremove=ind;
%    end
% end
% if isempty(indremove )  
foundcell=0;
while foundcell==0
  indremove=floor(length(matcells)*rand(1))+1; % randomly choose cell index to remove from
   % can't remove material where loads are attached
  if indremove~= indF
    foundcell=1;
    removecell=matcells(indremove);
    rindx=ceil(removecell/Nrows);
    rindy=removecell-((rindx-1)*Nrows);
    rindy=rindy+1; rindx=rindx+1; %account for boundary in XM
  end
end
%end


% find a new empty cell next to an existing cell where to add
foundcell=0; 
while foundcell==0;
    % pick an empty (void) cell at random
    indadd=floor(length(emptycells)*rand(1))+1;
    % check if  addcell is a valid cell
    addcell=emptycells(indadd);
    addcellx= ceil(addcell/Nrows);
    addcelly=addcell-((addcellx-1)*Nrows);
    addcelly=addcelly+1; addcellx=addcellx+1; %account for boundary in XM
    neighbors=sum([ XM(addcelly,addcellx-1) XM(addcelly,addcellx+1) XM(addcelly-1,addcellx) XM(addcelly+1,addcellx) ]);
    if addcellx>2 
        if  neighbors>0 
            foundcell=1;
        end
    else
        if neighbors>1
            foundcell=1;
        end
    end
end

% perturb x
xp=x;
xp(matcells(indremove))=0; % remove cell
removecell=matcells(indremove); 
rindx=mod(removecell,Ncolumns); 
rindy=Nrows-floor(removecell/Ncolumns);

% find symmetrical indremove
% only works for even number of rows! (not checked)
% indremovesymmetric=(rindy-1)*Ncolumns+rindx;
% xp(indremovesymmetric)=0; % remove symmetric cell

% add cells
% convert indadd to original index vector
addcellx=addcellx-1; addcelly=addcelly-1;
indaddo=(Nrows-addcelly)*Ncolumns+addcellx;
%indaddosymmetric=(addcelly-1)*Ncolumns+addcellx;
xp(indaddo)=1; 
%xp(indaddosymmetric)=1;

%keyboard
% end perturbator


   
   

