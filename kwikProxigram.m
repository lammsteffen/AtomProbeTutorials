function [proxi, binvector, species] = kwikProxigram(pos,interface,xrngORpos,bin)
% calculates a proxigram for the patch 'interface' for the atoms in 'pos'
% with a binwidth of bin. Alternatively, a second pos file can be parsed in
% xrngORpos, which is then the basis of the proxigram. Parsing of multiple
% posfiles through a cell array is allowed.
% percentages are with respect to all RANGED ATOMS!

DEBUG = false;

%%reads a pos file [x,y,z] converted to a Matlab variable and a vertex file
%%[x,y,z] and assigns every atom to the closest vertex.

%vap: m y z mc vert# disttovert distalongnormal( or to line element)
%shiftdistance
%vertex: x y z obj# nx ny nz d A(or l)

addpath('patch_normals','xml_tools');
addpath('./general/resources/pos_tools');
addpath('resources');


%% variable setup

% read pos file if not parsed 
if exist('pos','var')
    if ischar(pos)
        pos = readpos;
    end
else
    pos = readpos;
end


% read interface obj if not parsed
if ischar(interface)
    clear interface;
    [file, path] = uigetfile('*.obj','choose interface file');
    
    obj = read_wobj_v2([path file]);
    interface.vertices = obj.vertices;
    interface.faces = obj.objects{1}.vertices;
end


% choose if species is choosen from xrng or from parsed from sub-pos
if exist('xrngORpos','var')
    
    if ismatrix(xrngORpos)
        speciesPos{1} = xrngORpos;
        species{1} = 'external';
        unranged = [];
        
    elseif iscell(xrngOrpos)
        speciesPos = xrngORpos;
        for sp = 1:length(speciesPos)
            species{sp} = 'external';
        end
        unranged = [];
        
        
    else
        if ischar(xrngORpos)
            xrngORpos = readxrng;
        end
        [speciesPos, species, unranged] = chooseElement(pos,xrngORpos);
    end
    
else
    answ = questdlg('proxigram for raw pos file?','proxigram source data','choose pos file','select species','select species');
    switch answ
        case 'choose pos file'
            
            speciesPos{1} = readpos;
            species{1} = 'external';
            unranged = [];
            
        case 'select species'
            
            xrngORpos = readxrng;
            [speciesPos, species, unranged] = chooseElement(pos,xrngORpos);
    end
    
    
end

if ~exist('bin','var')
    answ = inputdlg('proxigram bin width:');
    bin = str2double(answ);
    
end


%numAtom = length(pos(:,1));
%numVerts = length(interface.vertices);

% distances are calculated along vertex normals.
normals = patchnormals(interface);


%% tessellation and distance calculation
% for overall pos file
% finding closest point for each atomic position
closest = dsearchn(interface.vertices,delaunayn(interface.vertices),pos(:,1:3));
distVec = pos(:,1:3) - interface.vertices(closest,:);
%distance through dot product
dist = sum(normals(closest,:) .* distVec,2);

% for unranged atoms
% finding closest point for each atomic position
closestU = dsearchn(interface.vertices,delaunayn(interface.vertices),unranged(:,1:3));
distVecU = unranged(:,1:3) - interface.vertices(closestU,:);
%distance through dot product
distU = sum(normals(closestU,:) .* distVecU,2);


% calculating bin centers
binvector = linspace(0,1000*bin,1001);
binvector = [fliplr(uminus(binvector(2:end))) binvector];
binvector(binvector<min(dist) | binvector>max(dist)) = [];

% number of atoms per bin
posHist = hist(dist,binvector) - hist(distU,binvector);


%% for element pos files
for sp = 1:length(speciesPos)
    % finding closest point for each atomic position
    closestS = dsearchn(interface.vertices,delaunayn(interface.vertices),speciesPos{sp}(:,1:3));
    
    distVecS = speciesPos{sp}(:,1:3) - interface.vertices(closestS,:);
    
    %distance through dot product
    distS = sum(normals(closestS,:) .* distVecS,2);
    
    % number of atoms per bin
    proxi{sp} = hist(distS,binvector)./posHist;
end


%% plotting
f = figure;
hold all;

for sp = 1:length(speciesPos)
    plot(binvector,proxi{sp}*100);
end

legend(species)

set(gcf,'Name','proximity histogram');
set(gcf,'Color',[1 1 1]);

set(get(gca,'XLabel'),'String','distance [nm]');
set(get(gca,'YLabel'),'String','concentration [%]');

end

















