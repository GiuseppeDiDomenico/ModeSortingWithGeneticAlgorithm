function T = extractFOS(Z, varargin)
numLeaves = size(Z, 1) + 1;
if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

[p, offset] = parsePAndOffset(nargin,varargin);
[check, color, leafOrder, obslabels, orientation, threshold] = ...
    parseArguments(nargin, varargin, offset, Z, numLeaves);

Z = transz(Z);
[T, Z, leafOrder, numLeaves] = decrowd(numLeaves, Z, p, leafOrder);
end

function [p, offset] = parsePAndOffset(narg,varg)
offset = 0;
p = 30;
if narg == 2
    p = varg{1};
elseif narg > 2
    if isnumeric(varg{1})
        p = varg{1};
        offset = 1;
    end
end
if ~isscalar(p) || p < 0 || p == 1
    error(message('stats:dendrogram:BadP'));
end
end

% Parsing functions
function [check, color, leafOrder, obslabels, orientation, threshold] = ...
    parseArguments(narg, varg, offset, Z, numLeaves)
color = false;
orientation = 't'; %default top
obslabels = [];
threshold = 0.7 * max(Z(:, 3));
leafOrder = [];
check = true;

if narg > 2
    pnames = {'orientation', 'colorthreshold', 'labels', 'reorder', 'checkcrossing'};
    dflts = {orientation, 'default', obslabels, leafOrder, true};
    [orientation, threshold, obslabels, leafOrder, check, setFlag] = ...
        internal.stats.parseArgs(pnames, dflts, varg{1 + offset:end});
    
    verifyCheck(check)
    [color, threshold] = parseColorThreshold(setFlag, threshold, Z);
    obslabels = parseObservationLabels(obslabels, numLeaves);
    leafOrder = parseLeafOrder(leafOrder, numLeaves);
    orientation = parseOrientation(orientation);
end
end

% Tree functions
function Z = transz(Z)
%TRANSZ Translate output of LINKAGE into another format.
%   This is a helper function used by DENDROGRAM and COPHENET.
%   For each node currently labeled numLeaves+k, replace its index by
%   min(i,j) where i and j are the nodes under node numLeaves+k.

%   In LINKAGE, when a new cluster is formed from cluster i & j, it is
%   easier for the latter computation to name the newly formed cluster
%   min(i,j). However, this definition makes it hard to understand
%   the linkage information. We choose to give the newly formed
%   cluster a cluster index M+k, where M is the number of original
%   observation, and k means that this new cluster is the kth cluster
%   to be formed. This helper function converts the M+k indexing into
%   min(i,j) indexing.

numLeaves = size(Z, 1) + 1;

for i = 1:(numLeaves - 1)
    if Z(i, 1) > numLeaves
        Z(i, 1) = traceback(Z, Z(i, 1));
    end
    
    if Z(i, 2) > numLeaves
        Z(i, 2) = traceback(Z, Z(i, 2));
    end
    
    if Z(i, 1) > Z(i, 2)
        Z(i, 1:2) = Z(i, [2 1]);
    end
end
end

function a = traceback(Z, b)
numLeaves = size(Z, 1) + 1;

if Z(b - numLeaves, 1) > numLeaves
    a = traceback(Z, Z(b - numLeaves, 1));
else
    a = Z(b - numLeaves, 1);
end

if Z(b - numLeaves, 2) > numLeaves
    c = traceback(Z, Z(b - numLeaves, 2));
else
    c = Z(b - numLeaves, 2);
end

a = min(a, c);
end

function [T, Z, leafOrder, numLeaves] = decrowd(numLeaves, Z, p, leafOrder)
% If there are more than p nodes, the dendrogram looks crowded.
% The following code will make the last p link nodes into leaf nodes,
% and only these p nodes will be visible.
T = (1:numLeaves)';

if (numLeaves > p) && (p ~= 0)
    Y = Z((numLeaves - p + 1):end, :); % get the last nodes
    [T, W] = getDecrowdedOrdering(T, Y, p);
    T = assignOriginalLeavesToNodes(numLeaves, T, Z, p);
    Z = createSmallerTree(W, Y);
    [leafOrder, numLeaves] = assignLeafOrderToNodes(T, leafOrder, p);
end
end

function [T, W] = getDecrowdedOrdering(T, Y, p)
R = unique(Y(:, 1:2));
Rlp = R(R <= p);
Rgp = R(R > p);
W(Rlp) = Rlp; % use current node number if <=p
W(Rgp) = setdiff(1:p, Rlp); % otherwise get unused numbers <=p
W = W(:);
T(R) = W(R);
end

function T = assignOriginalLeavesToNodes(numLeaves, T, Z, p)
for i = numLeaves - p:-1:1
    T(Z(i, 2)) = T(Z(i, 1));
end
end

function Z = createSmallerTree(W, Y)
Y(:, 1) = W(Y(:, 1));
Y(:, 2) = W(Y(:, 2));
% At this point, it's possible that Z(i,1) < Z(i,i) for some rows.
% The newly formed cluster will always be
% represented by the number in Z(i,1);
Z = Y;
end

function [leafOrder, numLeaves] = assignLeafOrderToNodes(T, leafOrder, p)
if ~isempty(leafOrder)
    leafOrder = T(leafOrder)';
    d = diff(leafOrder);
    d = [1 d];
    leafOrder = leafOrder(d ~= 0);
    
    if numel(leafOrder) ~= p
        error(message('stats:dendrogram:InvalidLeafOrder'));
    end
end
numLeaves = p; % reset the number of nodes to be p (row number = p-1).
end