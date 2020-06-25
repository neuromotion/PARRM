function Prediction = SlidingFourierRegression(Times,Values,Period,Support,NumBasis,EvalTimes)

if nargin < 6 || isempty(EvalTimes), EvalTimes = Times; end

if numel(Support)==1
    Support = [-Support; 0; 0;Support];
elseif numel(Support)==2
    Support = [-flipud(Support(:));Support(1:2)];
elseif numel(Support)==3 && Support(2) == 0
    Support = [Support(1);0;0;Support(3)];
elseif numel(Support) ~= 4
    error('Support incorrectly sized')
end
%if ~issorted(Support), error('Support must be sorted'); end

[Times,sortNdx] = sort(Times(:));
TimesN = numel(Times);
if numel(Values) ~= TimesN, error('Times and Values must be the same size'); end
Values = Values(sortNdx); if size(Values,1) ~= TimesN, Values = Values(:); end

[Prediction,sortNdx] = sort(EvalTimes(:));

t1 = 1;
t2 = 1;
r = max(abs(Support));

% suppress non-invertible warnings
warning off MATLAB:singularMatrix
warning off MATLAB:nearlySingularMatrix

for k = 1:numel(Prediction)

    t = Prediction(k);
    a = t + Support(1);
    b = t + Support(4);
    
    while t1 <= TimesN && Times(t1) < a, t1 = t1 + 1; end
    t2 = max(t2,t1);
    while t2 <= TimesN && Times(t2) <= b, t2 = t2 + 1; end
    
    n = t2-t1;
    if n > 0
        V = Times(t1:t2-1);
        
        X = ones(n,2*NumBasis+1);
        
        for j = 1:NumBasis
            jV = (2*pi*j/Period)*V;
            X(:,2*j) = sin(jV);
            X(:,2*j+1) = cos(jV);
        end
        
        V = V-t;
        W = sqrt((1 - (V/r).^2)).*(V <= Support(2) | V >= Support(3));
        
        X = X.*W;
                
        beta = (X.'*X)\X.'*(Values(t1:t2-1).*W);
        
        p = beta(1);
        for j = 1:NumBasis
            jt = (2*pi*j/Period)*t;
            p = p + beta(2*j)*sin(jt) + beta(2*j+1)*cos(jt);
        end
        
%         windowSize=50;
%         [sortedTimes,I]=sort(mod(Times(t1:t2-1),Period));
%         %allTimes=repmat(sortedTimes,3,1);
%         tempValues=Values(Times(t1:t2-1));
%         allSamps=repmat(tempValues(I),3,1);
%         allWeights=repmat(W(I),3,1);
%         movingAvg=movsum(allWeights.*allSamps,windowSize)./movsum(allWeights,windowSize);
%         movingAvg=movingAvg(length(I)+1:2*length(I));
        
        Prediction(k) = p;
        %Prediction(k)=movingAvg(sortedTimes==mod(t,Period));
    else
        Prediction(k) = NaN;
    end

end

% turn warnings back on
warning on MATLAB:singularMatrix
warning on MATLAB:nearlySingularMatrix

Prediction = reshape(Prediction(sortNdx),size(EvalTimes));
