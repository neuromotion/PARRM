function [y,count,c] = uniquecount(x,varargin)
% function [y,count,c] = uniquecount(x,varargin)
%
% varargin passes arguments to sort(x(:),varargin{:})
% 
% y are the unique elements of x, sorted as above
% count(i) is the number of times y(i) occurred in x
% c satisfies x(:) = y(c)

if nargout > 2
	[x,d] = sort(x(:),varargin{:});
else
	x = sort(x(:),varargin{:});
end

m = numel(x);

ndx = [true ; x(1:m-1) ~= x(2:m)];

y = x(ndx);

if nargout > 1
	
	n = numel(y);
	
	count = ones(n,1);
	
	j = 1;
	
	for k = 1:n-1
		
		j = j + 1;
		
		z = 1;
		
		while ~ndx(j)
			j = j + 1;
			z = z + 1;
		end
		
		count(k) = z;

	end

	count(n) = m - j + 1;

	if nargout > 2
		c = cumsum(ndx); 
		r = zeros(m,1); 
		r(d)=1:m; 
		c = c(r);
	end
end
