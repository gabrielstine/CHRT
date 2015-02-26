function closeParallel()
%DELETEPOOL close current matlab parallel pool if exist

% Copyright 2014 Jian Wang

delete(gcp('nocreate'));

end