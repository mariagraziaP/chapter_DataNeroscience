function [Jin,Jout,Jall] = matching_index_wei2(CIJw)

% computes a new version of the weighted matching index for a weighted,
% directed CIJw matrix

CIJ = double(CIJw>0);

N = size(CIJ,1);

% weighted Jaccard (Ruzicka) similarity - incoming
Jin = zeros(N,N);
for i=1:N
    for j=1:N
        c1 = CIJw(:,i);
        c2 = CIJw(:,j);
        use = setdiff(1:N,[i j]);
        Jin(i,j) = (sum(min(c1(use),c2(use)))./sum(max(c1(use),c2(use))));
    end
end

% weighted Jaccard (Ruzicka) similarity - outgoing
Jout = zeros(N,N);
for i=1:N
    for j=1:N
        c1 = CIJw(i,:)';
        c2 = CIJw(j,:)';
        use = setdiff(1:N,[i j]);
        Jout(i,j) = (sum(min(c1(use),c2(use)))./sum(max(c1(use),c2(use))));
    end
end

% all
Jall = zeros(N,N);
for i=1:N
    for j=1:N
        c1 = [CIJw(:,i) CIJw(i,:)'];
        c2 = [CIJw(:,j) CIJw(j,:)'];
        use = setdiff(1:2*N,[i j i+N j+N]);
        Jall(i,j) = (sum(min(c1(use),c2(use)))./sum(max(c1(use),c2(use))));
    end
end

