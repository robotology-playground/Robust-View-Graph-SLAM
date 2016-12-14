function a=luinv(a)

n=size(a,1);

%BEGIN DECOMPOSITION
% Diagonalization
for k=1:n-1
    for i=k+1:n
        a(i,k)=a(i,k)/a(k,k);
        for j=k+1:n
            a(i,j)=a(i,j)-a(i,k)*a(k,j);
        end
    end
end
%END DECOMPOSITION

for i = 1 : n
    for j = 1 : n
        if i == j
            b(j) = 1;
        else
            b(j) = 0;
        end
    end
    x(1)=b(1);
    for k = 2 : n
        s=0;
        for j = 1 : k-1
            s = s + a(k,j) * x(j);
        end
        x(k) = b(k) - s;
    end
    % back substitution to solve Ux=d
    x(n) = x(n) / a(n, n);
    for k = n-1:-1:1
        s=0;
        for j = k+1 : n
            s = s + a( k, j) * x(j);
        end
        x(k) = (x(k) - s)/a(k,k);
    end
    for j = 1: n
        a(j,i) = x(j);
    end
end