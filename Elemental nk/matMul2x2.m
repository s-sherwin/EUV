function AB = matMul2x2(A,B)
    s = size(A);
    A = reshape(A,2,2,[]);
    B = reshape(B,2,2,[]);
    AB = zeros(size(A));
    
    AB(1,1,:) = A(1,1,:).*B(1,1,:) + A(1,2,:).*B(2,1,:);
    AB(1,2,:) = A(1,1,:).*B(1,2,:) + A(1,2,:).*B(2,2,:);
    AB(2,1,:) = A(2,1,:).*B(1,1,:) + A(2,2,:).*B(2,1,:);
    AB(2,2,:) = A(2,1,:).*B(1,2,:) + A(2,2,:).*B(2,2,:);
    
    AB = reshape(AB,s);
end