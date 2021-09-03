function [Q, T] = func_part_Schur_factorize(H, SizePart1, SizePart2)
	if SizePart1 == 0
		[Q,T] = schur(H,'complex');
	else
		H1 = H(1:SizePart1,1:SizePart1);
		H12 = H(1:SizePart1,(SizePart1 + 1):end);
		H2 = H((SizePart1 + 1):end,(SizePart1 + 1):end);
		[Q2,T2] = schur(H2,'complex');
		Q = blkdiag(eye(SizePart1),Q2);
		H12T = H12*Q2;
		T = [H1,H12T;zeros(SizePart2,SizePart1),T2];
	end
end