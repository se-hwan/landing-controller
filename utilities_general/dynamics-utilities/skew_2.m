function  out = skew_2( in )

% skew, but only the mat to vec one.

% if all(size(in)==[3 3])			% do v = skew(A)
  out = 0.5 * [ in(3,2) - in(2,3);
		in(1,3) - in(3,1);
		in(2,1) - in(1,2) ];
% else					% do S = skew(v)
%   out = [  0,    -in(3),  in(2);
% 	   in(3),  0,    -in(1);
% 	  -in(2),  in(1),  0 ];
% end
