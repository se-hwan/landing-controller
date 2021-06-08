function  [o1,o2] = plux_2( i1 )

% plux, but only the X -> E,r case.

  o1 = i1(1:3,1:3);
  o2 = -skew_2(o1'*i1(4:6,1:3));

end
