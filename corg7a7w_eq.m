function fhandle=corg7a7w_eq(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb)

  fhandle = @corg;

  function Res7=corg(zo)
  [s,x,Rm,C1,D7,D6,D5,D4,D3,D2,D1,C2,C3,C4,C5,C6,C7,Res7]=...
      corg7a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zo);
  end

end

