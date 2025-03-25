function fhandle=corg4a7w_eq(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb)

  fhandle = @corg;

  function Res4=corg(zo)
  [s,x,Rm,C1,D4,D3,D2,D1,C2,C3,C4,Res4]=...
      corg4a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zo);
  end

end

