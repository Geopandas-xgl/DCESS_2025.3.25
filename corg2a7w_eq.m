function fhandle=corg2a7w_eq(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb)

  fhandle = @corg;

  function Res2=corg(zo)
  
    [s,x,Rm,C1,D2,D1,C2,Res2]=corg2a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zo);

  end

end

