function fhandle=corg1a7w_eq(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb)

  fhandle = @corg;

  function Res1=corg(zo)
  
    [s,x,Rm,C1,D1,Res1]=corg1a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zo);

  end

end

