newPackage(
	"Msolve",
	Version => "1.1", 
    	Date => "July 8, 2023",
    	Authors => {{Name => "Martin Helmer", 
		  Email => "mhelmer@ncsu.edu", 
		  HomePage => "http://martin-helmer.com/"}},
    	Headline => "An interface to the msolve package which computes Groebner Basis and does real root isolation",
    	DebuggingMode => false,
	Configuration => {"msolveBinaryFolder"=>"~/msolve-v0.5.0"},
	PackageImports=>{"Elimination","PrimaryDecomposition","Saturation","SegreClasses"}
    	);
--currently setup to use the  configation option to point to the folder where the msolve binary file lives
msolvePath = (options Msolve).Configuration#"msolveBinaryFolder";
if not instance(msolvePath,String) then error "expected configuration option msolveMainFolder to be a string."
export{
    "grobBasis",
    "realSolutions", 
    "rUR",
    "saturateByPoly",
     "Output"
    }
inputOkay=method(TypicalValue=>Boolean);
inputOkay(Ideal):=I->(
    R:=ring I;
    if not instance(R,PolynomialRing) then (print "input must be in a polynomial ring over a feild";return false;);
    kk:=coefficientRing(R);
    if instance(kk,InexactFieldFamily) then(
	print "input must be over the rationals or a finite feild of chactersitic between 2^16 and 2^32";
	return false;
	);
    if char(kk)>0 then (
	if (char(kk)<2^16) or (char(kk)>2^32) then(
	    print "input must be over the rationals or a finite feild of chactersitic between 2^16 and 2^32";
	    return false;
	    );
	);
    return true;
    );
grobBasis=method(TypicalValue=>Matrix);
grobBasis(Ideal):=(I)->(
    if not inputOkay(I) then return 0;
    mIn:=temporaryFileName()|".ms";
    R:=ring I;
    kk:=coefficientRing(R);
    gR:=toString(gens R);
    l1:=substring(1,#gR-2,gR);
    l2:=char R;
    gI:=toString flatten entries gens I;
    if (isField(kk) and (char(kk)==0)) then gI=replace("[)(]","",gI);
    Igens:=replace(" ",""|newline,substring(1,#gI-2,gI));
    inStr:=l1|newline|l2|newline|Igens;
    mIn<<inStr<<close;
    mOut:=temporaryFileName()|".ms";
    callStr:=msolvePath|"/msolve -t "|toString(max(maxAllowableThreads,allowableThreads))|" -g 2 -f "|mIn|" -o "|mOut;
    run(callStr);
    moutStrL:=separate(",",first separate("[]]",last separate("[[]",get(mOut))));
    M2Out:=for s in moutStrL list value(s);
    msolGB:=matrix {M2Out};
    return gens forceGB msolGB;
    );
saturateByPoly=method(TypicalValue=>Matrix);
saturateByPoly(Ideal,RingElement):=(I,f)->(
    if not inputOkay(I) then return 0;
    mIn:=temporaryFileName()|".ms";
    R:=ring I;
    kk:=coefficientRing(R);
    gR:=toString(gens R);
    l1:=substring(1,#gR-2,gR);
    l2:=char R;
    gI:=toString flatten entries gens I;
    if (isField(kk) and (char(kk)==0)) then gI=replace("[)(]","",gI);
    Igens:=replace(" ",""|newline,substring(1,#gI-2,gI));
    inStr:=l1|newline|l2|newline|Igens|","|newline|toString(f);
    mIn<<inStr<<close;
    mOut:=temporaryFileName()|".ms";
    callStr:=msolvePath|"/msolve -f "|mIn|" -S -g 2 -o "|mOut;
    run(callStr);
    moutStrL:=separate(",",first separate("[]]",last separate("[[]",get(mOut))));
    M2Out:=for s in moutStrL list value(s);
    msolGB:=matrix {M2Out};
    return gens forceGB msolGB;
    );

realSolutions=method(TypicalValue=>List,Options => {Output=>"rationalInterval"});
realSolutions(Ideal):=opts->(I)->(
    if not inputOkay(I) then return 0;
    mIn:=temporaryFileName()|".ms";
    R:=ring I;
    kk:=coefficientRing(R);
    gR:=toString(gens R);
    l1:=substring(1,#gR-2,gR);
    l2:=char R;
    gI:=toString flatten entries gens I;
    if (isField(kk) and (char(kk)==0)) then gI=replace("[)(]","",gI);
    Igens:=replace(" ",""|newline,substring(1,#gI-2,gI));
    inStr:=l1|newline|l2|newline|Igens;
    mIn<<inStr<<close;
    mOut:=temporaryFileName()|".ms";
    callStr:=msolvePath|"/msolve -t "|toString(max(maxAllowableThreads,allowableThreads))|" -f "|mIn|" -o "|mOut;
    run(callStr);
    mOutStr:=replace("[]]","}",replace("[[]","{",(separate("[:]",get(mOut)))_0));
    solsp:=value(mOutStr);
    print callStr;
    if solsp_0>0 then (print "Input ideal not zero dimensional, no solutions found."; return 0;);
    if (solsp_1)_0>1 then (print "unexpected msolve output, returning full output"; return solsp;);
    sols:=(solsp_1)_1;
    if opts.Output=="rationalInterval" then return sols;
    if opts.Output=="floatInterval" then return (1.0*sols);
    if opts.Output=="float" then return (for s in sols list(for s1 in s list sum(s1)/2.0));
    );
rUR=method(TypicalValue=>List);
rUR(Ideal):=(I)->(
    if not inputOkay(I) then return 0;
    mIn:=temporaryFileName()|".ms";
    R:=ring I;
    kk:=coefficientRing(R);
    gRm2:=gens R;
    gR:=toString(gRm2);
    l1:=substring(1,#gR-2,gR);
    l2:=char R;
    gI:=toString flatten entries gens I;
    if (isField(kk) and (char(kk)==0)) then gI=replace("[)(]","",gI);
    Igens:=replace(" ",""|newline,substring(1,#gI-2,gI));
    inStr:=l1|newline|l2|newline|Igens;
    mIn<<inStr<<close;
    mOut:=temporaryFileName()|".ms";
    run(msolvePath|"/msolve -t "|toString(max(maxAllowableThreads,allowableThreads))|" -P 2 -f "|mIn|" -o "|mOut);
    mOutStr:=replace("[]]","}",replace("[[]","{",(separate("[:]",get(mOut)))_0));
    solsp:=value replace("[']","\"",mOutStr);
    if first solsp!=0  then (print "Input ideal not zero dimensional, no solutions found."; return 0;);
    lc:=(solsp_1)_4;
    l:=sum(numgens R,i->lc_i*gRm2_i);
    RUR:= new MutableHashTable;
    T:= getSymbol("T");
    S2:=(coefficientRing(R))[T];
    T=first gens S2;
    RUR#"T"=l;
    RUR#"var"=T;
    RUR#"degree"=(solsp_1)_2;
    para:= ((solsp_1)_5)_1;
    W:=sum((para_0)_0+1,i->(T)^i*(((para_0)_1)_i));
    RUR#"findRootsUniPoly"=W;
    RUR#"denominator"=diff(T,W);
    vs:=last para;
    RUR#"numerator"=append(for f in vs list sum((f_0)_0+1,i->T^i*((f_0)_1)_i),T);
    return new HashTable from RUR;
    );

beginDocumentation()
multidoc ///
///	      
	      
TEST ///
-*  
    restart
    needsPackage "Msolve"
    installPackage "Msolve"
*-

R=ZZ/1073741827[x,y]
R=QQ[x,y]
n=2
f = product(n,i->x-i)
g = product(n,i->y-1/(i+1))
I = ideal (f,2*g)
grobBasis(I)
groebnerBasis I
time realSolutions(I)
rur= time rUR(I)
peek oo
time realSolutions(I,Output=>"float")




S = QQ[t12,t13,t14,t23,t24,t34]
S = ZZ/1073741827[t12,t13,t14,t23,t24,t34]
I = ideal(
(t13*t12-t23)*(1-t14)+(t14*t12-t24)*(1-t13) - (t12+1)*(1-t13)*(1-t14), 
(t23*t12-t13)*(1-t24)+(t24*t12-t14)*(1-t23) - (t12+1)*(1-t23)*(1-t24), 
(t14*t34-t13)*(1-t24)+(t24*t34-t23)*(1-t14) - (t34+1)*(1-t14)*(1-t24));

sat=(1-t24)
J1= saturate(I,sat)
J2=ideal saturateByPoly(I,sat)
J1==J2



///
