newPackage(
	"Msolve",
	Version => "1.2", 
    	Date => "August 6, 2023",
    	Authors => {{Name => "Martin Helmer", 
		  Email => "mhelmer@ncsu.edu", 
		  HomePage => "http://martin-helmer.com/"}},
    	Headline => "An interface to the msolve package (https://msolve.lip6.fr/) which computes Groebner Basis and does real root isolation",
    	DebuggingMode => false,
	Configuration => {"msolveBinaryFolder"=>"~/msolve-v0.5.0"}
    	);
--currently setup to use the  configation option to point to the folder where the msolve binary file lives
msolvePath = (options Msolve).Configuration#"msolveBinaryFolder";
if not instance(msolvePath,String) then error "expected configuration option msolveMainFolder to be a string."
export{
    "grobBasis",
    "realSolutions", 
    "rUR",
    "saturateByPoly",
    "eliminationIdeal",
    "leadingMonomials",
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
leadingMonomials=method(TypicalValue=>Matrix);
leadingMonomials(Ideal):=(I)->(
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
    callStr:=msolvePath|"/msolve -t "|toString(max(maxAllowableThreads,allowableThreads))|" -g 1 -f "|mIn|" -o "|mOut;
    run(callStr);
    moutStrL:=separate(",",first separate("[]]",last separate("[[]",get(mOut))));
    M2Out:=for s in moutStrL list value(s);
    msolGB:=matrix {M2Out};
    return gens forceGB msolGB;
    );
eliminationIdeal=method(TypicalValue=>Matrix);
eliminationIdeal(Ideal,RingElement):=(I,elimvar)->(
    return eliminationIdeal(I,{elimvar});
    );
eliminationIdeal(RingElement,Ideal):=(elimvar,I)->(
    return eliminationIdeal(I,{elimvar});
    );
eliminationIdeal(List,Ideal):=(elimvars,I)->(
    return eliminationIdeal(I,elimvars);
    );
eliminationIdeal(Ideal,List):=(J,elimvars)->(
    if not inputOkay(J) then return 0;
    mIn:=temporaryFileName()|".ms";
    S:=ring J;
    kk:=coefficientRing(S);
    keepVars:=toList(set(gens(S))-set(elimvars));
    elimNum:=length(elimvars);
    R:=kk[join(elimvars,keepVars),MonomialOrder=>{numgens(S):1}];
    I:=sub(J,R);
    gR:=toString(gens R);
    l1:=substring(1,#gR-2,gR);
    l2:=char R;
    gI:=toString flatten entries gens I;
    if (isField(kk) and (char(kk)==0)) then gI=replace("[)(]","",gI);
    Igens:=replace(" ",""|newline,substring(1,#gI-2,gI));
    inStr:=l1|newline|l2|newline|Igens;
    mIn<<inStr<<close;
    mOut:=temporaryFileName()|".ms";
    callStr:=msolvePath|"/msolve -t "|toString(max(maxAllowableThreads,allowableThreads))|" -e "|toString(elimNum)|" -g 2 -f "|mIn|" -o "|mOut;
    run(callStr);
    moutStrL:=separate(",",first separate("[]]",last separate("[[]",get(mOut))));
    M2Out:=for s in moutStrL list value(s);
    msolGB:=sub(selectInSubring(elimNum,matrix {M2Out}),S);
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

Node 
     Key
     	  Msolve
     Headline
     	  Macaulay2 interface for msolve; computes real solutions, Groebner basis (in GRevLex), saturations, and elimination ideals.
     Description
     	  Text
	      This package provides a Macaulay2 interface for the msolve package (https://msolve.lip6.fr/) developed by Jérémy Berthomieu, Christian Eder, and Mohab Safey El Din.  
	      
	      The package has functions to compute Groebner basis, in GRevLex order only, for ideals with rational or finite feild coefficents. Finite feild chacterisitics must be between 2^16 and 2^32. There are also functions to compute elimination ideals, for ideals with rational or finite feild coefficents. Finite feild chacterisitics must be between 2^16 and 2^32. 
	      
	      The saturation of an ideal by a single polynomial may be computed for ideals with finite feild coefficents, again with charcterisitc between 2^16 and 2^32.
	      
	      For zero dimensional polynomial ideals, with integer or rational coefficents, there are functions to compute all real solutions, and to compute a rational univariante representation of all (complex) solutions.
	      
	      The M2 interface assumes that the binary executable is named "msolve", to tell the M2 package where to find the executable you should use either needsPackage("Msolve", Configuration=>{"msolveBinaryFolder"=>"path_to_folder_with_msolve_binary"}), or installPackage("Msolve", Configuration=>{"msolveBinaryFolder"=>"path_to_folder_with_msolve_binary"}) once and needsPackage("Msolve") on subsequent usages; here "path_to_folder_with_msolve_binary" denotes the folder where your msolve exectuable is saved. E.g. if your msolve binary (named msolve) is in a folder called msolve-v0.5.0 in your home directory then you would use either needsPackage("Msolve", Configuration=>{"msolveBinaryFolder"=>"~/msolve-v0.5.0"}) or installPackage("Msolve", Configuration=>{"msolveBinaryFolder"=>"~/msolve-v0.5.0"}) and needsPackage("Msolve") subseqently   
	      References:
	      [1] msolve: https://msolve.lip6.fr/
	      
Node 
    Key
    	grobBasis
	(grobBasis, Ideal)
    Headline
    	Computes generators of a Groebner basis in GrevLex order.
    Usage
    	grobBasis(I)
    Inputs
    	I:Ideal
	    an ideal in a polynomial ring with GrevLex order and either rational coefficents, integer coefficents, or finite feild coefficents. For a finite feild the charcterisitc must be between 2^16 and 2^32. 
    Outputs
        GB:Matrix
	    a matrix whose columns form a Groebner basis for the input ideal I, in the GrevLex order.    
    Description 
        Text
    	    This functions uses the F4 implmentation in the msolve package to compute a Groebner basis, in GrevLex order, of a polynomial ideal with either rational coefficents or finite feild coefficents with charcterisitc between 2^16 and 2^32. If the input ideal is a polynomial ring with monomial order other than GrevLex a GrevLex basis is returned (and no warning is given). The input ideal may also be given in a ring with integer coefficents, in this case a Groebner basis for the given ideal over the rationals  will be computed, denominators will be cleared, and the output will be a Groebner basis over the rationals in GrevLex order with integer coefficents.
    	Text
	    First an example over a finite feild
	Example
	    R=ZZ/1073741827[z_1..z_3]
	    I=ideal(7*z_1*z_2+5*z_2*z_3+z_3^2+z_1+5*z_3+10,8*z_1^2+13*z_1*z_3+10*z_3^2+z_2+z_1)
	    gB=grobBasis I
	    lT=monomialIdeal leadTerm gB
	    degree lT
	    dim lT	    
	Text
	    Now the same example over the rationals. 
	Example 
	    R=QQ[z_1..z_3]
	    I=ideal(7*z_1*z_2+5*z_2*z_3+z_3^2+z_1+5*z_3+10,8*z_1^2+13*z_1*z_3+10*z_3^2+z_2+z_1)
	    gB=grobBasis I
	    (ideal gB)== ideal(groebnerBasis I)
	    lT=monomialIdeal leadTerm gB
	    degree lT
	    dim lT  
Node 
    Key
    	leadingMonomials
	(leadingMonomials, Ideal)
    Headline
    	Computes the leading monomials of a Groebner basis in GrevLex order.
    Usage
    	leadingMonomials(I)
    Inputs
    	I:Ideal
	    an ideal in a polynomial ring with GrevLex order and either rational coefficents, integer coefficents, or finite feild coefficents. For a finite feild the charcterisitc must be between 2^16 and 2^32. 
    Outputs
        GB:Matrix
	    a matrix whose columns are the leading monomials (of a Groebner basis for) the input ideal I, in the GrevLex order.    
    Description 
        Text
    	    This functions uses the F4 implmentation in the msolve package to compute leading monomials via a Groebner basis, in GrevLex order, of a polynomial ideal with either rational coefficents or finite feild coefficents with charcterisitc between 2^16 and 2^32. If the input ideal is a polynomial ring with monomial order other than GrevLex a GrevLex basis is returned (and no warning is given). The input ideal may also be given in a ring with integer coefficents, in this case a Groebner basis for the given ideal over the rationals  will be computed, denominators will be cleared, and the output will be a Groebner basis over the rationals in GrevLex order with integer coefficents.
    	Text
	    First an example over a finite feild
	Example
	    R=ZZ/1073741827[z_1..z_3]
	    I=ideal(7*z_1*z_2+5*z_2*z_3+z_3^2+z_1+5*z_3+10,8*z_1^2+13*z_1*z_3+10*z_3^2+z_2+z_1)
	    lm=monomialIdeal leadingMonomials I
	    degree lm
	    dim lm	    
	Text
	    Now the same example over the rationals; note over the rationals msolve first computes over a finite feild and when only the leading monomials are asked for the correct leading monomials will be returned but the full Groebner basis over Q will not be computed. Hence if only degree and dimension are desired this command will often be faster that the Groeber basis command.  
	Example 
	    R=QQ[z_1..z_3]
	    I=ideal(7*z_1*z_2+5*z_2*z_3+z_3^2+z_1+5*z_3+10,8*z_1^2+13*z_1*z_3+10*z_3^2+z_2+z_1)
	    lm=monomialIdeal leadingMonomials I
	    lt=monomialIdeal leadTerm groebnerBasis I
	    lm==lt
	    degree lm
	    dim lm

Node 
    Key
    	realSolutions
	(realSolutions, Ideal)
    Headline
    	Uses symbolic methods to compute a certified solution interval for each real solution to a zero dimentional polynomial system.
    Usage
    	realSolutions(I)
    Inputs
    	I:Ideal
	    a zero dimensional ideal in a polynomial ring with either rational or integer coefficents. 
    Outputs
        Sols:List
	    a list of lists, each entry in the list Sol consists of a list repsenting the coordinates of a solution. By default each solution coordinate value is also represetened by a two element list of rational numbers, {a, b}, this means that that coordinate of the solution has a value greater than or equal to a and less than or equal to b. This interval is computed symbolically and its correctness is gaurentted by exact methods.      
    Description 
        Text
    	    This functions uses the msolve package to compute the real solutions to a zero dimensional polynomial ideal with either integer or rational coefficents. The output is a list of lists, each entry in the list Sol consists of a list repsenting the coordinates of a solution. By default each solution coordinate value is also represetened by a two element list of rational numbers, {a, b}, this means that that coordinate of the solution has a value greater than or equal to a and less than or equal to b. This interval is computed symbolically and its correctness is gaurentted by exact methods. Note that using the option Output one may specify the output in terms of either a float inteval with "floatInterval" or an average of the interval endpoints as a single float with "float".      

    	Text
	    First an example over a finite feild
	Example
	    n=3
	    R=QQ[x_1..x_n]
	    f = (x_1-1)*x_1
	    g = (x_2-5/2)*(x_2-1/2)
	    h=(x_3-2)*(x_3-1)
	    I = ideal (f,g,h)
	    sols=realSolutions(I)
	    floatSolsInterval=realSolutions(I,Output=>"floatInterval")
	    floatAproxSols=realSolutions(I,Output=>"float")
	    	    
	Text
	    Note in cases where solutions have multiplicity this is not returned, but the presence of multiplcity also does not reduce accuracy or reliability in any way.   
	Example 
	    n=3
	    R=QQ[x_1..x_n]
	    f = (x_1-1)*x_1^3
	    g = (x_2-5/2)*(x_2-1/2)
	    h=(x_3-2)*(x_3-1)
	    I = ideal (f,g,h)
	    sols=realSolutions(I)
	    floatSolsInterval=realSolutions(I,Output=>"floatInterval")
	    floatAproxSols=realSolutions(I,Output=>"float")
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
leadingMonomials I
time realSolutions(I)
rur= time rUR(I)
peek oo
time realSolutions(I,Output=>"float")


S = QQ[t12,t13,t14,t23,t24,t34]
S = ZZ/1073741827[t12,t13,t14,t23,t24,t34]
gens S
I = ideal(
(t13*t12-t23)*(1-t14)+(t14*t12-t24)*(1-t13) - (t12+1)*(1-t13)*(1-t14), 
(t23*t12-t13)*(1-t24)+(t24*t12-t14)*(1-t23) - (t12+1)*(1-t23)*(1-t24), 
(t14*t34-t13)*(1-t24)+(t24*t34-t23)*(1-t14) - (t34+1)*(1-t14)*(1-t24));

sat=(1-t24)
J1= saturate(I,sat)
J2=ideal saturateByPoly(I,sat)
J1==J2

R = ZZ/1073741827[x,a,b,c,d]
f = x^2+a*x+b
g = x^2+c*x+d
time eliminate(x,ideal(f,g))
time eliminationIdeal(x,ideal(f,g))
time eliminationIdeal(ideal(f,g),x)

///
