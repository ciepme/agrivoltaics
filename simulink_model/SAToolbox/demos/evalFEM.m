function [compliance]=evalFEM(density);
% [compliance]=evalFEM(density)
% 
global Nrows Ncolumns
Nrows=8;
Ncolumns=9;

%fig_option=2;
fig_option=0;

% make nodes
counter=0;
for indrows=1:Nrows+1
    for indcolumns=1:Ncolumns+1
        counter=counter+1;
        xyz(counter,:)=[indcolumns-1 indrows-1 0];
    end
end

%make elements
counter=0;
for indrows=1:Nrows
    for indcolumns=1:Ncolumns
    counter=counter+1;
    ni(counter,:)=[counter indcolumns+(indrows-1)*(Ncolumns+1) indcolumns+1+(indrows-1)*(Ncolumns+1) ...
            indcolumns+1+(indrows)*(Ncolumns+1) indcolumns+(indrows)*(Ncolumns+1) counter];
end
end

%specify boundary conditions

bc=zeros(size(xyz,1),6);
for indrows=1:Nrows+1
    bc([(2+(indrows-1)*(Ncolumns+1)):(Ncolumns+1+(indrows-1)*(Ncolumns+1))],[1:2])=ones(Ncolumns,2);
end
    
% number of elements
nelt=size(ni,1);

% initialize m and k
[bco,ndof]=bcond(bc);
k=zeros(ndof,ndof);
m=k;

%define material properties
thick_min = 1.e-5;
thick_max = 0.03;
thickness=density.*(thick_max-thick_min) + thick_min;
thickness=thickness';
prop=repmat([7.1e9 2*7.1e9 2700 0.01 0],nelt,1);
%keyboard
prop(:,4)=thickness;

% compute mass and stiffness matrix
[ko,mo]=plate(ni,xyz,prop,bco,k,m);

%plot structure
if(fig_option==2) 
global labels
labels=[1 1];
viewfem(xyz,ni,4);
hold on
view([0 0 1])
axis off
axis equal
end

% apply loads
if mod(Nrows,2)
F=[(floor(Nrows/2)+1)*(Ncolumns+1) 2 1 -1000; ...
   (floor(Nrows/2)+2)*(Ncolumns+1) 2 1 -1000]; %tip load
else
 F=[(floor(Nrows/2)+1)*(Ncolumns+1) 2 1 -2000];
end

f=loads(F,bco);

% compute displacements
u=ko\f;
dout=disps(u,bco,1);

if(fig_option==2) 
blowup=5000;
xyzd=xyz+blowup*dout;
labels=[1 2];
viewfem(xyzd,ni,4);
end

% compute compliance
comply=sum(dot(f,u));

% Scale the objective function
compliance = 1000*comply;


function [bco,ndof] = bcond(bc)
%BCOND	bcond(bc) determines equation or degree of freedom numbers
%			given fixity inputs for nodes
%
%	bc is a (:,6) array of fixity indicaters for each node
%		bc(i,j)>0 means that degree of freedom will be active
%			in the model
%		otherwise the degree of freedom is assigned the value 0
%	ndof is the total number of degrees of freedom
%	bco is the (:,6) array of degree of freedom numbers
%
% Modified 3/29/95 by Bob Norton:
% Changed bco=zeros(bc) to bco=zeros(size(bc)) to accomodate
% a change by Matlab.
[n,m]=size(bc);
bco=zeros(size(bc));
ndof=0;
for i=1:n
	for j=1:m
		if bc(i,j)>0
			ndof=ndof+1;
			bc(i,j)=ndof;
		end
	end
end
bco = bc;

function [ko,mo]=plate(ni,xyz,prop,bc,k,m)
%PLATE  [ko,mo]=plate(ni,xyz,prop,bc,k,m)
%   computes and assembles the stiffness and mass matrices
%   for plate elements with bending and membrane stiffness.
%
%     ni is (:,5) nodal incidence array, where each row describes
%     an element.
%     ni(:,1) is the element number
%     ni(:,2) is the first node number
%     ni(:,3) is the second node number
%     ni(:,4) is the third node number
%     ni(:,5) is the fourth node number (quads only, 0 for triangles)
%     ni(:,6) is the property number
%
%     xyz is the (:,3) array of nodal coordinates
%
%     prop is a (:,5) array of plate properties
%     prop(:,1) is the modulus of elasticity, e
%     prop(:,2) is the shear modulus, g
%     prop(:,3) is the mass density, mass per unit volume
%     prop(:,4) is the plate thickness
%     prop(:,5) is the non-structural mass per unit area
%
%   This element formulation is taken from the Q-19 element
%   formulated by Clough and Felippa as summarized in
%   AFFDL-TR-74-120, Parts 1 and 2.  The adaptation to
%   the IMOS environment was done by Bob Norton, JPL.

% Copyright (C) 1992, California Institute of Technology.
% U.S. Government Sponsorship under NASA Contract NAS7-918
% is acknowledged.
%
% Revision 00001

%History
%  31Aug93 jmelody:   copied from Stuart Shaklan
%  31Aug93 jmelody:   node numbering capabilities added
%  17Jun94 jmelody:   removed element type from calls to nifix
%

[nxyz,mxyz]=size(xyz);
if mxyz == 4    %node numbering
  index=nodesort(xyz);
  xyz=xyz(:,2:4);               %make xyz as if no nodal numbering
  ni=nifix(index,ni);
end

[n,temp]=size(ni);
ipermq=[2 3 4 1];
iperm=[2 3 1];
mfr=[3 3 3 2 2];
c=zeros(3,3);
t=zeros(1,9);
t1=zeros(1,3);
t2=zeros(1,3);
t3=zeros(1,3);
td1=zeros(1,13);
td2=zeros(1,13);
td3=zeros(1,9);
loc=zeros(1,5);
x=zeros(1,5); y=x; z=x;
nc=zeros(1,3); a=nc; b=nc; u=nc;
p=zeros(21,12);
st=zeros(12,12);
q=zeros(3,6);
g=zeros(1,21);
ht=zeros(1,3); tx=ht; ty=ht;
tr1=zeros(1,9); tr2=tr1; tr3=tr1;

for iel=1:n
  if ni(iel,1) > 0
    np=ni(iel,6);
% Set up the material properties for this element
    [temp,nsm]=size(prop);
    if nsm>4
      nsm=prop(np,5);
    else
      nsm=0;
    end
    h=prop(np,4);
    e=prop(np,1);
    gg=prop(np,2);
    rho=prop(np,3);
    xnu=e/(2.0*gg)-1.0;
    c(1,1)=e/(1.0-xnu*xnu);
    c(1,2)=xnu*c(1,1);
    c(1,3)=0.0;
    c(2,1)=c(1,2);
    c(2,2)=c(1,1);
    c(2,3)=0.0;
    c(3,1)=0.0;
    c(3,2)=0.0;
    c(3,3)=gg;
% Set up the plate geometry
    x=zeros(1,5);
    y=x;
    z=x;
    for i=1:3
      x(i)=xyz(ni(iel,i+1),1);
      y(i)=xyz(ni(iel,i+1),2);
      z(i)=xyz(ni(iel,i+1),3);
    end
%    if ni(iel,5) ~= ni(iel,4)          % if (elt 5 .neq. elt 4)  SBSSBSBBSBBS
	 if ni(iel,5) ~= 0
      x(4)=xyz(ni(iel,5),1);
      y(4)=xyz(ni(iel,5),2);
      z(4)=xyz(ni(iel,5),3);
      x(5)=0.25*(x(1)+x(2)+x(3)+x(4));
      y(5)=0.25*(y(1)+y(2)+y(3)+y(4));
      z(5)=0.25*(z(1)+z(2)+z(3)+z(4));
      n=4;
      nen=4;
      id=zeros(1,24); 
      id(1:6)=bc(ni(iel,2),1:6);
      id(7:12)=bc(ni(iel,3),1:6);
      id(13:18)=bc(ni(iel,4),1:6);
      id(19:24)=bc(ni(iel,5),1:6);
      s=zeros(30,30);
      ek=zeros(24,24);
      em=ek;
    else
      n=1;
      nen=3;
      id=zeros(1,18);
      id(1:6)=bc(ni(iel,2),1:6);
      id(7:12)=bc(ni(iel,3),1:6);
      id(13:18)=bc(ni(iel,4),1:6);
      ek=zeros(18,18);
      s=ek;
      em=ek;
    end
    n3=2*nen-3;
    nef=6*nen;
    nsf=nef+6*(nen-3);
% Compute direction cosine matrix t0 of local element system 
    x1=x(2)+x(3)-x(n)-x(1);
    y1=y(2)+y(3)-y(n)-y(1);
    z1=z(2)+z(3)-z(n)-z(1);
    x2=x(3)+x(n)-x(1)-x(2);
    y2=y(3)+y(n)-y(1)-y(2);
    z2=z(3)+z(n)-z(1)-z(2);
    s1=x1^2+y1^2+z1^2;
    cc=(x1*x2+y1*y2+z1*z2)/s1;
    x2=x2 - cc*x1;
    y2=y2 - cc*y1;
    z2=z2 - cc*z1;
    s1=sqrt (s1);
    s2=sqrt (x2^2+y2^2+z2^2);
    x1=x1/s1;
    y1=y1/s1;
    z1=z1/s1;
    x2=x2/s2;
    y2=y2/s2;
    z2=z2/s2;
    t(1)=x1;
    t(2)=x2;
    t(3)=y1*z2-y2*z1;
    t(4)=y1;
    t(5)=y2;
    t(6)=z1*x2-z2*x1;
    t(7)=z1;
    t(8)=z2;
    t(9)=x1*y2-x2*y1;
% Loop over the triangle components
    for nt=1:n
      n1=nt;
      n2=ipermq(n1);
      nc(1)=n1;
      nc(2)=n2;
      nc(3)=n3;
      nod=3;
% Compute direction cosines of local triangle system
% and the triangle projections a,b onto it.
      a1=x(n1)-x(n3);
      b1=y(n1)-y(n3);
      c1=z(n1)-z(n3);
      a2=x(n2)-x(n3);
      b2=y(n2)-y(n3);
      c2=z(n2)-z(n3);
      t3(1)=b1*c2-b2*c1;
      t3(2)=c1*a2-c2*a1;
      t3(3)=a1*b2-a2*b1;
      ss=sqrt (t3(1)^2+t3(2)^2+t3(3)^2);
      t3(1)=t3(1)/ss;
      t3(2)=t3(2)/ss;
      t3(3)=t3(3)/ss;
      t1(1)=t3(3)*t(5)-t3(2)*t(8);
      t1(2)=t3(1)*t(8)-t3(3)*t(2);
      t1(3)=t3(2)*t(2)-t3(1)*t(5);
      ss=sqrt (t1(1)^2+t1(2)^2+t1(3)^2);
      t1(1)=t1(1)/ss;
      t1(2)=t1(2)/ss;
      t1(3)=t1(3)/ss;
      t2(1)=t1(3)*t3(2)-t1(2)*t3(3);
      t2(2)=t1(1)*t3(3)-t1(3)*t3(1);
      t2(3)=t1(2)*t3(1)-t1(1)*t3(2);
      a(1)=-t1(1)*a2-t1(2)*b2-t1(3)*c2;
      a(2)= t1(1)*a1+t1(2)*b1+t1(3)*c1;
      b(1)= t2(1)*a2+t2(2)*b2+t2(3)*c2;
      b(2)=-t2(1)*a1-t2(2)*b1-t2(3)*c1;
      a(3)=-a(1)-a(2);
      b(3)=-b(1)-b(2);
% Set up inputs for triangle subroutines
      for i=1:3
	l=nc(i);
	loc(i)=6*(l-1);
      end
% Form transformations between element and nodal systems.
%
%     first for translation
%
      for i=1:3
	j=i + 3;
	kk=i + 6;
	td1(i)=t1(i);
	td1(j)=t1(i);
	td2(i)=t2(i);
	td2(j)=t2(i);
	td3(i)=t3(i);
	td3(j)=t3(i);
	if nen == 4
	  ci=t(i);
	  cj=t(j);
	  ck=t(kk);
	  td1(kk)=t1(1)*ci    + t1(2)*cj    + t1(3)*ck;
	  td2(kk)=t2(1)*ci    + t2(2)*cj    + t2(3)*ck;
	  td3(kk)=t3(1)*ci    + t3(2)*cj    + t3(3)*ck;
	else
	  td1(kk)=t1(i);
	  td2(kk)=t2(i);
	  td3(kk)=t3(i);
	end 
      end 
%
%     now for rotation
%
      for i=1:3
	j=i + 3;
	kk=i + 6;
	tr1(i)=t1(i);
	tr1(j)=t1(i);
	tr2(i)=t2(i);
	tr2(j)=t2(i);
	tr3(i)=t3(i);
	tr3(j)=t3(i) ;
	if nen == 4
	  ci=t(i);
	  cj=t(j);
	  ck=t(kk);
	  tr1(kk)=t1(1)*ci    + t1(2)*cj    + t1(3)*ck;
	  tr2(kk)=t2(1)*ci    + t2(2)*cj    + t2(3)*ck;
	  tr3(kk)=t3(1)*ci    + t3(2)*cj    + t3(3)*ck;
	else
	  tr1(kk)=t1(i);
	  tr2(kk)=t2(i);
	  tr3(kk)=t3(i);
	end 
      end
      for i=7:8
	td1(i+3)=td1(i);
	td1(i+5)=td1(i);
	td2(i+3)=td2(i);
	td2(i+5)=td2(i);
      end 
      loc(4)=nsf + 3*(n2-1);
      loc(5)=nsf + 3*(n1-1);
      n4=loc(4) + 3;
      n5=loc(5) + 3;
%
%     Membrane contribution 
%
      ndf=6;
      area=a(3)*b(2)-a(2)*b(3);
      fac=h/(2.*area);
      c11=c(1,1)*fac;
      c22=c(2,2)*fac;
      c33=c(3,3)*fac;
      c12=c(1,2)*fac;
      c13=c(1,3)*fac;
      c23=c(2,3)*fac;
      for j=1:3;
	l=j+j;
	for i=1:j
	  kk=i + i;
	  aa=a(i)*a(j);
	  bb=b(i)*b(j);
	  ab=a(i)*b(j);
	  ba=b(i)*a(j);
	  aba=ab+ba;
	  st(kk-1,l-1)=c11*bb + c33*aa + c13*aba;
	  st(kk  ,l  )=c22*aa + c33*bb + c23*aba;
	  cab=c13*bb + c23*aa;
	  st(kk-1,l  )=c12*ba + c33*ab + cab;
	  st(kk  ,l-1)=c12*ab + c33*ba + cab;
	end
      end
      for i=2:ndf
	for j=1:i
	  st(i,j)=st(j,i);
	end 
      end 
      lt=0;
      for jj=1:nod
	j=jj + jj;
	mm=loc(jj);
	ll=mfr(jj);
	for l=1:ll
	  mm=mm + 1;
	  lt=lt + 1;
	  c1=td1(lt);
	  c2=td2(lt);
	  kt=0;
	  for ii=1:jj
	    i=ii + ii;
	    kk=mfr(ii);
	    if ii == jj  kk=l; end
	    h1=st(i-1,j-1)*c1 + st(i-1,j)*c2;
	    h2=st(i  ,j-1)*c1 + st(i  ,j)*c2;
	    n=loc(ii);
	    for k1=1:kk
	      n=n + 1;
	      kt=kt + 1;
	      sq=s(n,mm) + td1(kt)*h1 + td2(kt)*h2;
	      s(n,mm)=sq;
	      s(mm,n)=sq;
	    end 
	  end 
	end 
      end 
%
%     Plate bending contribution 
%
      ndf=9;
      area=a(3)*b(2)-a(2)*b(3);
      fac=h^3*area/864.;
      for i=1:3
	j=iperm(i);
	xx=a(i)^2+b(i)^2;
	u(i)=-(a(i)*a(j)+b(i)*b(j))/xx;
	xx=sqrt(xx);
	yy=2.*area/xx;
	ht(i)= 2.*yy;
	tx(i)= yy*a(i)/xx;
	ty(i)=-yy*b(i)/xx;
	a1=a(i)/area;
	a2=a(j)/area;
	b1=b(i)/area;
	b2=b(j)/area;
	q(1,i) =b1*b1;
	q(2,i) =a1*a1;
	q(3,i) =2.*a1*b1;
	q(1,i+3)=2.*b1*b2;
	q(2,i+3)=2.*a1*a2;
	q(3,i+3)=2.*(a1*b2+a2*b1);
      end 
      for i=1:3
	j=iperm(i);
	k1=iperm(j);
	ii=3*i;
	jj=3*j;
	kk=3*k1;
	a1=a(i);
	a2=a(j);
	a3=a(k1);
	b1=b(i);
	b2=b(j);
	b3=b(k1);
	u1=u(i);
	u2=u(j);
	u3=u(k1);
	w1=1.-u1;
	w2=1.-u2;
	w3=1.-u3;
	b1d=b1 + b1;
	b2d=b2 + b2;
	b3d=b3 + b3;
	a1d=a1 + a1;
	a2d=a2 + a2;
	a3d=a3 + a3;
	c21=b1-b3*u3         + tx(k1);
	c22=-b1d+b2*w2+b3*u3 + tx(j)-tx(k1);
	c31=a1-a3*u3         + ty(k1);
	c22=-b1d+b2*w2+b3*u3 + tx(j)-tx(k1);
	c31=a1-a3*u3         + ty(k1);
	c32=-a1d+a2*w2+a3*u3 + ty(j)-ty(k1);
	c51=b3*w3-b2         + tx(k1);
	c52=b2d-b3*w3-b1*u1  + tx(i)-tx(k1);
	c61=a3*w3-a2         + ty(k1);
	c62=a2d-a3*w3-a1*u1  + ty(i)-ty(k1);
	c81=b3-b2d-b2*u2     + tx(j);
	c82=b1d-b3+b1*w1     + tx(i);
	c91=a3-a2d-a2*u2     + ty(j);
	c92=a1d-a3+a1*w1     + ty(i);
	u37=7.*u3;
	w27=7.*w2;
	w24=4.*w2;
	u34=4.*u3;
	c1=54.+w27;
	c2=54.+u37;
	c3=15.+w24;
	txs=tx(j)+tx(k1);
	tys=ty(j)+ty(k1);
	for n=1:3
	  l=6*(i-1) + n;
	  q11=q(n,i);
	  q22=q(n,j);
	  q33=q(n,k1);
	  q12=q(n,i+3);
	  q23=q(n,j+3);
	  q31=q(n,k1+3);
	  q2333=q23-q33;
	  q3133=q31-q33;
	  p(l   ,ii-2)=6.*(-q11+w2*q33+u3*q2333);
	  p(l   ,ii-1)=c21*q23+c22*q33-b3d*q12+b2d*q31;
	  p(l   ,ii  )=c31*q23+c32*q33-a3d*q12+a2d*q31;
	  p(l   ,jj-2)=6.*(q22+w3*q2333);
	  p(l   ,jj-1)=c51*q2333+b3d*q22;
	  p(l   ,jj  )=c61*q2333+a3d*q22;
	  p(l   ,kk-2)=6.*(1.+u2)*q33;
	  p(l   ,kk-1)=c81*q33;
	  p(l   ,kk  )=c91*q33;
	  p(l   ,i+9 )=0.;
	  p(l   ,j+9 )=ht(j)*q33;
	  p(l   ,k1+9)=ht(k1)*q2333;
	  p(l+3 ,ii-2)=6.*(q11+u3*q3133);
	  p(l+3 ,ii-1)=c21*q3133-b3d*q11;
	  p(l+3 ,ii  )=c31*q3133-a3d*q11;
	  p(l+3 ,jj-2)=6.*(-q22+u1*q33+w3*q3133);
	  p(l+3 ,jj-1)=c51*q31+c52*q33+b3d*q12-b1d*q23;
	  p(l+3 ,jj  )=c61*q31+c62*q33+a3d*q12-a1d*q23;
	  p(l+3 ,kk-2)=6.*(1.+w1)*q33;
	  p(l+3 ,kk-1)=c82*q33;
	  p(l+3 ,kk  )=c92*q33;
	  p(l+3 ,i+9 )=ht(i)*q33;
	  p(l+3 ,j+9 )=0.;
	  p(l+3 ,k1+9 )=ht(k1)*q3133;
	  p(n+18,ii-2)=2.*(q11+u3*q12+w2*q31);
	  p(n+18,kk-1)=((b1d-b2d)*q33+c82*q23+c81*q31)/3.;
	  p(n+18,kk  )=((a1d-a2d)*q33+c92*q23+c91*q31)/3.;
	  p(n+18,k1+9 )=ht(k1)*q12/3.;
	end 
      end 
      for j=1:ndf
	for l=1:3
	  ii=l;
	  kk=l + 18;
	  pp3=p(kk,j);
	  g(kk)=0.;
	  for n=1:3
	    i=iperm(n);
	    jj=ii + 3;
	    pp1=p(ii,j);
	    pp2=p(jj,j);
	    sum=pp1 + pp2 + pp3;
	    g1=sum + pp1;
	    g2=sum + pp2;
	    g3=sum + pp3;
	    g(ii)=g1;
	    g(jj)=g2;
	    g(kk)=g3 + g(kk);
	    ii=ii + 6;
	  end 
	end 
	for n=1:3:19
	  g1=g(n);
	  g2=g(n+1);
	  g3=g(n+2);
	  g(n) =c(1,1)*g1 + c(1,2)*g2 + c(1,3)*g3;
	  g(n+1)=c(1,2)*g1 + c(2,2)*g2 + c(2,3)*g3;
	  g(n+2)=c(1,3)*g1 + c(2,3)*g2 + c(3,3)*g3;
	end 
	for i=1:j
	  xx=0.;
	  for n=1:21
	    xx=xx + g(n)*p(n,i);
	  end 
	  xx=xx*fac;
	  st(i,j)=xx;
	  st(j,i)=xx;
	end 
      end 
      for jj=1:3
	jt=3*jj-3;
	j=jt + 1;
	for l=1:6
	  mm=loc(jj)+l;
	  l3=l - 3;
	  if l3 <= 0
	    c3=td3(jt+l);
	  end
	end
	for ii=1:jj
	  it=3*ii-3;
	  i=it+1;
	  kk=6;
	  for l=1:6
	    if ii == jj 
	      kk=l;
	    end
	    mm=loc(jj)+l;
	    l3=l-3;
	    if l3 <= 0
	      c3=td3(jt+l);
	      h1=st(i  ,j)*c3;
	      h2=st(i+1,j)*c3;
	      h3=st(i+2,j)*c3;
	    else
	      c1=tr1(jt+l3);
	      c2=tr2(jt+l3);
	      h1=st(i  ,j+1)*c1 + st(i  ,j+2)*c2;
	      h2=st(i+1,j+1)*c1 + st(i+1,j+2)*c2;
	      h3=st(i+2,j+1)*c1 + st(i+2,j+2)*c2;
	    end 
	    nn=loc(ii);
	    for k4=1:kk
	      nn=nn + 1;
	      k3=k4 - 3;
	      k1=it + k4;
	      k2=it + k3;
	      if k3 <= 0
		sq=s(nn,mm) + td3(k1)*h1;
	      else
		sq=s(nn,mm) + tr1(k2)*h2 + tr2(k2)*h3;
	      end 
	      s(nn,mm)=sq;
	      s(mm,nn)=sq;
	    end 
	  end
	end 
      end 
    end 
%
%     Condensation of internal degrees of freedom for the quads
%
    if nen==4
      for n=1:6
	kk=30-n; 
	l=kk+1;
	pivot=s(l,l);
	if pivot > 0.0
	  for i=1:kk
	    gg=s(i,l);
	    if gg ~= 0.0
	      gg=gg/pivot;
	      s(i,l)=gg;
	      for j=i:kk;
		sq=s(i,j)-gg*s(l,j);
		s(i,j)=sq;
		s(j,i)=sq;
	      end
	    end 
	  end 
	end
      end 
      ek=s(1:24,1:24);
    else
      ek=s(1:18,1:18);
    end 
    em=plate_lumped(x,y,z,nen,em,rho,h,nsm);
    [k,m]=addkm(ek,em,id,k,m);
  end 
end 
ko=k;
mo=m;

function em=plate_lumped(x,y,z,nen,em,rho,h,nsm)
%PLATE_LUMPED em=plate_lumped(x,y,z,nen,em,rho,h,nsm)
%   Computes and assembles the lumped mass matrix for plate elements.
%   Normally is called by plate.
%
%   x,y,z are vectors storing the nodal coordinates
%   nen is the number of external nodes, 3 for triangles, 4 for quads
%   em is the element mass matrix
%   rho is the plate material density (per unit volume)
%   h is the plate thickness
%   nsm is the non-structural mass (per unit area)
%
%   This function only calculates the mass in the translational
%   degrees of freedom.  The rotational degrees of freedom must
%   use Guyan reduction unless other elements are attached which
%   provide values for the rotational degrees of freedom.
%
%   Written by Bob Norton, 10/6/92.

% Copyright (C) 1992, California Institute of Technology.
% U.S. Government Sponsorship under NASA Contract NAS7-918
% is acknowledged.
%
% Revision 00001

%History
%       31Aug93 jmelody - copied from Stuart Shaklan
%       31Aug93 jmelody - fixed apparent error with the cross.m
%                         function calls

if nen==3
%get the vectors from node 1 to 2 and 1 to 3
  v12=[x(2)-x(1);y(2)-y(1);z(2)-z(1)];
  v13=[x(3)-x(1);y(3)-y(1);z(3)-z(1)];
  v=cross(v12,v13);
  area=norm(v)/2;
  mass=area*(rho*h+nsm);
  mass3=mass/3;
  for i=1:6:13
    em(i,i)=mass3; em(i+1,i+1)=mass3; em(i+2,i+2)=mass3;
  end

else
%get the vectors from 1 to 2, 1 to 4, 2 to 3, and 3 to 4
  v14=[x(4)-x(1);y(4)-y(1);z(4)-z(1)];
  v12=[x(2)-x(1);y(2)-y(1);z(2)-z(1)];
  v23=[x(3)-x(2);y(3)-y(2);z(3)-z(2)];
  v34=[x(4)-x(3);y(4)-y(3);z(4)-z(3)];
  v1=cross(v14,v12);
  v2=cross(v12,v23);
  v3=cross(v23,v34);
  v4=cross(v34,v14);
%get the mass of each of the subtriangles
  m1=(rho*h+nsm)*norm(v1)/2;
  m2=(rho*h+nsm)*norm(v2)/2;
  m3=(rho*h+nsm)*norm(v3)/2;
  m4=(rho*h+nsm)*norm(v4)/2;
%sum for each node
  node=[(m1+m2+m4)/6;(m1+m2+m3)/6;(m2+m3+m4)/6;(m1+m3+m4)/6];
  for i=1:4
    j=6*(i-1)+1;
    em(j,j)=node(i);em(j+1,j+1)=node(i);em(j+2,j+2)=node(i);
  end
end

function [ko,mo]=addkm(ek,em,eid,k,m)
%addkm  [ko,mo]= addkm(ek,em,eid,k,m) assembles element k&m 
%				into system k&m
%
%	ek,em are the element stiffness and mass matrices
%	eid is a vector of degree of freedom numbers for the element.
%		Inactive degrees of freedom with eid(i) == 0 are allowed.
%	The system matrices k,m are returned in ko,mo.

[temp,ndof]=size(eid);
	for i=1:ndof
	ni=eid(1,i);
	if ni>0
		for j=1:ndof
		nj=eid(1,j);
		if nj> 0
			k(ni,nj)=k(ni,nj)+ek(i,j);
			if nargin > 4,
				m(ni,nj)=m(ni,nj)+em(i,j);
			end
		end
		end
	end
	end
ko=k;
if nargin > 4,
	mo=m;
end

function viewfem(xyz,ni,nn)
%VIEWFEM	plots a finite element geometry
%
%		viewfem(xyz,ni,nn)
%
%		xyz 	nodal coordinates
%		ni 	connectivities matrix
%		nn	number of nodes (see below)
%
%	The global variable labels is (1,2) and controls
%the number printing.  When using labels, "global labels"
%must be typed before a value is assigned to labels.
%
%		If labels(1,1) = 1 then node numbers are printed.
%		If labels(1,1) = 2 then the nodes are plotted.
%		If labels(1,2) = 1 then element numbers in ni(*,1)
%				   are printed near the mid point
%				   of the line segment.
%		If labels(1,2) = 2 then dashed lines are used.
%		If labels(1,2) = 3 then element property numbers
%				   are printed.
%
%Caution: the placement of the fixed size font labels can 
%	be confusing if many labels are printed.  The lower
%	left corner of the first character cell is the place
%	the labels are printed and the characters string out
%	to the right from there.  This is the middle of a
%	line segment or the node.
%
%	A third argument (nn) has been added to indicate the 
%number of nodes that define the outline of the element.
%If the argument is not used, a line segment with two end
%nodes will be assumed.  If there are more than two arguments
%the third will be used to define the number of nodes in a
%closed polygon.  In this case extra columns in ni will be
%interpretted as the node numbers.

% Copyright (C) 1992, California Institute of Technology.
% U.S. Government Sponsorship under NASA Contract NAS7-918
% is acknowledged.
%
% Revision 00001

%History
%  18Jun93 jmelody:   declaration of labels as a global variable
%                     is included as required by matlab4
%  29Sep92 jmelody:   dashed line capabilities added
%		      node plotting capabilities added
%  10Oct93 jmelody:   arbitrary node numbering capabilities added
%  19Oct93 jmelody:   played around with global definition of 
%		      labels, to try and make things better (not much)
%  25Jan94 jmelody:   add ability to plot element property numbers
%		      for labels(2)=3
%  17Jun94 jmelody:   removed element type from calls to nifix
%   7Nov97 gmosier:   changed line color from white to black,
%                     since v5 uses white background
%

global labels
if (exist('labels') ~= 1)  %if labels doesn't exist as a global variable
  labels=[0 0];
end

[nxyz,mxyz]=size(xyz);
if (mxyz == 3)		%if no node labels in xyz
  if nargin < 3
    [lx,ly]=list(xyz(:,1:2),ni);
  else
    [lx,ly]=list(xyz(:,1:2),ni,nn);
  end

  if labels(2) == 2
    plot(lx,ly,'--k')
  else
    plot(lx,ly,'-k')
  end

  if any(labels(1,1:2)) == 0
    return
  end

  [n,m]=size(ni);
  if nargin < 3
    nn = 2;
  end

  if labels(1) == 1
    hold on
    for i=1:n
      for j=2:nn+1
        text(xyz(ni(i,j),1),xyz(ni(i,j),2),int2str(ni(i,j)));
        plot(xyz(ni(i,j),1),xyz(ni(i,j),2),'+');
      end
    end
    hold off
  end

  if labels(1) == 2
    hold on
    for i=1:n
      for j=2:nn+1
        plot(xyz(ni(i,j),1),xyz(ni(i,j),2),'+');
      end
    end
    hold off
  end

  if labels(2) == 1	%element numbering
    for i=1:n
      cg=[0 0];
      for j=2:nn+1
        cg=cg+xyz(ni(i,j),1:2);
      end;
      cg=cg/nn;
      text(cg(1),cg(2),int2str(ni(i,1)));
    end
  end

  if labels(2) == 3	%property numbering
    for i=1:n
      cg=[0 0];
      for j=2:nn+1
        cg=cg+xyz(ni(i,j),1:2);
      end;
      cg=cg/nn;
      text(cg(1),cg(2),int2str(ni(i,m)));
    end
  end

elseif mxyz==4		%if arbitrary node numbering
  index=nodesort(xyz);

  if nargin < 3
    oldni=ni;
    ni=nifix(index,ni);
    [lx,ly]=list(xyz(:,2:3),ni);
  else
    oldni=ni;
    ni=nifix(index,ni);
    [lx,ly]=list(xyz(:,2:3),ni,nn);
  end

  if labels(2) == 2
    plot(lx,ly,'--k')
  else
    plot(lx,ly,'-k')
  end

  if any(labels(1,1:2)) == 0
    return
  end

  [n,m]=size(ni);
  if nargin < 3
    nn = 2;
  end

  if labels(1) == 1
    hold on
    for i=1:n
      for j=2:nn+1
        %this must be unfixed ni (oldni) to print the correct node numbers
        text(xyz(ni(i,j),2),xyz(ni(i,j),3),int2str(oldni(i,j)));
        plot(xyz(ni(i,j),2),xyz(ni(i,j),3),'+');
      end
    end
    hold off
  end

  if labels(1) == 2
    hold on
    for i=1:n
      for j=2:nn+1
        plot(xyz(ni(i,j),2),xyz(ni(i,j),3),'+');
      end
    end
    hold off
  end

  if labels(2) == 1	%element numbering
    for i=1:n
      cg=[0 0];
      for j=2:nn+1
        cg=cg+xyz(ni(i,j),2:3);
      end;
      cg=cg/nn;
      text(cg(1),cg(2),int2str(ni(i,1)));
    end
  end

  if labels(2) == 3	%property numbering
    for i=1:n
      cg=[0 0];
      for j=2:nn+1
        cg=cg+xyz(ni(i,j),2:3);
      end;
      cg=cg/nn;
      text(cg(1),cg(2),int2str(ni(i,m)));
    end
  end
else
  disp('nodal matrix is the wrong size')
end

function [lx,ly] = list(xy,ni,nn)
% LIST    [lx,ly] = list(xy,ni,nn)    creates ordered pairs for plotting 
%       with VIEWFEM
%
%   xy is a (:,2) array of x,y coordinate pairs in the image plane
%   ni is a (:,3) array of line segments where one end is xy(ni(*,2))
%       and the other end is xy(ni(*,3)).
%   nn is the number of sides in a polygon.  nn is optional

% Copyright (C) 1992, California Institute of Technology.
% U.S. Government Sponsorship under NASA Contract NAS7-918
% is acknowledged.
%
% Revision 00001

%History
%   8Sep95 jmelody:   make list.m capable of handling a 0 in the node list
%   7Nov97 gmosier:   initialized lx and ly as null matrices for v5 compatibility

lx = [];
ly = [];

[n,temp]=size(ni);
if nargin < 3
  for i=1:n
    lx=[lx,[xy(ni(i,2),1);xy(ni(i,3),1)]];
    ly=[ly,[xy(ni(i,2),2);xy(ni(i,3),2)]];
  end
else
  for i=1:n
    n=ni(i,[1:nn]+1);	%nodes
    n=n(find(n~= 0));	%throw away "zero" nodes
    for j=1:(length(n)-1)
      lx=[lx,[xy(n(j),1);xy(n(j+1),1)]];
      ly=[ly,[xy(n(j),2);xy(n(j+1),2)]];
    end
    lx=[lx,[xy(n(length(n)),1);xy(n(1),1)]];
    ly=[ly,[xy(n(length(n)),2);xy(n(1),2)]];
  end
end

function f=loads(F,bc);
% LOADS    f=loads(F,bc)    prepare load vectors from a table of forces
%       Multiple load cases are allowed.
%
%       Forces on inactive node degrees of freedom are allowed
%       but ignored.
%
%   F is a (:,4) table of force descriptions
%     F(:,1)is the node number and
%     F(:,2)is the degree of freedom number <=6
%     F(:,3)is the load case number
%     F(:,4)is the force magnitude
%   bc is the array of degree of freedom numbers
%   f is the load vectors sized (ndof,# load cases)

% Copyright (C) 1992, California Institute of Technology.
% U.S. Government Sponsorship under NASA Contract NAS7-918
% is acknowledged.
%
% Revision 00001

%History
% 13May96 briggs:   allow multiple loads on a node so that element
%		    by element load generators will work.  The change
%		    here is to sum the F's into f.


[nf,temp]=size(F);
ndof=max(max(bc));
nl=max(F(:,3));
f=zeros(ndof,nl);

%F (node, dir, loadcase, val)
for i=1:nf
  dof = bc(F(i,1),F(i,2));
  if (dof > 0)
    f(dof,F(i,3))=f(dof,F(i,3))+F(i,4);
  end
end

function dout=disps(din,bc,lc,tindex,tform);
% DISPS    dout=disps(din,bc,lc,tindex,tform)    prepares an array of
%       displacements for a given load case with only translations.
%
%   din is an array of displacement solutions sized (ndof,# load cases)
%   bc is the array of degree of freedom numbers
%   lc is the load case to construct displacements for.
%
%   The lack of multiple dimensioned arrays restricts this to 
%   extraction of a single load case.
%
%   dout will be sized (# of nodes,3) and filled with displacements.
%
%   A typical use for dout is to multiply it by a view factor and
%   add it to the original nodal coordinates to find the deformed
%   shape of the mesh.
%
%   xyzdeformed = xyz +blowup_factor*dout;
%
%   The function works for static solution vectors as well as dynamic
%   mode shapes in eigenvectors.
%
%   The remaining two arguments are optional and provide data for
%   node local coordinate systems.
%
%   tindex is #-of-nodes by 1 and contains the index for local
%   coordinate transformations.
%
%   tform is #-of-transforms by 9 and holds the 3x3 transformation
%   matrices stored by rows.
%

% Copyright (C) 1992, California Institute of Technology.
% U.S. Government Sponsorship under NASA Contract NAS7-918
% is acknowledged.
%
% Revision 00001

[ndof,nl]=size(din);
nodof=max(bc);
[nn,ndofn]=size(bc);
dout=zeros(nn,3);

for i=1:nn
  flag=0;
  if nargin > 3
    tni = tindex(i);
    tf=eye(3);
    if tni > 0
      tf(1,1:3)=tform(tni,1:3);
      tf(2,1:3)=tform(tni,4:6);
      tf(3,1:3)=tform(tni,7:9);
      flag=1;
    end
  end

  for j=1:3
    if nodof(j)>0,
      dof = bc(i,j);
      if dof > 0,
        dout(i,j)=din(dof,lc);
      end
    end
  end
  %apply the dof xform xposed so its local to global
  %dout(i,:)'=tf'*dout(i,:)';
  %temp=tf'*dout(i,:)';
  %dout(i,:)=temp';
  if flag==1
    dout(i,:)=dout(i,:)*tf;
  end
end
