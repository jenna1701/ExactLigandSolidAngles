(* ::Package:: *)


(* FindSolidAngle package for Mathematica, written by Wesley D. Allen*)
(* August 25, 2013 version *)

(* INPUT *)
(* XAtoms = array containing atomic symbols and Cartesian coordinates (in \[CapitalARing]) in the form {{"symbol1",{x1,y1,z1}},{"symbol2",{x2,y2,z2}},...};  If "symbol" is a number rather than an atomic symbol, then that number is used for the corresponding atomic radius instead of a standard value;
Apex = integer specifying which is the apex atom in XAtoms;
Ligands = integer array in the form {{i1,i2,i3,...},{j1,j2,j3,...},...} specifying which atoms in XAtoms correspond to each ligand for which a cone angle is to be computed;
kPrint = output option as specified below;
kR != 2 to use internal van der Waals atomic radii; kR = 2 to use internal zero energy point atomic radii *)

(* 
kPrint = 0  No printing within package;
kPrint = 1  Print total solid angle and dissection of loop contributions for each ligand; 
kPrint = 2  For each ligand print a 3D plot highlighting the border arcs and showing the shadow cones intersecting with the unit sphere;
kPrint = 3  For each ligand print a table of van der Waals radii, vertex angles, and Cartesian coordinates;
kPrint = 4  For each ligand print the aforementioned 3D plot with labels on the arcs and for comparison also tabulate the dissection of loop contributions for each ligand;
kPrint = 5  For each ligand print a 3D plot showing the border arcs and filling gray into the shadow cones intersecting with the unit sphere; *)    

(* OUTPUT *)

(*RvdW = Default van der Waals radii (in \[CapitalARing]) from Bondi (1964);
RZ = zero energy point atomic radii (in \[CapitalARing]) from Guzei and Wendt (2006);
NThreshold = numerical threshold for determining if an atom is contained in a cone;
nLigand = number of atoms in ligand k;
XLigand = array to hold Cartesian coordinates of ligand k;
RvdWLigand = atomic radii of atoms in ligand k;
XApex = Cartesian coordinates of apex atom;
rL = distances from apex to ligand atoms;
mL = unit vectors from apex to ligand atoms;
\[Beta]L = vertex angles (rad) of ligand atoms;  *)

FindSolidAngle[XAtoms_,Apex_,Ligands_,kPrint_,kR_]:=Module[{ConeAtoms,ConeAxis,ConeAngle,RvdW,RZ,Ratoms,nLigand,XLigand,RvdWLigand,XApex,kvdW,rL,mL,\[Beta]L,\[DoubledGamma]L,\[CapitalPhi],\[Theta],\[Phi],\[Chi],x1,x2,x3,y1,y2,y3,s1,s2,nvec,r,CrossXY,Cjk,F12,F13,gsphere,D2jk,NThreshold,Ligandtable,tconeL,tpos,Cone1,\[Phi]cone,Cones,Lc,Arcs,BorderArcs,m,p,q,i,j,k,kL,\[Theta]k,\[Phi]k,\[Alpha]k,\[Alpha]jk,nk,\[Theta]j,\[Phi]j,\[Alpha]j,nj,ax,bx,cx,dx,\[Phi]a,\[Phi]b,Xborder,disjoint,\[Phi]arcs,\[Phi]arcs0,Xmidarc,Arctest,\[Phi]intervals,turnangles,CosCones,\[CapitalOmega]disjoint,\[CapitalOmega]loops,ksort,Nloops,nloop,nborder,nArcs,nconvex,molscale,\[CapitalOmega]total,\[CapitalOmega],\[CapitalOmega]angle,\[CapitalOmega]data,Xloop,Xplot,arcstyle,\[Phi]plot,ConeArcPlot,ConeSurfaces,ConeSurfacePlot,ShadowPlot,BorderPointPlot,SurfaceFill,SolidAngleResults,PlotSpecs,ThickLine,ThinLine,AtomColor},
RvdW={{"H","He","C","N","O","F","Ne","Si","P","S","Cl","Ar","As","Se","Br","Kr","Te","I","Xe"},{1.20,1.40,1.70,1.55,1.52,1.47,1.54,2.10,1.80,1.80,1.75,1.88,1.85,1.90,1.85,2.02,2.06,1.98,2.16},
{White,White,Gray,Blue,Red,Yellow,White,Gray,Orange,Yellow,Green,White,White,White,White,White,White,White,White}};
RZ={{"H","He","C","N","O","F","Ne","Si","P","S","Cl","Ar","As","Se","Br","Kr","Te","I","Xe"},{1.000,1.311,1.539,1.521,1.470,1.413,1.350,1.834,1.801,1.757,1.599,1.649,1.879,1.861,1.845,1.831,1.955,1.941,1.928},
{White,White,Gray,Blue,Red,Yellow,White,Gray,Orange,Yellow,Green,White,White,White,White,White,White,White,White}};
PlotSpecs={PlotRange-> All,ImageSize-> 500,BoxRatios->{1,1,1},Boxed-> False,Axes-> False,Mesh-> False,PlotPoints-> Automatic};
ThickLine=Thickness[0.004];ThinLine=Thickness[0.001];
Clear[CrossXY,Cjk,D2jk,F12,F13,\[Phi],x1,x2,x3,y1,y2,y3,\[CapitalPhi],nvec,gsphere,\[Theta],\[Phi],\[Chi]];
(* Use this explicit function instead of inefficient Mathematica Cross *)
CrossXY[{x1_,x2_,x3_},{y1_,y2_,y3_}]={x2 y3-x3 y2,x3 y1-x1 y3,x1 y2-x2 y1};
(* Direction cosine matrix *)
\[CapitalPhi][\[Theta]_,\[Phi]_,\[Chi]_]={{Cos[\[Theta]] Cos[\[Phi]] Cos[\[Chi]]-Sin[\[Phi]] Sin[\[Chi]],Cos[\[Theta]] Sin[\[Phi]]Cos[\[Chi]]+Cos[\[Phi]]Sin[\[Chi]],-Sin[\[Theta]]Cos[\[Chi]]},{-Cos[\[Theta]] Cos[\[Phi]]Sin[\[Chi]]-Sin[\[Phi]] Cos[\[Chi]],-Cos[\[Theta]] Sin[\[Phi]]Sin[\[Chi]]+Cos[\[Phi]] Cos[\[Chi]],Sin[\[Theta]]Sin[\[Chi]]},{Sin[\[Theta]] Cos[\[Phi]],Sin[\[Theta]] Sin[\[Phi]],Cos[\[Theta]]}};
(* Unit vector along (\[Theta],\[Phi]) *)
nvec[\[Theta]_,\[Phi]_]={Sin[\[Theta]] Cos[\[Phi]],Sin[\[Theta]] Sin[\[Phi]],Cos[\[Theta]]};
Cjk[\[Theta]_,\[Phi]_,\[Chi]_]=Cos[\[Phi]]-Cos[\[Theta]] Cos[\[Chi]];
D2jk[\[Theta]_,\[Phi]_,\[Chi]_]=1-(Cos[\[Chi]])^2-(Cos[\[Theta]])^2-(Cos[\[Phi]])^2+2Cos[\[Theta]]Cos[\[Phi]]Cos[\[Chi]];
F12[\[Theta]_,\[Phi]_,\[Chi]_]=Sin[\[Theta]] Cos[\[Phi]] Cos[\[Chi]]-Cos[\[Theta]]Sin[\[Phi]];
F13[\[Theta]_,\[Phi]_]=Sin[\[Theta]]  Sin[\[Phi]];
gsphere[r_,x1_,y1_,x2_,y2_]=(r^2-(x1-x2)^2-(y1-y2)^2)^(1/2);
NThreshold=10^-10;
XApex=XAtoms[[Apex,2]];
SolidAngleResults={};
(* Loop over all ligands *)
Do[
nLigand=Length[Ligands[[kL]]];
XLigand=Table[XAtoms[[Ligands[[kL,j]],2]]-XApex,{j,nLigand}]//N;
(* Collect radii of ligand atoms *)
Ratoms=RvdW;
If[kR== 2,Ratoms=RZ];
kvdW=Table[Flatten[Position[Ratoms,XAtoms[[Ligands[[kL,j]],1]]]],{j,nLigand}];
RvdWLigand=Table[XAtoms[[Ligands[[kL,j]],1]],{j,nLigand}];
Do[If[Length[kvdW[[j]]]!= 0,If[kvdW[[j,1]]==1,RvdWLigand[[j]]=Ratoms[[2,kvdW[[j,2]]]]]],{j,Length[kvdW]}];
(* Vectors and distances from apex to ligand atoms *)
rL=Table[Norm[XLigand[[j]]],{j,nLigand}];
mL=Table[XLigand[[j]]/rL[[j]],{j,nLigand}];
(* Vertex angles (rad) of ligand atoms *)
\[Beta]L=Table[ArcSin[RvdWLigand[[j]]/rL[[j]]],{j,nLigand}];
If[(kPrint== 3),Print["Ligand atom van der Waals radii (r), vertex angles (\[Beta], deg), and Cartesian coordinates (x,y,z)"];
Ligandtable=Table[Flatten[{XAtoms[[Ligands[[kL,j]],1]],Ligands[[kL,j]],RvdWLigand[[j]],\[Beta]L[[j]]/Degree,XAtoms[[Ligands[[kL,j]],2]]}],{j,nLigand}];
Print[TableForm[Join[{{"Atom","Number ","  r","  \[Beta]","  x","  y","  z"}},Ligandtable]]]];
(* Angles (rad) between direction vectors to ligand atoms *)
\[DoubledGamma]L=Table[0,{i1,nLigand},{i2,nLigand}];
Do[\[DoubledGamma]L[[i1,i2]]=ArcCos[Dot[mL[[i1]],mL[[i2]]]];\[DoubledGamma]L[[i2,i1]]=\[DoubledGamma]L[[i1,i2]],{i1,2,nLigand},{i2,i1-1}];
(* tconeL(i,j) test of whether atom i lies completely within cone of atom j *)
tconeL=Table[\[Beta]L[[i2]]>=\[Beta]L[[i1]]+\[DoubledGamma]L[[i1,i2]] ,{i1,nLigand},{i2,nLigand}];
Do[tconeL[[i1,i1]]=False,{i1,nLigand}];
tpos=Position[tconeL,True];
(* Ligand atoms that lie in shadow of another ligand atom *)
q=Union[Table[tpos[[i1,1]],{i1,Length[tpos]}]];
(* Cone1 has active atoms for determining solid angle; Lc = number of active atoms *)
Cone1=Complement[Range[nLigand],q];
Lc=Length[Cone1];
\[Phi]cone=Table[{ax,bx,cx}=mL[[Cone1[[j]]]];
If[Abs[ax]>= NThreshold,\[Phi]=ArcTan[bx/ax],\[Phi]=(\[Pi]/2)Sign[bx]];
If[ax<0,\[Phi]=\[Pi]+\[Phi]];\[Phi],{j,Lc}];
(* Cones contains (\[Theta],\[Phi]) angles of cone axes and \[Alpha] cone angle for active atoms *)
Cones=Table[{ArcCos[mL[[Cone1[[j]],3]]],\[Phi]cone[[j]],\[Beta]L[[Cone1[[j]]]]},{j,Lc}];
(* \[Phi]arcs[[k]] contains a list of {\[Phi],j} values for the local azimuthal angle of all intersections of cone k with cone j *)
\[Phi]arcs={};
Do[{\[Theta]k,\[Phi]k,\[Alpha]k}=Cones[[k]];\[Phi]arcs0={};
Do[{\[Theta]j,\[Phi]j,\[Alpha]j}=Cones[[j]];\[Alpha]jk=Re[ArcCos[Dot[nvec[\[Theta]j,\[Phi]j],nvec[\[Theta]k,\[Phi]k]]]];
If[(D2jk[\[Alpha]j,\[Alpha]k,\[Alpha]jk]>= 0)&&(j!= k),
{ax,cx}=Csc[\[Alpha]k](Csc[\[Alpha]jk])^2 Cjk[\[Alpha]k,\[Alpha]j,\[Alpha]jk]{F12[\[Theta]j,\[Theta]k,\[Phi]j-\[Phi]k],F13[\[Theta]j,\[Phi]j-\[Phi]k]};
{bx,dx}=Csc[\[Alpha]k](Csc[\[Alpha]jk])^2 Sqrt[D2jk[\[Alpha]j,\[Alpha]k,\[Alpha]jk]]{F13[\[Theta]j,\[Phi]j-\[Phi]k],-F12[\[Theta]j,\[Theta]k,\[Phi]j-\[Phi]k]};
If[cx+dx>= 0,\[Phi]a=ArcCos[ax+bx],\[Phi]a=2\[Pi]-ArcCos[ax+bx]];
If[cx-dx>= 0,\[Phi]b=ArcCos[ax-bx],\[Phi]b=2\[Pi]-ArcCos[ax-bx]];
\[Phi]arcs0=Join[\[Phi]arcs0,{{\[Phi]a,j},{\[Phi]b,j}}]],{j,Lc}];
\[Phi]arcs=Join[\[Phi]arcs,{\[Phi]arcs0}],{k,Lc}];
(* nArcs contains number of intersection points for each cone k *)
nArcs=Table[Length[\[Phi]arcs[[k]]],{k,Lc}];
(* Compute solid angle contributions for all disjoint cones that have no intersection points. *)
disjoint={};\[CapitalOmega]disjoint={};
Do[If[nArcs[[k]]==0,disjoint=Join[disjoint,{k}];
{\[Theta]k,\[Phi]k,\[Alpha]k}=Cones[[k]];ax=2\[Pi] (1-Cos[\[Alpha]k]);
\[CapitalOmega]disjoint=Join[\[CapitalOmega]disjoint,{ax}]],{k,Lc}];
\[Phi]arcs0=Table[Sort[Table[\[Phi]arcs[[k,j,1]],{j,nArcs[[k]]}]],{k,Lc}];
ksort=Table[Take[Flatten[Position[\[Phi]arcs,\[Phi]arcs0[[k,j]]]],2][[2]],{k,Lc},{j,nArcs[[k]]}];
(* Sort each \[Phi]arcs list by increasing \[Phi] value to define arc sequence for loops *)
\[Phi]arcs=Table[\[Phi]arcs[[k,ksort[[k,j]]]],{k,Lc},{j,Length[ksort[[k]]]}];
(* Compute midpoints of all arcs of intersection *)
Xmidarc=Table[{\[Theta]k,\[Phi]k,\[Alpha]k}=Cones[[k]];
Table[If[j<nArcs[[k]],\[Phi]=(\[Phi]arcs[[k,j,1]]+\[Phi]arcs[[k,j+1,1]])/2,\[Phi]=\[Pi]+(\[Phi]arcs[[k,j,1]]+\[Phi]arcs[[k,1,1]])/2];
Dot[Transpose[\[CapitalPhi][\[Theta]k,\[Phi]k,0]],nvec[\[Alpha]k,\[Phi]]],{j,nArcs[[k]]}],{k,Lc}];
(* Test Xmidarc points to determine if a given arc is a border arc or inside a loop *)
Arctest=Table[ax=True;k=0;
While[ax&&(k<Lc),k=k+1;{\[Theta]k,\[Phi]k,\[Alpha]k}=Cones[[k]];
ax=Sign[Chop[ArcCos[Dot[nvec[\[Theta]k,\[Phi]k],Xmidarc[[m,j]]]]-\[Alpha]k,NThreshold]]>= 0];ax,{m,Lc},{j,nArcs[[m]]}];
(* Collect unsequenced border arc indices {m,k} corresponding to arc k of atom m in \[Phi]arcs. *)
BorderArcs={};
Do[If[Arctest[[m,k]],BorderArcs=Join[BorderArcs,{{m,k}}]],{m,Lc},{k,nArcs[[m]]}];
nborder=Length[BorderArcs];
(* Xborder contains Cartesian coordinates of link points of border arcs *)
Xborder=Table[
{k,j}=BorderArcs[[m]];\[Phi]a=\[Phi]arcs[[k,j,1]];{\[Theta]k,\[Phi]k,\[Alpha]k}=Cones[[k]];
If[j<nArcs[[k]],\[Phi]b=\[Phi]arcs[[k,j+1,1]],\[Phi]b=\[Phi]arcs[[k,1,1]]];
{Dot[Transpose[\[CapitalPhi][\[Theta]k,\[Phi]k,0]],{Cos[\[Phi]a]Sin[\[Alpha]k],Sin[\[Phi]a]Sin[\[Alpha]k],Cos[\[Alpha]k]}],
Dot[Transpose[\[CapitalPhi][\[Theta]k,\[Phi]k,0]],{Cos[\[Phi]b]Sin[\[Alpha]k],Sin[\[Phi]b]Sin[\[Alpha]k],Cos[\[Alpha]k]}]},{m,nborder}];
(* Use Xborder values to sequence the border arcs in one or more closed loops *)
ksort={};m=0;q=1;
While[m<nborder,k=q;p=0;ax={};
While[p!= q,m=m+1;bx=Table[Chop[Norm[Xborder[[k,2]]-Xborder[[j,1]]],NThreshold],{j,Length[Xborder]}];
p=Position[bx,0][[1,1]];k=p;ax=Join[ax,{p}]];
ksort=Join[ksort,{ax}];
If[m<nborder,q=Complement[Range[nborder],Flatten[ksort]][[1]]]];
(* Nloops = number of non-disjoint loops; nloop[[p]] = number of arcs in loop p *)
Nloops=Length[ksort];
nloop=Table[Length[ksort[[p]]],{p,Nloops}];
ksort=Table[ksort[[p,nloop[[p]]-Mod[nloop[[p]]+1-m,nloop[[p]]]]],{p,Nloops},{m,nloop[[p]]}];
(* Sort BorderArcs so that (p,m) element contains indices for arc m of loop p *)
BorderArcs=Table[BorderArcs[[ksort[[p,m]]]],{p,Nloops},{m,nloop[[p]]}];
(* \[Phi]intervals[[p,m]] has change in local azimuthal angles for border arc m of loop p. *)
\[Phi]intervals=
Table[{k,j}=BorderArcs[[p,m]];
If[j<Length[\[Phi]arcs[[k]]],\[Phi]=\[Phi]arcs[[k,j+1,1]]-\[Phi]arcs[[k,j,1]],\[Phi]=2\[Pi]+\[Phi]arcs[[k,1,1]]-\[Phi]arcs[[k,j,1]]];\[Phi],{p,Nloops},{m,nloop[[p]]}];
(* CosCones[[p,m]] cosine of cone angle corresponding to border arc m of loop p. *)
CosCones=Table[Cos[Cones[[BorderArcs[[p,m,1]],3]]]//N,{p,Nloops},{m,nloop[[p]]}];
(* turnangles[[p,m]] has turn angle for the end of arc m of loop p. *)
turnangles=Table[{k,j}={BorderArcs[[p,m,1]],BorderArcs[[p,1+Mod[m,nloop[[p]]],1]]};
{\[Theta]k,\[Phi]k,\[Alpha]k}=Cones[[k]];nk=nvec[\[Theta]k,\[Phi]k];{\[Theta]j,\[Phi]j,\[Alpha]j}=Cones[[j]];nj=nvec[\[Theta]j,\[Phi]j];
ArcCos[Csc[\[Alpha]j]Csc[\[Alpha]k](Dot[nk,nj]-Cos[\[Alpha]j]Cos[\[Alpha]k])],{p,Nloops},{m,nloop[[p]]}]; 
(* Add the solid angle contributions for non-disjoint loops to existing sum for disjoint loops *)
\[CapitalOmega]loops={};
If[Nloops>0,\[CapitalOmega]loops=Table[2\[Pi]+Sum[turnangles[[p,m]]-CosCones[[p,m]] \[Phi]intervals[[p,m]],{m,nloop[[p]]}],{p,Nloops}]];
\[CapitalOmega]=Join[\[CapitalOmega]disjoint,\[CapitalOmega]loops];\[CapitalOmega]total=Tr[\[CapitalOmega]];
nconvex=Floor[Abs[\[CapitalOmega]total]/(4 \[Pi])];
If[\[CapitalOmega]total>4\[Pi],\[CapitalOmega]total=\[CapitalOmega]total-4 \[Pi] nconvex];
\[CapitalOmega]angle=2ArcCos[1-\[CapitalOmega]total/(2\[Pi])]/Degree;
SolidAngleResults=Join[SolidAngleResults,{\[CapitalOmega]total,\[CapitalOmega]angle}];
(* Report results of solid angle computations *)
If[kPrint!= 0,
Print["Number of loops whose concave domain is inside a convex domain = ",nconvex];
Print["Number of disjoint loops = ",Length[disjoint]];
Print["\[CapitalOmega] convex to loops = ",\[CapitalOmega],"   \[CapitalOmega] concave to loops = ",4\[Pi]-\[CapitalOmega]];
Print["Net \[CapitalOmega] of ligand shadow = ",\[CapitalOmega]total,"  Solid cone angle (deg) = ",\[CapitalOmega]angle]];
If[(kPrint== 1)||(kPrint== 4),
\[CapitalOmega]data= Table[j=Cone1[[BorderArcs[[p,m,1]]]];
{XAtoms[[Ligands[[kL,j]],1]],Ligands[[kL,j]],\[Phi]intervals[[p,m]]/Degree,turnangles[[p,m]]/Degree},{p,Nloops},{m,nloop[[p]]}];
Do[Print["Ligand = ",kL,"   Loop = ",p];Print[TableForm[Join[{{"Cone Atom","Number","Arc \[Phi] (deg)","Turn Angle (deg)"}},\[CapitalOmega]data[[p]]]]],{p,Nloops}];
Do[Print["Disjoint shadow cone formed by atom ",Ligands[[kL,disjoint[[k]]]],"."],{k,Length[disjoint]}]];
(* Plot section *)
If[(kPrint== 2)||(kPrint== 4)||(kPrint== 5),
If[kPrint== 4,BorderPointPlot=True,BorderPointPlot=False];
If[kPrint== 5,SurfaceFill=True,SurfaceFill=False];
Clear[\[Phi],\[Theta]];Xplot={};
ConeSurfaces=Table[{\[Theta]k,\[Phi]k,\[Alpha]k}=Cones[[k]];nk=nvec[\[Theta]k,\[Phi]k];Dot[Transpose[\[CapitalPhi][\[Theta]k,\[Phi]k,0]],nvec[\[Theta],\[Phi]]],{k,Lc}];
Arcs=Table[Evaluate[ConeSurfaces[[k]]]/.{\[Theta]->Cones[[k,3]]},{k,Lc}];
If[Lc== Length[disjoint],BorderPointPlot=False];
If[BorderPointPlot,
Xloop=Flatten[Table[{k,j}=BorderArcs[[p,m]];{Ligands[[kL,Cone1[[k]]]],Xmidarc[[k,j]]},{p,Nloops},{m,nloop[[p]]}],1];
Xplot=Graphics3D[Join[Table[{Text[Style[ToString[Xloop[[m,1]]],FontSize-> 14],Xloop[[m,2]],Background->White]},{m,Length[Xloop]}]]]];
If[Length[BorderArcs]> 0,
arcstyle=Table[{ThinLine},{p,Lc}];
Do[If[nArcs[[p]]>0,
bx=Table[If[Arctest[[p,k]],ax={ThickLine},ax={ThinLine}];ax,{k,nArcs[[p]]}];arcstyle[[p]]=bx],{p,Lc}];
If[Length[disjoint]>0,
Do[arcstyle[[disjoint[[m]]]]={ThickLine},{m,Length[disjoint]}]];
\[Phi]plot=Table[{0,2\[Pi]},{p,Lc},{k,Max[nArcs[[p]],1]}];
Do[
If[nArcs[[p]]>0,
Do[\[Phi]plot[[p,k]]={\[Phi]arcs[[p,k,1]],\[Phi]arcs[[p,k+1,1]]},{k,nArcs[[p]]-1}];
\[Phi]plot[[p,-1]]={\[Phi]arcs[[p,-1,1]],\[Phi]arcs[[p,1,1]]+2 \[Pi]}],{p,Lc}];
ConeArcPlot=Flatten[Table[
ParametricPlot3D@@Join[{Arcs[[p]],{\[Phi],\[Phi]plot[[p,k,1]],\[Phi]plot[[p,k,2]]},PlotStyle-> arcstyle[[p,k]]},PlotSpecs],{p,Lc},{k,Length[\[Phi]plot[[p]]]}],1];

(*ProjArcPlot=Flatten[Table[
ParametricPlot@@Join[{ArcCos[Arcs[[p,3]]]/Sqrt[1-(Arcs[[p,3]])^2] {Arcs[[p,1]],Arcs[[p,2]]},{\[Phi],\[Phi]plot[[p,k,1]],\[Phi]plot[[p,k,2]]},PlotRange-> All,Frame-> True,Axes-> False}],{p,Lc},{k,Length[\[Phi]plot[[p]]]}],1],*)

ConeArcPlot=Table[ParametricPlot3D@@Join[{Arcs[[p]],{\[Phi],0,2\[Pi]},PlotStyle->ThickLine},PlotSpecs]
,{p,Lc}]];
If[SurfaceFill,ax=0.98;
ConeSurfacePlot=Flatten[Table[
ParametricPlot3D@@Join[{ConeSurfaces[[p]],{\[Theta],0,Cones[[p,3]]},{\[Phi],0,2\[Pi]},PlotStyle-> {Gray,Lighting->"Neutral",Opacity[1.0]}},PlotSpecs],{p,Lc}],1],ax=1.0;ConeSurfacePlot={}];

AtomColor=Table[Transparent,{j,nLigand}];
Do[If[Length[kvdW[[j]]]!= 0,If[kvdW[[j,1]]==1,AtomColor[[j]]=RvdW[[3,kvdW[[j,2]]]]]],{j,nLigand}];
(*ProjPlot=Show[ProjArcPlot];*)
ShadowPlot=Show@@Join[ConeArcPlot,ConeSurfacePlot,{Graphics3D[{Lighting->"Neutral",Yellow,Opacity[1.0],Sphere[{0,0,0},ax]}],Xplot}];
Print[ShadowPlot](*;Print[ProjPlot]*)],{kL,Length[Ligands]}];SolidAngleResults];


(* ::Input:: *)
(**)
(**)
