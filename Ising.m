(* ::Package:: *)

(*           Two dimensional simple ising model simulation 1993.8.4                             *)
(*                  NVT-const : one-spin flop Monte Carlo simulation                            *)

(*                  K. Nishidate, K. Nishikawa(i),                              *)
(*                  M. Baba, T. Ikeda                                           *)
(* Dept. of Electronic Engineering, Faculty of Engineering, Iwate University, Morioka 020       *)
(* (i) Dept. of Chemistory, Faculty of Science,Kanazawa University, Kanazawa 920                *)
(*   References :
      # D. W. Heermann : Computer Simulation Methods in Theoretical Physics
                         (2nd. edn) , Springer-Verlag Berlin Heidelbelg 1986 and 1990 
      # S. E. Koonin: Computational Physics, Benjamin Cummings, 1986                            *)

(*-----------------<<<<<<<<<<((((((     Ising spin system.     ))))))>>>>>>>>>>-----------------*)

(* This Mathematica packege was debelopped at HP-730 WS *)
(* initial settiong *)
Maxstep=5000;          (* maximum step *)
Gstep=500;             (* graphics output *)
LOstep=10;             (* LongOrder count step *)
L=50;                  (* length *)
L2=L L;                (* total number *)
Bkt= 0;                (* = external magnetic field /kT *)
f[x_]:= x (-1);        (* spin flop opereter *)
Jkt= 0.25;             (* = J/kT *) 

(* random spin system *)
SeedRandom[];
ml=Table[Random[Integer],{L},{L}];
For[j=1,j<=Length[ml],j++,
   For[i=1,i<=Length[ml],i++,
      If[ml[[j,i]]==0,ml[[j,i]]=-1 
      ]
   ]
];
LongOrder=Abs[Apply[Plus,Flatten[ml]]/L2];
gra=DensityGraphics[ListDensityPlot[ml,Mesh->False]];

(* Metropolis methods with cyclic boudary condition *)
Do[
   i1 = Random[Integer,{1,L}]; i2=Random[Integer,{1,L}];
      j1 = i1+1;j2=i2+1; k1=i1-1;k2=i2-1;
          If[j1>L,j1=1]; If[j2>L,j2=1];
          If[k1<1,k1=L]; If[k2<1,k2=L];
   dG = ml[[j1,i2]]+ml[[k1,i2]]+ml[[i1,j2]]+ml[[i1,k2]];
   dE =-2f[ml[[i1,i2]]] ( Jkt dG + Bkt );
   dW = N[Exp[ -dE ]];
   If[dW<1 && dW>Random[] || dE<0,
                            ml=MapAt[f,ml,{{i1,i2}}]]; 
   If[Mod[i,LOstep]==0,
          LongOrder=Flatten[{LongOrder,
              Abs[Apply[Plus,Flatten[ml]]]/L2}]
   ];
   If[Mod[i,Gstep]==0,
      Show[DensityGraphics[
                 ListDensityPlot[ml,Mesh->False]]] ],
{i,1,Maxstep}]
 ListPlot[LongOrder,PlotLabel->"Long-range order parameter"];


(* Professor Richard J. Gaylord of the Department of Materials Science at the
   University of Illinois pointed out a correction to the program.
   Please enjoy it.

       nisidate@hep.s.kanazawa-u.ac.jp
       kiyoshi@hep.s.kanazawa-u.ac.jp
                                                                                *)

