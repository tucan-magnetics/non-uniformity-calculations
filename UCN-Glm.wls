#!/usr/bin/env wolframscript
(* ::Package:: *)

c[l_, m_, phi_] = (l-1)!(-2)^Abs[m]/(l+Abs[m])!If[m>=0, Cos[Abs[m]*phi], Sin[Abs[m]*phi]];


sigma[l_, m_, r_, phi_, theta_] = c[l, m, phi]*r^l*LegendreP[l, Abs[m], Cos[theta]];


sigmac[l_, m_, x_, y_, z_]=TransformedField["Spherical"->"Cartesian", sigma[l, m, r, phi, theta], {r, theta, phi}->{x, y, z}];


pi[l_, m_, x_, y_, z_]={D[sigmac[l+1, m, x, y, z], x], D[sigmac[l+1, m, x, y, z], y], D[sigmac[l+1, m, x, y, z], z]};


maxL=2
glm=Table[Table[Symbol["$g"<>ToString@(i-1)<>StringReplace[ToString@(j-i-1), "-"->"n"]],{j,2i+1}],{i,maxL+1}]
g[l_, m_]=Indexed[glm, {l+1,m+l+2}]
g[l_]=Indexed[glm,l+1]
(*For[i=0, i \[LessEqual]maxL, i++, For[j = -i-1, j \[LessEqual]i+1, j++, SetAttributes[g[i, j], Constant]]*)
(*Map[SetAttributes[#, Constant]&, glm, {2}]*)


polyExpand=Simplify[TrigExpand[{bx[x, y, z], by[x, y, z], bz[x, y, z]}==Sum[Sum[g[l,m]*pi[l, m, x, y, z], {m, -l-1, l+1}], {l, 0, maxL}]], x^2+y^2+z^2>0]
(*Simplify[polyExpand, {x\[Equal]0, y\[Equal]0, z\[Equal]0}]*)


equations = {polyExpand};
Do[
    Print[Flatten[Normal[Solve[Simplify[equations[[i+1]], {x==0, y==0, z==0}], g[i], MaxExtraConditions->Automatic]]]];
    AppendTo[equations,Flatten[{D[equations[[i+1]], x], D[equations[[i+1]], y], D[equations[[i+1]], z]}]]
    (*Print[i];
    Print[equations[[i+1]]];*)
, {i, 0, maxL}]
