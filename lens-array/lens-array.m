(* ::Package:: *)

(*p=None;
size=None*p;
\[Lambda]=None;
\[Theta]x=None Degree;
\[Theta]y=None Degree;
f=None;
padRatio=None;
T2=None;
*)



rotL[data_]:=RotateLeft[data, Floor@Dimensions[data][[#]]/2&/@Range[Length@Dimensions[data]]];
(* corner at orign to center at origin *)
rotR[data_]:=RotateRight[data, Floor@Dimensions[data][[#]]/2&/@Range[Length@Dimensions[data]]];
(* {x,y} means {x,y}/\[Lambda], {kx,ky} means {kx,ky}/k *)
FFT[f_] :=rotR[Fourier[f, FourierParameters->{1, -1}]];
iFFT[f_]:=InverseFourier[rotL[f], FourierParameters->{1, -1}];

centerRange[n_,d_]:=Table[(i-n/2-.5)*d,{i,n}];

(*plot[plotFunction_,func_, args_]:=Map[func/*plotFunction, args, {2}]//GraphicsGrid//Rasterize;*)
Options[plot]={plotFunction->ArrayPlot, cropRatio->1, rasterize->True};
plot[func_, args_, OptionsPattern[]]:=
Module[{},
	Map[func/*(ArrayPad[#,-Floor[(1-OptionValue[cropRatio])/2*Dimensions@#]]&)/*OptionValue[plotFunction], 
	args, {2}]//GraphicsGrid
	]


k0=2\[Pi]/\[Lambda]
{nx,ny}={Floor[size/p],Floor[size/p]};
{lx,ly}={nx*p,ny*p};
{dx,dy}={p,p};
{dkx, dky}={(2\[Pi])/(p*nx),(2\[Pi])/(p*ny)};
{padx,pady}= {(padRatio-1)/2 nx,(padRatio-1)/2 ny};
k = k0*Sin/@({
 {{-\[Theta]x, -\[Theta]y}, {0, -\[Theta]y}, {\[Theta]x, -\[Theta]y}},
 {{-\[Theta]x, 0}, {0, 0}, {\[Theta]x, 0}},
 {{-\[Theta]x, \[Theta]y}, {0, \[Theta]y}, {\[Theta]x, \[Theta]y}}
})//N;
r0 = -p*({
 {{-nx, -ny}, {0, -ny}, {nx, -ny}},
 {{-nx, 0}, {0, 0}, {nx, 0}},
 {{-nx, ny}, {0, ny}, {nx, ny}}
})//N;
{kx,ky}={Flatten@Part[k, ;;,;;,1],Flatten@Part[k, ;;,;;,2]};
{x0,y0}={Flatten@Part[r0, ;;,;;,1],Flatten@Part[r0, ;;,;;,2]};
{fx,fy}={(f*kx)/Sqrt[k0^2-kx^2],(f*ky)/Sqrt[k0^2-ky^2]};
{x,y}={centerRange[nx,dx],centerRange[ny,dy]};
Echo[{nx p,ny p},"Sample size in um: "];
Echo[{nx,ny},"Sample size in pixel: "];


ellipticalMask:=Map[If[#<(lx ly)/4,1,0]&,Outer[Plus,y^2,x^2],{2}];
globalEllipticalMask := KroneckerProduct[({
 {1, 1, 1},
 {1, 1, 1},
 {1, 1, 1}
}), ellipticalMask];
eGlobalLens[arg_]:= ArrayFlatten[Map[eLocalLens,arg,{2}]];
eGlobalInput[arg_]:=ArrayFlatten[Map[#*ConstantArray[1,{ny,nx}]&,arg,{2}]];
H[A_,z_,{dx_,dy_}] := Module[{nx, ny, kx, ky},
	{nx, ny}=Dimensions[A][[#]]&/@{2,1};
	{kx,ky}=centerRange@@#&/@{{nx, (2\[Pi])/(dx nx)},{ny,(2\[Pi])/(dy ny)}};
	Exp[I *z*Sqrt[k0^2 - Outer[Plus,ky^2,kx^2]]]
]; 
forwardPropagate[A_,z_,{dx_, dy_}]:=iFFT[FFT[A] *H[A,z,{dx, dy}]];


convertPair[meta_] := Module[{\[Phi]1, \[Phi]2},
	\[Phi]1 = Arg[#]+ArcCos[Abs[#]]&;
	\[Phi]2 = Arg[#]-ArcCos[Abs[#]]&;
	Exp[I Outer[#2[#1]&,meta/Max[Abs[meta]],({
 {\[Phi]1, \[Phi]2},
 {\[Phi]2, \[Phi]1}
})]]//ArrayFlatten
];

convertArg[meta_] := Module[{\[Phi]0},
	\[Phi]0 = Arg[#]&;
	Exp[I Outer[#2[#1]&,meta/Max[Abs[meta]],{{\[Phi]0}}]]//ArrayFlatten
];
