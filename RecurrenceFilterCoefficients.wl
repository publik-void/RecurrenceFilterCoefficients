(* ::Package:: *)

BeginPackage["RecurrenceFilterCoefficients`"]

ClearAll[R5Lowpass1Coefficients, RBJLowpassCoefficients, MVHighpassCoefficients,
  MVBandshelfCoefficients];
R5Lowpass1Coefficients::usage="\!\(\*RowBox[{\"R5Lowpass1Coefficients\"}]\) \
computes coefficients for a 1st-order Reaktor 5 style lowpass filter."
RBJLowpassCoefficients::usage="\!\(\*RowBox[{\"RBJLowpassCoefficients\"}]\) \
computes coefficients for a 2nd-order Robert Bristow Johnson style lowpass \
filter."
MVHighpassCoefficients::usage"\!\(\*RowBox[{\"MVHighpassCoefficients\"}]\) \
computes coefficients for a 2nd-order Martin Vicanek style highpass filter."
MVBandshelfCoefficients::usage"\!\(\*RowBox[{\"MVBandshelfCoefficients\"}]\) \
computes coefficients for a 2nd-order Martin Vicanek style bandshelf filter."

Begin["`Private`"]

ClearAll[compilationTarget, runtimeOptions];
compilationTarget = "C";
runtimeOptions = {"Speed",
  "CatchMachineOverflow" -> False,
  "CatchMachineIntegerOverflow" -> False,
  "CompareWithTolerance" -> False,
  "EvaluateSymbolically" -> False};

With[{ro = runtimeOptions}, MVHighpassCoefficients =
  Compile[{{w, _Real}, {Q, _Real}},
   With[{q = .5/Q, phi1 = With[{x = Sin[2 Pi w*.5]}, x x]}, 
    With[{sqrta2 = Exp[w*q*(-2 Pi)], phi0 = 1 - phi1}, 
     With[{phi2 = 4*phi0*phi1, 
       a1 = sqrta2*
         If[q <= 1, -2*(Cos[2 Pi w *Sqrt[1 - q*q]]), 
          With[{e = Exp[2 Pi w *Sqrt[(q*q) - 1]]}, (-1/e) - e]], 
       a2 = sqrta2*sqrta2}, 
      With[{b0 = ((Q*.25)/
            phi1)*(Sqrt[(phi0*(With[{x = (a2 + 1 + a1)}, 
                 x x])) + (phi1*
               With[{x = (a2 + 1 - a1)}, x x]) + (a2*(-4)*
               phi2)])}, {{1, a1, a2}, {b0, -2 b0, b0}}]]]], 
  CompilationTarget -> compilationTarget,
  RuntimeOptions -> ro]];

With[{ro = runtimeOptions}, MVBandshelfCoefficients = 
  Compile[{{w, _Real}, {Q, _Real}, {A, _Real}}, 
   With[{q = .5/(Q Sqrt[A]), phi1 = With[{x = Sin[Pi w]}, x x]}, 
    With[{sqrta2 = Exp[w*q*(-2 Pi)], phi0 = 1 - phi1},
     With[{phi2 = phi0 phi1 4,
       a1 = sqrta2*
         If[q <= 1, -2*(Cos[2 Pi w *Sqrt[1 - q*q]]),
          With[{e = Exp[2 Pi w *Sqrt[(q*q) - 1]]}, (-1/e) - e]],
       a2 = sqrta2*sqrta2},
      With[{A0 = With[{x = a2 + 1 + a1}, x x],
        A1 = With[{x = a2 + 1 - a1}, x x], A2 = -4 a2},
       With[{R1 = A A (A0 phi0 + A1 phi1 + A2 phi2),
         R2 = A A (A1 - A0 - A2 4 (phi1 - phi0))},
        With[{B0 = A0, B2 = (R1 - R2 phi1 - A0)/(phi1 phi1 4)},
         With[{B1 = R2 + B0 + B2 4 (phi1 - phi0)},
          With[{sqrtB0 = Sqrt[B0], sqrtB1 = Sqrt[B1]},
           With[{W = .5 (sqrtB0 + sqrtB1)},
            With[{b0 = .5 (W + Sqrt[B2 + W W]),
              b1 = .5 (sqrtB0 - sqrtB1)},
             With[{b2 = B2/(-4 b0)},
             {{1, a1, a2}, {b0, b1, b2}}]]]]]]]]]]],
  CompilationTarget -> compilationTarget,
  RuntimeOptions -> ro]];

With[{ro = runtimeOptions}, RBJLowpassCoefficients =
  Compile[{{w, _Real}, {Q, _Real}},
   With[{cs = Cos[2 Pi w], alpha = (.5/Q) Sin[2 Pi w]},
    With[{b1 = 1 - cs},
     With[{b0 = .5 b1, a0 = alpha + 1, a1 = -2 cs,
       a2 = 1 - alpha}, {{a0, a1, a2}, {b0, b1, b0}}]]],
  CompilationTarget -> compilationTarget,
  RuntimeOptions -> ro]];

With[{ro = runtimeOptions}, R5Lowpass1Coefficients =
  Compile[{{w, _Real}},
    With[{t = Tan[Pi w]}, With[{a1 = (t-1)/(t+1), b0 = t/(t+1)}, {{1, a1},
      {b0, b0}}]],
  CompilationTarget -> compilationTarget,
  RuntimeOptions -> ro]];

End[]

EndPackage[]
