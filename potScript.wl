(* ::Package:: *)

PacletManager`PacletDirectoryAdd["/gscratch/chem/b3m2a1/Mathematica/Paclets"];
$Path=DeleteDuplicates@Append[$Path, "/gscratch/chem/b3m2a1/H5+/Common"];
Quiet@CreateDirectory["logs"];
Quiet@CreateDirectory["outs"];
$Output=OpenWrite[FileNameJoin@{"outs", "out_script.txt"}]//Prepend[$Output, #]&; 
$Messages=OpenWrite[FileNameJoin@{"logs", "log_script.txt"}]//Prepend[$Messages, #]&; ; 


<<SpectrumScript`


r1r2Grid=
  $H2DVR[
    "Grid",
    "PotentialOptimize"->False
    ]@"Points";


r1r2Pot=
  With[{H5FullVecPot=H5FullVecPot, r1r2Grid=r1r2Grid},
    ParallelTable[
      {a, s,
        Min[
          H5FullVecPot@
            Join[
              ConstantArray[a,   {Length@r1r2Grid, 1}],
              ConstantArray[s,   {Length@r1r2Grid, 1}],
              r1r2Grid,
              2
              ]
          ]
        },
      {a, -2, 2,   1},
      {s, .5, 2.5, 1}
      ]
     ];


Export["test.mx", r1r2Pot]
