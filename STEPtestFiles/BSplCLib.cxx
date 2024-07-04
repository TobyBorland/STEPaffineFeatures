gitprojects / occt.git / blob

commit
 ? search: 
 re
summary | shortlog | log | commit | commitdiff | tree
blame | history | raw | HEAD
0033261: Data Exchange, Step Import - Empty shape after reading process
[occt.git] / src / BSplCLib / BSplCLib.cxx
   1 // Created on: 1991-08-09
   2 // Created by: JCV
   3 // Copyright (c) 1991-1999 Matra Datavision
   4 // Copyright (c) 1999-2014 OPEN CASCADE SAS
   5 //
   6 // This file is part of Open CASCADE Technology software library.
   7 //
   8 // This library is free software; you can redistribute it and/or modify it under
   9 // the terms of the GNU Lesser General Public License version 2.1 as published
  10 // by the Free Software Foundation, with special exception defined in the file
  11 // OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
  12 // distribution for complete text of the license and disclaimer of any warranty.
  13 //
  14 // Alternatively, this file may be used under the terms of Open CASCADE
  15 // commercial license or contractual agreement.
  16 
  17 // Modified RLE 9 Sep 1993
  18 // pmn : modified 28-01-97  : fixed a mistake in LocateParameter (PRO6973)
  19 // pmn : modified 4-11-96   : fixed a mistake in BuildKnots (PRO6124)
  20 // pmn : modified 28-Jun-96 : fixed a mistake in AntiBoorScheme
  21 // xab : modified 15-Jun-95 : fixed a mistake in IsRational
  22 // xab : modified 15-Mar-95 : removed Epsilon comparison in IsRational
  23 //                            added RationalDerivatives.
  24 // xab : 30-Mar-95 : fixed coupling with lti in RationalDerivatives
  25 // xab : 15-Mar-96 : fixed a typo in Eval with extrapolation
  26 // jct : 15-Apr-97 : added TangExtendToConstraint
  27 // jct : 24-Apr-97 : correction on computation of Tbord and NewFlatKnots
  28 //                   in TangExtendToConstraint; Continuity can be equal to 0
  29 
  30 #include <BSplCLib.hxx>
  31 #include <ElCLib.hxx>
  32 #include <gp_Pnt.hxx>
  33 #include <math_Matrix.hxx>
  34 #include <NCollection_LocalArray.hxx>
  35 #include <PLib.hxx>
  36 #include <Precision.hxx>
  37 #include <Standard_NotImplemented.hxx>
  38 #include <math_Vector.hxx>
  39 
  40 typedef gp_Pnt Pnt;
  41 typedef gp_Vec Vec;
  42 typedef TColgp_Array1OfPnt Array1OfPnt;
  43 typedef TColStd_Array1OfReal Array1OfReal;
  44 typedef TColStd_Array1OfInteger Array1OfInteger;
  45 
  46 //=======================================================================
  47 //class : BSplCLib_LocalMatrix
  48 //purpose: Auxiliary class optimizing creation of matrix buffer for
  49 //         evaluation of bspline (using stack allocation for main matrix)
  50 //=======================================================================
  51 
  52 class BSplCLib_LocalMatrix : public math_Matrix 
  53 {
  54 public:
  55   BSplCLib_LocalMatrix (Standard_Integer DerivativeRequest, Standard_Integer Order)
  56     : math_Matrix (myBuffer, 1, DerivativeRequest + 1, 1, Order)
  57   {
  58     Standard_OutOfRange_Raise_if (DerivativeRequest > BSplCLib::MaxDegree() ||
  59         Order > BSplCLib::MaxDegree()+1 || BSplCLib::MaxDegree() > 25,
  60         "BSplCLib: bspline degree is greater than maximum supported");
  61   }
  62 
  63  private:
  64   // local buffer, to be sufficient for addressing by index [Degree+1][Degree+1]
  65   // (see math_Matrix implementation)
  66   Standard_Real myBuffer[27*27];
  67 };
  68 
  69 //=======================================================================
  70 //function : Hunt
  71 //purpose  : 
  72 //=======================================================================
  73 
  74 void BSplCLib::Hunt (const TColStd_Array1OfReal& theArray,
  75                      const Standard_Real theX,
  76                      Standard_Integer&   theXPos)
  77 {
  78   // replaced by simple dichotomy (RLE)
  79   if (theArray.First() > theX)
  80   {
  81     theXPos = theArray.Lower() - 1;
  82     return;
  83   }
  84   else if (theArray.Last() < theX)
  85   {
  86     theXPos = theArray.Upper() + 1;
  87     return;
  88   }
  89 
  90   theXPos = theArray.Lower();
  91   if (theArray.Length() <= 1)
  92   {
  93     return;
  94   }
  95 
  96   Standard_Integer aHi = theArray.Upper();
  97   while (aHi - theXPos != 1)
  98   {
  99     const Standard_Integer aMid = (aHi + theXPos) / 2;
 100     if (theArray.Value (aMid) < theX)
 101     {
 102       theXPos = aMid;
 103     }
 104     else
 105     {
 106       aHi = aMid;
 107     }
 108   }
 109 }
 110 
 111 //=======================================================================
 112 //function : FirstUKnotIndex
 113 //purpose  : 
 114 //=======================================================================
 115 
 116 Standard_Integer BSplCLib::FirstUKnotIndex (const Standard_Integer Degree,
 117                                    const TColStd_Array1OfInteger& Mults)
 118 { 
 119   Standard_Integer Index = Mults.Lower();
 120   Standard_Integer SigmaMult = Mults(Index);
 121 
 122   while (SigmaMult <= Degree) {
 123     Index++;
 124     SigmaMult += Mults (Index);
 125   }
 126   return Index;
 127 }
 128 
 129 //=======================================================================
 130 //function : LastUKnotIndex
 131 //purpose  : 
 132 //=======================================================================
 133 
 134 Standard_Integer BSplCLib::LastUKnotIndex  (const Standard_Integer Degree,
 135                                    const Array1OfInteger& Mults) 
 136 { 
 137    Standard_Integer Index = Mults.Upper();
 138    Standard_Integer SigmaMult = Mults(Index);
 139 
 140    while (SigmaMult <= Degree) {
 141      Index--;
 142      SigmaMult += Mults.Value (Index);
 143    }
 144    return Index;
 145 }
 146 
 147 //=======================================================================
 148 //function : FlatIndex
 149 //purpose  : 
 150 //=======================================================================
 151 
 152 Standard_Integer  BSplCLib::FlatIndex
 153   (const Standard_Integer Degree,
 154    const Standard_Integer Index,
 155    const TColStd_Array1OfInteger& Mults,
 156    const Standard_Boolean Periodic)
 157 {
 158   Standard_Integer i, index = Index;
 159   const Standard_Integer MLower = Mults.Lower();
 160   const Standard_Integer *pmu = &Mults(MLower);
 161   pmu -= MLower;
 162 
 163   for (i = MLower + 1; i <= Index; i++)
 164     index += pmu[i] - 1;
 165   if ( Periodic)
 166     index += Degree;
 167   else
 168     index += pmu[MLower] - 1;
 169   return index;
 170 }
 171 
 172 //=======================================================================
 173 //function : LocateParameter
 174 //purpose  : Processing of nodes with multiplicities
 175 //pmn  28-01-97 -> compute eventual of the period.
 176 //=======================================================================
 177 
 178 void BSplCLib::LocateParameter
 179 (const Standard_Integer          , //Degree,
 180  const Array1OfReal&    Knots,
 181  const Array1OfInteger& , //Mults,
 182  const Standard_Real             U,
 183  const Standard_Boolean          IsPeriodic,
 184  const Standard_Integer          FromK1,
 185  const Standard_Integer          ToK2,
 186  Standard_Integer&               KnotIndex,
 187  Standard_Real&                  NewU)
 188 {
 189   Standard_Real uf = 0, ul=1;
 190   if (IsPeriodic) {
 191     uf = Knots(Knots.Lower());
 192     ul = Knots(Knots.Upper());
 193   }
 194   BSplCLib::LocateParameter(Knots,U,IsPeriodic,FromK1,ToK2,
 195                             KnotIndex,NewU, uf, ul);
 196 }
 197 
 198 //=======================================================================
 199 //function : LocateParameter
 200 //purpose  : For plane nodes 
 201 //   pmn  28-01-97 -> There is a need of the degre to calculate
 202 //   the eventual period
 203 //=======================================================================
 204 
 205 void BSplCLib::LocateParameter
 206 (const Standard_Integer          Degree,
 207  const Array1OfReal&    Knots,
 208  const Standard_Real             U,
 209  const Standard_Boolean          IsPeriodic,
 210  const Standard_Integer          FromK1,
 211  const Standard_Integer          ToK2,
 212  Standard_Integer&               KnotIndex,
 213  Standard_Real&                  NewU)
 214 { 
 215   if (IsPeriodic)
 216     BSplCLib::LocateParameter(Knots, U, IsPeriodic, FromK1, ToK2,
 217                               KnotIndex, NewU,
 218                               Knots(Knots.Lower() + Degree),
 219                               Knots(Knots.Upper() - Degree));
 220   else 
 221     BSplCLib::LocateParameter(Knots, U, IsPeriodic, FromK1, ToK2,
 222                               KnotIndex, NewU,
 223                               0.,
 224                               1.);
 225 }
 226 
 227 //=======================================================================
 228 //function : LocateParameter
 229 //purpose  : Effective computation
 230 // pmn 28-01-97 : Add limits of the period as input argument,  
 231 //                as it is impossible to produce them at this level.
 232 //=======================================================================
 233 
 234 void BSplCLib::LocateParameter 
 235 (const TColStd_Array1OfReal& Knots,
 236  const Standard_Real         U,
 237  const Standard_Boolean      IsPeriodic,
 238  const Standard_Integer      FromK1,
 239  const Standard_Integer      ToK2,
 240  Standard_Integer&           KnotIndex,
 241  Standard_Real&              NewU,
 242  const Standard_Real         UFirst,
 243  const Standard_Real         ULast)
 244 {
 245   /*
 246   Let Knots are distributed as follows (the array is sorted in ascending order):
 247     
 248       K1, K1,..., K1, K1, K2, K2,..., K2, K2,..., Kn, Kn,..., Kn
 249            M1 times             M2 times             Mn times
 250 
 251   NbKnots = sum(M1+M2+...+Mn)
 252   If U <= K1 then KnotIndex should be equal to M1.
 253   If U >= Kn then KnotIndex should be equal to NbKnots-Mn-1.
 254   If Ki <= U < K(i+1) then KnotIndex should be equal to sum (M1+M2+...+Mi).
 255   */
 256 
 257   Standard_Integer First,Last;
 258   if (FromK1 < ToK2) {
 259     First = FromK1;
 260     Last  = ToK2;
 261   }
 262   else {
 263     First = ToK2;
 264     Last  = FromK1;
 265   }
 266   Standard_Integer Last1 = Last - 1;
 267   NewU = U;
 268   if (IsPeriodic && (NewU < UFirst || NewU > ULast))
 269     NewU = ElCLib::InPeriod(NewU, UFirst, ULast);
 270   
 271   BSplCLib::Hunt (Knots, NewU, KnotIndex);
 272   
 273   Standard_Real val;
 274   const Standard_Integer  KLower = Knots.Lower(),
 275                           KUpper = Knots.Upper();
 276 
 277   const Standard_Real Eps = Epsilon(Min(Abs(Knots(KUpper)), Abs(U)));
 278 
 279   const Standard_Real *knots = &Knots(KLower);
 280   knots -= KLower;
 281   if ( KnotIndex < Knots.Upper()) {
 282     val = NewU - knots[KnotIndex + 1];
 283     if (val < 0) val = - val;
 284     // <= to be coherent with Segment where Eps corresponds to a bit of error.
 285     if (val <= Eps) KnotIndex++; 
 286   }
 287   if (KnotIndex < First) KnotIndex = First;
 288   if (KnotIndex > Last1) KnotIndex = Last1;
 289   
 290   if (KnotIndex != Last1) {
 291     Standard_Real K1 = knots[KnotIndex];
 292     Standard_Real K2 = knots[KnotIndex + 1];
 293     val = K2 - K1;
 294     if (val < 0) val = - val;
 295 
 296     while (val <= Eps) {
 297       KnotIndex++;
 298 
 299       if(KnotIndex >= Knots.Upper())
 300         break;
 301 
 302       K1 = K2;
 303       K2 = knots[KnotIndex + 1];
 304       val = K2 - K1;
 305       if (val < 0) val = - val;
 306     }
 307   }
 308 }
 309 
 310 //=======================================================================
 311 //function : LocateParameter
 312 //purpose  : the index is recomputed only if out of range
 313 //pmn  28-01-97 -> eventual computation of the period.
 314 //=======================================================================
 315 
 316 void BSplCLib::LocateParameter 
 317 (const Standard_Integer         Degree,
 318  const TColStd_Array1OfReal&    Knots,
 319  const TColStd_Array1OfInteger* Mults,
 320  const Standard_Real            U,
 321  const Standard_Boolean         Periodic,
 322  Standard_Integer&              KnotIndex,
 323  Standard_Real&                 NewU) 
 324 {
 325   Standard_Integer first,last;
 326   if (Mults) {
 327     if (Periodic) {
 328       first = Knots.Lower();
 329       last  = Knots.Upper();
 330     }
 331     else {
 332       first = FirstUKnotIndex(Degree,*Mults);
 333       last  = LastUKnotIndex (Degree,*Mults);
 334     }
 335   }
 336   else {
 337     first = Knots.Lower() + Degree;
 338     last  = Knots.Upper() - Degree;
 339   }
 340   if ( KnotIndex < first || KnotIndex > last)
 341     BSplCLib::LocateParameter(Knots, U, Periodic, first, last,
 342                               KnotIndex, NewU, Knots(first), Knots(last));
 343   else
 344     NewU = U;
 345 }
 346 
 347 //=======================================================================
 348 //function : MaxKnotMult
 349 //purpose  : 
 350 //=======================================================================
 351 
 352 Standard_Integer BSplCLib::MaxKnotMult
 353 (const Array1OfInteger& Mults,
 354  const Standard_Integer          FromK1,
 355  const Standard_Integer          ToK2)
 356 {
 357   Standard_Integer MLower = Mults.Lower();
 358   const Standard_Integer *pmu = &Mults(MLower);
 359   pmu -= MLower;
 360   Standard_Integer MaxMult = pmu[FromK1];
 361 
 362   for (Standard_Integer i = FromK1; i <= ToK2; i++) {
 363     if (MaxMult < pmu[i]) MaxMult = pmu[i];
 364   }
 365   return MaxMult;
 366 }
 367 
 368 //=======================================================================
 369 //function : MinKnotMult
 370 //purpose  : 
 371 //=======================================================================
 372 
 373 Standard_Integer BSplCLib::MinKnotMult
 374 (const Array1OfInteger& Mults,
 375  const Standard_Integer          FromK1,
 376  const Standard_Integer          ToK2)
 377 {
 378   Standard_Integer MLower = Mults.Lower();
 379   const Standard_Integer *pmu = &Mults(MLower);
 380   pmu -= MLower;
 381   Standard_Integer MinMult = pmu[FromK1];
 382 
 383   for (Standard_Integer i = FromK1; i <= ToK2; i++) {
 384     if (MinMult > pmu[i]) MinMult = pmu[i];
 385   }
 386   return MinMult;
 387 }
 388 
 389 //=======================================================================
 390 //function : NbPoles
 391 //purpose  : 
 392 //=======================================================================
 393 
 394 Standard_Integer BSplCLib::NbPoles(const Standard_Integer Degree,
 395                                    const Standard_Boolean Periodic,
 396                                    const TColStd_Array1OfInteger& Mults)
 397 {
 398   Standard_Integer i,sigma = 0;
 399   Standard_Integer f = Mults.Lower();
 400   Standard_Integer l = Mults.Upper();
 401   const Standard_Integer * pmu = &Mults(f);
 402   pmu -= f;
 403   Standard_Integer Mf = pmu[f];
 404   Standard_Integer Ml = pmu[l];
 405   if (Mf <= 0) return 0;
 406   if (Ml <= 0) return 0;
 407   if (Periodic) {
 408     if (Mf > Degree) return 0;
 409     if (Ml > Degree) return 0;
 410     if (Mf != Ml   ) return 0;
 411     sigma = Mf;
 412   }
 413   else {
 414     Standard_Integer Deg1 = Degree + 1;
 415     if (Mf > Deg1) return 0;
 416     if (Ml > Deg1) return 0;
 417     sigma = Mf + Ml - Deg1;
 418   }
 419     
 420   for (i = f + 1; i < l; i++) {
 421     if (pmu[i] <= 0    ) return 0;
 422     if (pmu[i] > Degree) return 0;
 423     sigma += pmu[i];
 424   }
 425   return sigma;
 426 }
 427 
 428 //=======================================================================
 429 //function : KnotSequenceLength
 430 //purpose  : 
 431 //=======================================================================
 432 
 433 Standard_Integer BSplCLib::KnotSequenceLength
 434 (const TColStd_Array1OfInteger& Mults,
 435  const Standard_Integer         Degree,
 436  const Standard_Boolean         Periodic)
 437 {
 438   Standard_Integer i,l = 0;
 439   Standard_Integer MLower = Mults.Lower();
 440   Standard_Integer MUpper = Mults.Upper();
 441   const Standard_Integer * pmu = &Mults(MLower);
 442   pmu -= MLower;
 443 
 444   for (i = MLower; i <= MUpper; i++)
 445     l += pmu[i];
 446   if (Periodic) l += 2 * (Degree + 1 - pmu[MLower]);
 447   return l;
 448 }
 449 
 450 //=======================================================================
 451 //function : KnotSequence
 452 //purpose  : 
 453 //=======================================================================
 454 
 455 void BSplCLib::KnotSequence 
 456 (const TColStd_Array1OfReal&    Knots,
 457  const TColStd_Array1OfInteger& Mults,
 458  TColStd_Array1OfReal&          KnotSeq,
 459  const Standard_Boolean         Periodic)
 460 {
 461   BSplCLib::KnotSequence(Knots,Mults,0,Periodic,KnotSeq);
 462 }
 463 
 464 //=======================================================================
 465 //function : KnotSequence
 466 //purpose  : 
 467 //=======================================================================
 468 
 469 void BSplCLib::KnotSequence 
 470 (const TColStd_Array1OfReal&    Knots,
 471  const TColStd_Array1OfInteger& Mults,
 472  const Standard_Integer         Degree,
 473  const Standard_Boolean         Periodic,
 474  TColStd_Array1OfReal&          KnotSeq)
 475 {
 476   Standard_Real K;
 477   Standard_Integer Mult;
 478   Standard_Integer MLower = Mults.Lower();
 479   const Standard_Integer * pmu = &Mults(MLower);
 480   pmu -= MLower;
 481   Standard_Integer KLower = Knots.Lower();
 482   Standard_Integer KUpper = Knots.Upper();
 483   const Standard_Real * pkn = &Knots(KLower);
 484   pkn -= KLower;
 485   Standard_Integer M1 = Degree + 1 - pmu[MLower];  // for periodic
 486   Standard_Integer i,j,index = Periodic ? M1 + 1 : 1;
 487 
 488   for (i = KLower; i <= KUpper; i++) {
 489     Mult = pmu[i];
 490     K    = pkn[i];
 491 
 492     for (j = 1; j <= Mult; j++) { 
 493       KnotSeq (index) = K;   
 494       index++;
 495     }
 496   }
 497   if (Periodic) {
 498     Standard_Real period = pkn[KUpper] - pkn[KLower];
 499     Standard_Integer m;
 500     m = 1;
 501     j = KUpper - 1;
 502 
 503     for (i = M1; i >= 1; i--) {
 504       KnotSeq(i) = pkn[j] - period;
 505       m++;
 506       if (m > pmu[j]) {
 507         j--;
 508         m = 1;
 509       }
 510     }
 511     m = 1;
 512     j = KLower + 1;
 513 
 514     for (i = index; i <= KnotSeq.Upper(); i++) {
 515       KnotSeq(i) = pkn[j] + period;
 516       m++;
 517       if (m > pmu[j]) {
 518         j++;
 519         m = 1;
 520       }
 521     }
 522   }
 523 }
 524 
 525 //=======================================================================
 526 //function : KnotsLength
 527 //purpose  : 
 528 //=======================================================================
 529  Standard_Integer BSplCLib::KnotsLength(const TColStd_Array1OfReal& SeqKnots,
 530 //                                      const Standard_Boolean Periodic)
 531                                         const Standard_Boolean )
 532 {
 533   Standard_Integer sizeMult = 1; 
 534   Standard_Real val = SeqKnots(1);
 535   for (Standard_Integer jj=2;
 536        jj<=SeqKnots.Length();jj++)
 537     {
 538       // test on strict equality on nodes
 539       if (SeqKnots(jj)!=val)
 540         {
 541           val = SeqKnots(jj);
 542           sizeMult++;
 543         }
 544     }
 545   return sizeMult;
 546 }
 547 
 548 //=======================================================================
 549 //function : Knots
 550 //purpose  : 
 551 //=======================================================================
 552 void BSplCLib::Knots(const TColStd_Array1OfReal& SeqKnots, 
 553                      TColStd_Array1OfReal &knots,
 554                      TColStd_Array1OfInteger &mult,
 555 //                   const Standard_Boolean Periodic)
 556                      const Standard_Boolean )
 557 {
 558   Standard_Real val = SeqKnots(1);
 559   Standard_Integer kk=1;
 560   knots(kk) = val;
 561   mult(kk)  = 1;
 562 
 563   for (Standard_Integer jj=2;jj<=SeqKnots.Length();jj++)
 564     {
 565       // test on strict equality on nodes
 566       if (SeqKnots(jj)!=val)
 567         {
 568           val = SeqKnots(jj);
 569           kk++;
 570           knots(kk) = val;
 571           mult(kk)  = 1;
 572         }
 573       else
 574         {
 575           mult(kk)++;
 576         }
 577     }
 578 }
 579 
 580 //=======================================================================
 581 //function : KnotForm
 582 //purpose  : 
 583 //=======================================================================
 584 
 585 BSplCLib_KnotDistribution BSplCLib::KnotForm
 586 (const Array1OfReal& Knots,
 587  const Standard_Integer       FromK1,
 588  const Standard_Integer       ToK2)
 589 {
 590    Standard_Real DU0,DU1,Ui,Uj,Eps0,val;
 591    BSplCLib_KnotDistribution  KForm = BSplCLib_Uniform;
 592 
 593    if (FromK1 + 1 > Knots.Upper())
 594    {
 595      return BSplCLib_Uniform;
 596    }
 597 
 598    Ui  = Knots(FromK1);
 599    if (Ui < 0) Ui = - Ui;
 600    Uj  = Knots(FromK1 + 1);
 601    if (Uj < 0) Uj = - Uj;
 602    DU0 = Uj - Ui;
 603    if (DU0 < 0) DU0 = - DU0;
 604    Eps0 = Epsilon (Ui) + Epsilon (Uj) + Epsilon (DU0);
 605    Standard_Integer i = FromK1 + 1;
 606 
 607    while (KForm != BSplCLib_NonUniform && i < ToK2) {
 608      Ui = Knots(i);
 609      if (Ui < 0) Ui = - Ui;
 610      i++;
 611      Uj = Knots(i);
 612      if (Uj < 0) Uj = - Uj;
 613      DU1 = Uj - Ui;
 614      if (DU1 < 0) DU1 = - DU1;
 615      val = DU1 - DU0;
 616      if (val < 0) val = -val;
 617      if (val > Eps0) KForm = BSplCLib_NonUniform;
 618      DU0 = DU1;
 619      Eps0 = Epsilon (Ui) + Epsilon (Uj) + Epsilon (DU0);
 620    }
 621    return KForm;
 622 }
 623 
 624 //=======================================================================
 625 //function : MultForm
 626 //purpose  : 
 627 //=======================================================================
 628 
 629 BSplCLib_MultDistribution BSplCLib::MultForm
 630 (const Array1OfInteger& Mults,
 631  const Standard_Integer          FromK1,
 632  const Standard_Integer          ToK2)
 633 {
 634   Standard_Integer First,Last;
 635   if (FromK1 < ToK2) {
 636     First = FromK1;
 637     Last  = ToK2;
 638   }
 639   else {
 640     First = ToK2;
 641     Last  = FromK1;
 642   }
 643   if (First + 1 > Mults.Upper())
 644   {
 645     return BSplCLib_Constant;
 646   }
 647 
 648   Standard_Integer FirstMult = Mults(First);
 649   BSplCLib_MultDistribution MForm = BSplCLib_Constant;
 650   Standard_Integer i    = First + 1;
 651   Standard_Integer Mult = Mults(i);
 652   
 653 //  while (MForm != BSplCLib_NonUniform && i <= Last) { ???????????JR????????
 654   while (MForm != BSplCLib_NonConstant && i <= Last) {
 655     if (i == First + 1) {  
 656       if (Mult != FirstMult)      MForm = BSplCLib_QuasiConstant;
 657     }
 658     else if (i == Last)  {
 659       if (MForm == BSplCLib_QuasiConstant) {
 660         if (FirstMult != Mults(i))  MForm = BSplCLib_NonConstant;
 661       }
 662       else {
 663         if (Mult != Mults(i))       MForm = BSplCLib_NonConstant;
 664       }
 665     }
 666     else {
 667       if (Mult != Mults(i))         MForm = BSplCLib_NonConstant;
 668       Mult = Mults(i);
 669     }
 670     i++;
 671   }
 672   return MForm;
 673 }
 674 
 675 //=======================================================================
 676 //function : KnotAnalysis
 677 //purpose  : 
 678 //=======================================================================
 679 
 680 void BSplCLib::KnotAnalysis (const Standard_Integer         Degree,
 681                              const Standard_Boolean         Periodic,
 682                              const TColStd_Array1OfReal&    CKnots,
 683                              const TColStd_Array1OfInteger& CMults,
 684                              GeomAbs_BSplKnotDistribution&  KnotForm,
 685                              Standard_Integer&              MaxKnotMult)
 686 {
 687   KnotForm = GeomAbs_NonUniform;
 688 
 689   BSplCLib_KnotDistribution KSet = 
 690     BSplCLib::KnotForm (CKnots, 1, CKnots.Length());
 691   
 692 
 693   if (KSet == BSplCLib_Uniform) {
 694     BSplCLib_MultDistribution MSet =
 695       BSplCLib::MultForm (CMults, 1, CMults.Length());
 696     switch (MSet) {
 697     case BSplCLib_NonConstant   :       
 698       break;
 699     case BSplCLib_Constant      : 
 700       if (CKnots.Length() == 2) {
 701         KnotForm = GeomAbs_PiecewiseBezier;
 702       }
 703       else {
 704         if (CMults (1) == 1)  KnotForm = GeomAbs_Uniform;   
 705       }
 706       break;
 707     case BSplCLib_QuasiConstant :   
 708       if (CMults (1) == Degree + 1) {
 709         Standard_Real M = CMults (2);
 710         if (M == Degree )   KnotForm = GeomAbs_PiecewiseBezier;
 711         else if  (M == 1)   KnotForm = GeomAbs_QuasiUniform;
 712       }
 713       break;
 714     }
 715   }
 716 
 717   Standard_Integer FirstKM = 
 718     Periodic ? CKnots.Lower() :  BSplCLib::FirstUKnotIndex (Degree,CMults);
 719   Standard_Integer LastKM =
 720     Periodic ? CKnots.Upper() :  BSplCLib::LastUKnotIndex (Degree,CMults);
 721   MaxKnotMult = 0;
 722   if (LastKM - FirstKM != 1) {
 723     Standard_Integer Multi;
 724     for (Standard_Integer i = FirstKM + 1; i < LastKM; i++) {
 725       Multi = CMults (i);
 726       MaxKnotMult = Max (MaxKnotMult, Multi);
 727     }
 728   }
 729 }
 730 
 731 //=======================================================================
 732 //function : Reparametrize
 733 //purpose  : 
 734 //=======================================================================
 735 
 736 void BSplCLib::Reparametrize
 737 (const Standard_Real      U1,
 738  const Standard_Real      U2,
 739  Array1OfReal&   Knots)
 740 {
 741   Standard_Integer Lower  = Knots.Lower();
 742   Standard_Integer Upper  = Knots.Upper();
 743   Standard_Real UFirst    = Min (U1, U2);
 744   Standard_Real ULast     = Max (U1, U2);
 745   Standard_Real NewLength = ULast - UFirst;
 746   BSplCLib_KnotDistribution KSet = BSplCLib::KnotForm (Knots, Lower, Upper);
 747   if (KSet == BSplCLib_Uniform) {
 748     Standard_Real DU = NewLength / (Upper - Lower);
 749     Knots (Lower) = UFirst;
 750 
 751     for (Standard_Integer i = Lower + 1; i <= Upper; i++) {
 752       Knots (i) = Knots (i-1) + DU;
 753     }
 754   }
 755   else {
 756     Standard_Real K2;
 757     Standard_Real Ratio;
 758     Standard_Real K1 = Knots (Lower);
 759     Standard_Real Length = Knots (Upper) - Knots (Lower);
 760     Knots (Lower) = UFirst;
 761 
 762     for (Standard_Integer i = Lower + 1; i <= Upper; i++) {
 763       K2 = Knots (i);
 764       Ratio = (K2 - K1) / Length;
 765       Knots (i) = Knots (i-1) + (NewLength * Ratio);
 766 
 767       //for CheckCurveData
 768       Standard_Real Eps = Epsilon( Abs(Knots(i-1)) );
 769       if (Knots(i) - Knots(i-1) <= Eps)
 770         Knots(i) = NextAfter (Knots(i-1) + Eps, RealLast());
 771 
 772       K1 = K2;
 773     }
 774   }
 775 }
 776 
 777 //=======================================================================
 778 //function : Reverse
 779 //purpose  : 
 780 //=======================================================================
 781 
 782 void  BSplCLib::Reverse(TColStd_Array1OfReal& Knots)
 783 {
 784   Standard_Integer first = Knots.Lower();
 785   Standard_Integer last  = Knots.Upper();
 786   Standard_Real kfirst = Knots(first);
 787   Standard_Real klast = Knots(last);
 788   Standard_Real tfirst = kfirst;
 789   Standard_Real tlast  = klast;
 790   first++;
 791   last--;
 792 
 793   while (first <= last) {
 794     tfirst += klast - Knots(last);
 795     tlast  -= Knots(first) - kfirst;
 796     kfirst = Knots(first);
 797     klast  = Knots(last);
 798     Knots(first) = tfirst;
 799     Knots(last)  = tlast;
 800     first++;
 801     last--;
 802   }
 803 }
 804 
 805 //=======================================================================
 806 //function : Reverse
 807 //purpose  : 
 808 //=======================================================================
 809 
 810 void  BSplCLib::Reverse(TColStd_Array1OfInteger& Mults)
 811 {
 812   Standard_Integer first = Mults.Lower();
 813   Standard_Integer last  = Mults.Upper();
 814   Standard_Integer temp;
 815 
 816   while (first < last) {
 817     temp = Mults(first);
 818     Mults(first) = Mults(last);
 819     Mults(last) = temp;
 820     first++;
 821     last--;
 822   }
 823 }
 824 
 825 //=======================================================================
 826 //function : Reverse
 827 //purpose  : 
 828 //=======================================================================
 829 
 830 void  BSplCLib::Reverse(TColStd_Array1OfReal& Weights,
 831                         const Standard_Integer L)
 832 {
 833   Standard_Integer i, l = L;
 834   l = Weights.Lower()+(l-Weights.Lower())%(Weights.Upper()-Weights.Lower()+1);
 835 
 836   TColStd_Array1OfReal temp(0,Weights.Length()-1);
 837 
 838   for (i = Weights.Lower(); i <= l; i++)
 839     temp(l-i) = Weights(i);
 840 
 841   for (i = l+1; i <= Weights.Upper(); i++)
 842     temp(l-Weights.Lower()+Weights.Upper()-i+1) = Weights(i);
 843 
 844   for (i = Weights.Lower(); i <= Weights.Upper(); i++)
 845     Weights(i) = temp(i-Weights.Lower());
 846 }
 847 
 848 //=======================================================================
 849 //function : IsRational
 850 //purpose  : 
 851 //=======================================================================
 852 
 853 Standard_Boolean  BSplCLib::IsRational(const TColStd_Array1OfReal& Weights,
 854                                        const Standard_Integer I1,
 855                                        const Standard_Integer I2,
 856 //                                     const Standard_Real Epsi)
 857                                        const Standard_Real )
 858 {
 859   Standard_Integer i, f = Weights.Lower(), l = Weights.Length();
 860   Standard_Integer I3 = I2 - f;
 861   const Standard_Real * WG = &Weights(f);
 862   WG -= f;
 863 
 864   for (i = I1 - f; i < I3; i++) {
 865     if (WG[f + (i % l)] != WG[f + ((i + 1) % l)]) return Standard_True;
 866   }
 867   return Standard_False ;
 868 }
 869 
 870 //=======================================================================
 871 //function : Eval
 872 //purpose  : evaluate point and derivatives
 873 //=======================================================================
 874 
 875 void  BSplCLib::Eval(const Standard_Real U,
 876                      const Standard_Integer Degree,
 877                      Standard_Real& Knots, 
 878                      const Standard_Integer Dimension, 
 879                      Standard_Real& Poles)
 880 {
 881   Standard_Integer step,i,Dms,Dm1,Dpi,Sti;
 882   Standard_Real X, Y, *poles, *knots = &Knots;
 883   Dm1 = Dms = Degree;
 884   Dm1--;
 885   Dms++;
 886   switch (Dimension) { 
 887 
 888   case 1 : {
 889     
 890     for (step = - 1; step < Dm1; step++) {
 891       Dms--;
 892       poles = &Poles;
 893       Dpi   = Dm1;
 894       Sti   = step;
 895       
 896       for (i = 0; i < Dms; i++) {
 897         Dpi++;
 898         Sti++;
 899         X = (knots[Dpi] - U) / (knots[Dpi] - knots[Sti]);
 900         Y = 1 - X;
 901         poles[0] *= X; poles[0] += Y * poles[1];
 902         poles += 1;
 903       }
 904     }
 905     break;
 906   }
 907   case 2 : {
 908     
 909     for (step = - 1; step < Dm1; step++) {
 910       Dms--;
 911       poles = &Poles;
 912       Dpi   = Dm1;
 913       Sti   = step;
 914       
 915       for (i = 0; i < Dms; i++) {
 916         Dpi++;
 917         Sti++;
 918         X = (knots[Dpi] - U) / (knots[Dpi] - knots[Sti]);
 919         Y = 1 - X;
 920         poles[0] *= X; poles[0] += Y * poles[2];
 921         poles[1] *= X; poles[1] += Y * poles[3];
 922         poles += 2;
 923       }
 924     }
 925     break;
 926   }
 927   case 3 : {
 928     
 929     for (step = - 1; step < Dm1; step++) {
 930       Dms--;
 931       poles = &Poles;
 932       Dpi   = Dm1;
 933       Sti   = step;
 934       
 935       for (i = 0; i < Dms; i++) {
 936         Dpi++;
 937         Sti++;
 938         X = (knots[Dpi] - U) / (knots[Dpi] - knots[Sti]);
 939         Y = 1 - X;
 940         poles[0] *= X; poles[0] += Y * poles[3];
 941         poles[1] *= X; poles[1] += Y * poles[4];
 942         poles[2] *= X; poles[2] += Y * poles[5];
 943         poles += 3;
 944       }
 945     }
 946     break;
 947   }
 948   case 4 : {
 949     
 950     for (step = - 1; step < Dm1; step++) {
 951       Dms--;
 952       poles = &Poles;
 953       Dpi   = Dm1;
 954       Sti   = step;
 955       
 956       for (i = 0; i < Dms; i++) {
 957         Dpi++;
 958         Sti++;
 959         X = (knots[Dpi] - U) / (knots[Dpi] - knots[Sti]);
 960         Y = 1 - X;
 961         poles[0] *= X; poles[0] += Y * poles[4];
 962         poles[1] *= X; poles[1] += Y * poles[5];
 963         poles[2] *= X; poles[2] += Y * poles[6];
 964         poles[3] *= X; poles[3] += Y * poles[7];
 965         poles += 4;
 966       }
 967     }
 968     break;
 969   }
 970     default : {
 971       Standard_Integer k;
 972       
 973       for (step = - 1; step < Dm1; step++) {
 974         Dms--;
 975         poles = &Poles;
 976         Dpi   = Dm1;
 977         Sti   = step;
 978         
 979         for (i = 0; i < Dms; i++) {
 980           Dpi++;
 981           Sti++;
 982           X = (knots[Dpi] - U) / (knots[Dpi] - knots[Sti]);
 983           Y = 1 - X;
 984           
 985           for (k = 0; k < Dimension; k++) {
 986             poles[k] *= X;
 987             poles[k] += Y * poles[k + Dimension];
 988           }
 989           poles += Dimension;
 990         }
 991       }
 992     }
 993   }
 994 }
 995 
 996 //=======================================================================
 997 //function : BoorScheme
 998 //purpose  : 
 999 //=======================================================================
1000 
1001 void  BSplCLib::BoorScheme(const Standard_Real U,
1002                            const Standard_Integer Degree,
1003                            Standard_Real& Knots, 
1004                            const Standard_Integer Dimension, 
1005                            Standard_Real& Poles, 
1006                            const Standard_Integer Depth, 
1007                            const Standard_Integer Length)
1008 {
1009   //
1010   // Compute the values
1011   //
1012   //  P(i,j) (U).
1013   //
1014   // for i = 0 to Depth, 
1015   // j = 0 to Length - i
1016   //
1017   //  The Boor scheme is :
1018   //
1019   //  P(0,j) = Pole(j)
1020   //  P(i,j) = x * P(i-1,j) + (1-x) * P(i-1,j+1)
1021   //
1022   //    where x = (knot(i+j+Degree) - U) / (knot(i+j+Degree) - knot(i+j))
1023   //
1024   //
1025   //  The values are stored in the array Poles
1026   //  They are alternatively written if the odd and even positions.
1027   //
1028   //  The successives contents of the array are
1029   //   ***** means unitialised, l = Degree + Length
1030   //
1031   //  P(0,0) ****** P(0,1) ...... P(0,l-1) ******** P(0,l)
1032   //  P(0,0) P(1,0) P(0,1) ...... P(0,l-1) P(1,l-1) P(0,l)
1033   //  P(0,0) P(1,0) P(2,0) ...... P(2,l-1) P(1,l-1) P(0,l)
1034   //
1035 
1036   Standard_Integer i,k,step;
1037   Standard_Real *knots = &Knots;
1038   Standard_Real *pole, *firstpole = &Poles - 2 * Dimension;
1039   // the steps of recursion
1040 
1041   for (step = 0; step < Depth; step++) {
1042     firstpole += Dimension;
1043     pole = firstpole;
1044     // compute the new row of poles
1045 
1046     for (i = step; i < Length; i++) {
1047       pole += 2 * Dimension;
1048       // coefficient
1049       Standard_Real X = (knots[i+Degree-step] - U) 
1050         / (knots[i+Degree-step] - knots[i]);
1051       Standard_Real Y = 1. - X;
1052       // Value
1053       // P(i,j) = X * P(i-1,j) + (1-X) * P(i-1,j+1)
1054 
1055       for (k = 0; k < Dimension; k++)
1056         pole[k] = X * pole[k - Dimension] + Y * pole[k + Dimension];
1057     }
1058   }
1059 }
1060 
1061 //=======================================================================
1062 //function : AntiBoorScheme
1063 //purpose  : 
1064 //=======================================================================
1065 
1066 Standard_Boolean  BSplCLib::AntiBoorScheme(const Standard_Real    U,
1067                                            const Standard_Integer Degree,
1068                                            Standard_Real&         Knots, 
1069                                            const Standard_Integer Dimension, 
1070                                            Standard_Real&         Poles, 
1071                                            const Standard_Integer Depth, 
1072                                            const Standard_Integer Length,
1073                                            const Standard_Real    Tolerance)
1074 {
1075   // do the Boor scheme reverted.
1076 
1077   Standard_Integer i,k,step, half_length;
1078   Standard_Real *knots = &Knots;
1079   Standard_Real z,X,Y,*pole, *firstpole = &Poles + (Depth-1) * Dimension;
1080 
1081   // Test the special case length = 1 
1082   // only verification of the central point
1083 
1084   if (Length == 1) {
1085     X = (knots[Degree] - U) / (knots[Degree] - knots[0]);
1086     Y = 1. - X;
1087 
1088     for (k = 0; k < Dimension; k++) {
1089       z = X * firstpole[k] + Y * firstpole[k+2*Dimension];
1090       if (Abs(z - firstpole[k+Dimension]) > Tolerance) 
1091         return Standard_False;
1092     }
1093     return Standard_True;
1094   }
1095 
1096   // General case
1097   // the steps of recursion
1098 
1099   for (step = Depth-1; step >= 0; step--) {
1100     firstpole -= Dimension;
1101     pole = firstpole;
1102 
1103     // first step from left to right
1104 
1105     for (i = step; i < Length-1; i++) {
1106       pole += 2 * Dimension;
1107 
1108       X = (knots[i+Degree-step] - U) / (knots[i+Degree-step] - knots[i]);
1109       Y = 1. - X;
1110 
1111       for (k = 0; k < Dimension; k++)
1112         pole[k+Dimension] = (pole[k] - X*pole[k-Dimension]) / Y;
1113 
1114     }
1115 
1116     // second step from right to left
1117     pole += 4* Dimension;
1118     half_length = (Length - 1 + step) / 2  ;
1119     //
1120     // only do half of the way from right to left 
1121     // otherwise it start degenerating because of 
1122     // overflows
1123     // 
1124 
1125     for (i = Length-1; i > half_length ; i--) {
1126       pole -= 2 * Dimension;
1127 
1128       // coefficient
1129       X = (knots[i+Degree-step] - U) / (knots[i+Degree-step] - knots[i]);
1130       Y = 1. - X;
1131 
1132       for (k = 0; k < Dimension; k++) {
1133         z = (pole[k] - Y * pole[k+Dimension]) / X;
1134         if (Abs(z-pole[k-Dimension]) > Tolerance) 
1135           return Standard_False;
1136         pole[k-Dimension] += z;
1137         pole[k-Dimension] /= 2.;
1138       }
1139     }
1140   }
1141   return Standard_True;
1142 }
1143 
1144 //=======================================================================
1145 //function : Derivative
1146 //purpose  : 
1147 //=======================================================================
1148 
1149 void  BSplCLib::Derivative(const Standard_Integer Degree, 
1150                            Standard_Real& Knots, 
1151                            const Standard_Integer Dimension, 
1152                            const Standard_Integer Length, 
1153                            const Standard_Integer Order, 
1154                            Standard_Real& Poles)
1155 {
1156   Standard_Integer i,k,step,span = Degree;
1157   Standard_Real *knot = &Knots;
1158 
1159   for (step = 1; step <= Order; step++) {
1160     Standard_Real* pole = &Poles;
1161 
1162     for (i = step; i < Length; i++) {
1163       Standard_Real coef = - span / (knot[i+span] - knot[i]);
1164 
1165       for (k = 0; k < Dimension; k++) {
1166         pole[k] -= pole[k+Dimension];
1167         pole[k] *= coef;
1168       }
1169       pole += Dimension;
1170     }
1171     span--;
1172   }
1173 }
1174 
1175 //=======================================================================
1176 //function : Bohm
1177 //purpose  : 
1178 //=======================================================================
1179 
1180 void  BSplCLib::Bohm(const Standard_Real U,
1181                      const Standard_Integer Degree,
1182                      const Standard_Integer N,
1183                      Standard_Real& Knots,
1184                      const Standard_Integer Dimension,
1185                      Standard_Real& Poles)
1186 {
1187   // First phase independent of U, compute the poles of the derivatives
1188   Standard_Integer i,j,iDim,min,Dmi,DDmi,jDmi,Degm1;
1189   Standard_Real *knot = &Knots, *pole, coef, *tbis, *psav, *psDD, *psDDmDim;
1190   psav     = &Poles;
1191   if (N < Degree) min = N;
1192   else            min = Degree;
1193   Degm1 = Degree - 1;
1194   DDmi = (Degree << 1) + 1;
1195   switch (Dimension) { 
1196   case 1 : {
1197     psDD     = psav + Degree;
1198     psDDmDim = psDD - 1;
1199     
1200     for (i = 0; i < Degree; i++) {
1201       DDmi--;
1202       pole = psDD;
1203       tbis = psDDmDim;
1204       jDmi = DDmi;
1205       
1206       for (j = Degm1; j >= i; j--) {
1207         jDmi--;
1208         *pole -= *tbis;
1209   *pole = (knot[jDmi] == knot[j]) ? 0.0 :  *pole / (knot[jDmi] - knot[j]);
1210         pole--;
1211         tbis--;
1212       }
1213     }
1214     // Second phase, dependant of U
1215     iDim = - 1;
1216     
1217     for (i = 0; i < Degree; i++) {
1218       iDim += 1;
1219       pole  = psav + iDim;
1220       tbis  = pole + 1;
1221       coef  = U - knot[i];
1222       
1223       for (j = i; j >= 0; j--) {
1224         *pole += coef * (*tbis);
1225         pole--;
1226         tbis--;
1227       }
1228     }
1229     // multiply by the degrees
1230     coef = Degree;
1231     Dmi  = Degree;
1232     pole = psav + 1;
1233     
1234     for (i = 1; i <= min; i++) {
1235       *pole *= coef; pole++;
1236       Dmi--;
1237       coef  *= Dmi;
1238     }
1239     break;
1240   }
1241   case 2 : {
1242     psDD     = psav + (Degree << 1);
1243     psDDmDim = psDD - 2;
1244     
1245     for (i = 0; i < Degree; i++) {
1246       DDmi--;
1247       pole = psDD;
1248       tbis = psDDmDim;
1249       jDmi = DDmi;
1250       
1251       for (j = Degm1; j >= i; j--) {
1252         jDmi--;
1253         coef   = (knot[jDmi] == knot[j]) ? 0.0 : 1. / (knot[jDmi] - knot[j]);
1254         *pole -= *tbis; *pole *= coef; pole++; tbis++;
1255         *pole -= *tbis; *pole *= coef;
1256         pole  -= 3;
1257         tbis  -= 3;
1258       }
1259     }
1260     // Second phase, dependant of U
1261     iDim = - 2;
1262     
1263     for (i = 0; i < Degree; i++) {
1264       iDim += 2;
1265       pole  = psav + iDim;
1266       tbis  = pole + 2;
1267       coef  = U - knot[i];
1268       
1269       for (j = i; j >= 0; j--) {
1270         *pole += coef * (*tbis); pole++; tbis++;
1271         *pole += coef * (*tbis);
1272         pole  -= 3;
1273         tbis  -= 3;
1274       }
1275     }
1276     // multiply by the degrees
1277     coef = Degree;
1278     Dmi  = Degree;
1279     pole = psav + 2;
1280     
1281     for (i = 1; i <= min; i++) {
1282       *pole *= coef; pole++;
1283       *pole *= coef; pole++;
1284       Dmi--;
1285       coef  *= Dmi;
1286     }
1287     break;
1288   }
1289   case 3 : {
1290     psDD     = psav + (Degree << 1) + Degree;
1291     psDDmDim = psDD - 3;
1292     
1293     for (i = 0; i < Degree; i++) {
1294       DDmi--;
1295       pole = psDD;
1296       tbis = psDDmDim;
1297       jDmi = DDmi;
1298       
1299       for (j = Degm1; j >= i; j--) {
1300         jDmi--;
1301         coef   = (knot[jDmi] == knot[j]) ? 0.0 : 1. / (knot[jDmi] - knot[j]);
1302         *pole -= *tbis; *pole *= coef; pole++; tbis++;
1303         *pole -= *tbis; *pole *= coef; pole++; tbis++;
1304         *pole -= *tbis; *pole *= coef;
1305         pole  -= 5;
1306         tbis  -= 5;
1307       }
1308     }
1309     // Second phase, dependant of U
1310     iDim = - 3;
1311     
1312     for (i = 0; i < Degree; i++) {
1313       iDim += 3;
1314       pole  = psav + iDim;
1315       tbis  = pole + 3;
1316       coef  = U - knot[i];
1317       
1318       for (j = i; j >= 0; j--) {
1319         *pole += coef * (*tbis); pole++; tbis++;
1320         *pole += coef * (*tbis); pole++; tbis++;
1321         *pole += coef * (*tbis);
1322         pole  -= 5;
1323         tbis  -= 5;
1324       }
1325     }
1326     // multiply by the degrees
1327     coef = Degree;
1328     Dmi  = Degree;
1329     pole = psav + 3;
1330     
1331     for (i = 1; i <= min; i++) {
1332       *pole *= coef; pole++;
1333       *pole *= coef; pole++;
1334       *pole *= coef; pole++;
1335       Dmi--;
1336       coef  *= Dmi;
1337     }
1338     break;
1339   }
1340   case 4 : {
1341     psDD     = psav + (Degree << 2);
1342     psDDmDim = psDD - 4;
1343     
1344     for (i = 0; i < Degree; i++) {
1345       DDmi--;
1346       pole = psDD;
1347       tbis = psDDmDim;
1348       jDmi = DDmi;
1349       
1350       for (j = Degm1; j >= i; j--) {
1351         jDmi--;
1352         coef   = (knot[jDmi]  == knot[j]) ? 0.0 : 1. /(knot[jDmi] - knot[j]) ;
1353         *pole -= *tbis; *pole *= coef; pole++; tbis++;
1354         *pole -= *tbis; *pole *= coef; pole++; tbis++;
1355         *pole -= *tbis; *pole *= coef; pole++; tbis++;
1356         *pole -= *tbis; *pole *= coef;
1357         pole  -= 7;
1358         tbis  -= 7;
1359       }
1360     }
1361     // Second phase, dependant of U
1362     iDim = - 4;
1363     
1364     for (i = 0; i < Degree; i++) {
1365       iDim += 4;
1366       pole  = psav + iDim;
1367       tbis  = pole + 4;
1368       coef  = U - knot[i];
1369       
1370       for (j = i; j >= 0; j--) {
1371         *pole += coef * (*tbis); pole++; tbis++;
1372         *pole += coef * (*tbis); pole++; tbis++;
1373         *pole += coef * (*tbis); pole++; tbis++;
1374         *pole += coef * (*tbis);
1375         pole  -= 7;
1376         tbis  -= 7;
1377       }
1378     }
1379     // multiply by the degrees
1380     coef = Degree; 
1381     Dmi  = Degree;
1382     pole = psav + 4;
1383    
1384     for (i = 1; i <= min; i++) {
1385       *pole *= coef; pole++;
1386       *pole *= coef; pole++;
1387       *pole *= coef; pole++;
1388       *pole *= coef; pole++;
1389       Dmi--;
1390       coef  *= Dmi;
1391     }
1392     break;
1393   }
1394     default : {
1395       Standard_Integer k;
1396       Standard_Integer Dim2 = Dimension << 1;
1397       psDD     = psav + Degree * Dimension;
1398       psDDmDim = psDD - Dimension;
1399       
1400       for (i = 0; i < Degree; i++) {
1401         DDmi--;
1402         pole = psDD;
1403         tbis = psDDmDim;
1404         jDmi = DDmi;
1405         
1406         for (j = Degm1; j >= i; j--) {
1407           jDmi--;
1408           coef = (knot[jDmi] == knot[j]) ? 0.0 : 1. / (knot[jDmi] - knot[j]);
1409           
1410           for (k = 0; k < Dimension; k++) {
1411             *pole -= *tbis; *pole *= coef; pole++; tbis++;
1412           }
1413           pole -= Dim2;
1414           tbis -= Dim2;
1415         }
1416       }
1417       // Second phase, dependant of U
1418       iDim = - Dimension;
1419       
1420       for (i = 0; i < Degree; i++) {
1421         iDim += Dimension;
1422         pole  = psav + iDim;
1423         tbis  = pole + Dimension;
1424         coef  = U - knot[i];
1425         
1426         for (j = i; j >= 0; j--) {
1427           
1428           for (k = 0; k < Dimension; k++) {
1429             *pole += coef * (*tbis); pole++; tbis++;
1430           }
1431           pole -= Dim2;
1432           tbis -= Dim2;
1433         }
1434       }
1435       // multiply by the degrees
1436       coef = Degree;
1437       Dmi  = Degree;
1438       pole = psav + Dimension;
1439       
1440       for (i = 1; i <= min; i++) {
1441         
1442         for (k = 0; k < Dimension; k++) {
1443           *pole *= coef; pole++;
1444         }
1445         Dmi--;
1446         coef *= Dmi;
1447       }
1448     }
1449   }
1450 }
1451 
1452 //=======================================================================
1453 //function : BuildKnots
1454 //purpose  : 
1455 //=======================================================================
1456 
1457 void BSplCLib::BuildKnots(const Standard_Integer         Degree,
1458                           const Standard_Integer         Index,
1459                           const Standard_Boolean         Periodic,
1460                           const TColStd_Array1OfReal&    Knots,
1461                           const TColStd_Array1OfInteger* Mults,
1462                           Standard_Real&                 LK)
1463 {
1464   Standard_Integer KLower = Knots.Lower();
1465   const Standard_Real * pkn = &Knots(KLower);
1466   pkn -= KLower;
1467   Standard_Real *knot = &LK;
1468   if (Mults == NULL) {
1469     switch (Degree) {
1470     case 1 : {
1471       Standard_Integer j = Index    ;
1472       knot[0] = pkn[j]; j++;
1473       knot[1] = pkn[j];
1474       break;
1475     }
1476     case 2 : {
1477       Standard_Integer j = Index - 1;
1478       knot[0] = pkn[j]; j++;
1479       knot[1] = pkn[j]; j++;
1480       knot[2] = pkn[j]; j++;
1481       knot[3] = pkn[j];
1482       break;
1483     }
1484     case 3 : {
1485       Standard_Integer j = Index - 2;
1486       knot[0] = pkn[j]; j++;
1487       knot[1] = pkn[j]; j++;
1488       knot[2] = pkn[j]; j++;
1489       knot[3] = pkn[j]; j++;
1490       knot[4] = pkn[j]; j++;
1491       knot[5] = pkn[j];
1492       break;
1493     }
1494     case 4 : {
1495       Standard_Integer j = Index - 3;
1496       knot[0] = pkn[j]; j++;
1497       knot[1] = pkn[j]; j++;
1498       knot[2] = pkn[j]; j++;
1499       knot[3] = pkn[j]; j++;
1500       knot[4] = pkn[j]; j++;
1501       knot[5] = pkn[j]; j++;
1502       knot[6] = pkn[j]; j++;
1503       knot[7] = pkn[j];
1504       break;
1505     }
1506     case 5 : {
1507       Standard_Integer j = Index - 4;
1508       knot[0] = pkn[j]; j++;
1509       knot[1] = pkn[j]; j++;
1510       knot[2] = pkn[j]; j++;
1511       knot[3] = pkn[j]; j++;
1512       knot[4] = pkn[j]; j++;
1513       knot[5] = pkn[j]; j++;
1514       knot[6] = pkn[j]; j++;
1515       knot[7] = pkn[j]; j++;
1516       knot[8] = pkn[j]; j++;
1517       knot[9] = pkn[j];
1518       break;
1519     }
1520     case 6 : {
1521       Standard_Integer j = Index - 5;
1522       knot[ 0] = pkn[j]; j++;
1523       knot[ 1] = pkn[j]; j++;
1524       knot[ 2] = pkn[j]; j++;
1525       knot[ 3] = pkn[j]; j++;
1526       knot[ 4] = pkn[j]; j++;
1527       knot[ 5] = pkn[j]; j++;
1528       knot[ 6] = pkn[j]; j++;
1529       knot[ 7] = pkn[j]; j++;
1530       knot[ 8] = pkn[j]; j++;
1531       knot[ 9] = pkn[j]; j++;
1532       knot[10] = pkn[j]; j++;
1533       knot[11] = pkn[j];
1534       break;
1535     }
1536       default : {
1537         Standard_Integer i,j;
1538         Standard_Integer Deg2 = Degree << 1;
1539         j = Index - Degree;
1540         
1541         for (i = 0; i < Deg2; i++) {
1542           j++;
1543           knot[i] = pkn[j];
1544         }
1545       }
1546     }
1547   }
1548   else {
1549     Standard_Integer i;
1550     Standard_Integer Deg1 = Degree - 1;
1551     Standard_Integer KUpper = Knots.Upper();
1552     Standard_Integer MLower = Mults->Lower();
1553     Standard_Integer MUpper = Mults->Upper();
1554     const Standard_Integer * pmu = &(*Mults)(MLower);
1555     pmu -= MLower;
1556     Standard_Real dknot = 0;
1557     Standard_Integer ilow = Index    , mlow = 0;
1558     Standard_Integer iupp = Index + 1, mupp = 0;
1559     Standard_Real loffset = 0., uoffset = 0.;
1560     Standard_Boolean getlow = Standard_True, getupp = Standard_True;
1561     if (Periodic) {
1562       dknot = pkn[KUpper] - pkn[KLower];
1563       if (iupp > MUpper) {
1564         iupp = MLower + 1;
1565         uoffset = dknot;
1566       }
1567     }
1568     // Find the knots around Index
1569 
1570     for (i = 0; i < Degree; i++) {
1571       if (getlow) {
1572         mlow++;
1573         if (mlow > pmu[ilow]) {
1574           mlow = 1;
1575           ilow--;
1576           getlow =  (ilow >= MLower);
1577           if (Periodic && !getlow) {
1578             ilow = MUpper - 1;
1579             loffset = dknot;
1580             getlow = Standard_True;
1581           }
1582         }
1583         if (getlow)
1584           knot[Deg1 - i] = pkn[ilow] - loffset;
1585       }
1586       if (getupp) {
1587         mupp++;
1588         if (mupp > pmu[iupp]) {
1589           mupp = 1;
1590           iupp++;
1591           getupp = (iupp <= MUpper);
1592           if (Periodic && !getupp) {
1593             iupp = MLower + 1;
1594             uoffset = dknot;
1595             getupp = Standard_True;
1596           }
1597         }
1598         if (getupp)
1599           knot[Degree + i] = pkn[iupp] + uoffset;
1600       }
1601     }
1602   } 
1603 }
1604 
1605 //=======================================================================
1606 //function : PoleIndex
1607 //purpose  : 
1608 //=======================================================================
1609 
1610 Standard_Integer BSplCLib::PoleIndex(const Standard_Integer         Degree,
1611                                      const Standard_Integer         Index,
1612                                      const Standard_Boolean         Periodic,
1613                                      const TColStd_Array1OfInteger& Mults)
1614 {
1615   Standard_Integer i, pindex = 0;
1616 
1617   for (i = Mults.Lower(); i <= Index; i++)
1618     pindex += Mults(i);
1619   if (Periodic)
1620     pindex -= Mults(Mults.Lower());
1621   else
1622     pindex -= Degree + 1;
1623 
1624   return pindex;
1625 }
1626 
1627 //=======================================================================
1628 //function : BuildBoor
1629 //purpose  : builds the local array for boor
1630 //=======================================================================
1631 
1632 void  BSplCLib::BuildBoor(const Standard_Integer         Index,
1633                           const Standard_Integer         Length,
1634                           const Standard_Integer         Dimension,
1635                           const TColStd_Array1OfReal&    Poles,
1636                           Standard_Real&                 LP)
1637 {
1638   Standard_Real *poles = &LP;
1639   Standard_Integer i,k, ip = Poles.Lower() + Index * Dimension;
1640   
1641   for (i = 0; i < Length+1; i++) {
1642 
1643     for (k = 0; k < Dimension; k++) {
1644       poles[k] = Poles(ip);
1645       ip++;
1646       if (ip > Poles.Upper()) ip = Poles.Lower();
1647     }
1648     poles += 2 * Dimension;
1649   }
1650 }
1651 
1652 //=======================================================================
1653 //function : BoorIndex
1654 //purpose  : 
1655 //=======================================================================
1656 
1657 Standard_Integer  BSplCLib::BoorIndex(const Standard_Integer Index,
1658                                       const Standard_Integer Length,
1659                                       const Standard_Integer Depth)
1660 {
1661   if (Index <= Depth)  return Index;
1662   if (Index <= Length) return 2 * Index - Depth;
1663   return                      Length + Index - Depth;
1664 }
1665 
1666 //=======================================================================
1667 //function : GetPole
1668 //purpose  : 
1669 //=======================================================================
1670 
1671 void  BSplCLib::GetPole(const Standard_Integer         Index,
1672                         const Standard_Integer         Length,
1673                         const Standard_Integer         Depth,
1674                         const Standard_Integer         Dimension,
1675                         Standard_Real&                 LP,
1676                         Standard_Integer&              Position,
1677                         TColStd_Array1OfReal&          Pole)
1678 {
1679   Standard_Integer k;
1680   Standard_Real* pole = &LP + BoorIndex(Index,Length,Depth) * Dimension;
1681 
1682   for (k = 0; k < Dimension; k++) {
1683     Pole(Position) = pole[k];
1684     Position++;
1685   }
1686   if (Position > Pole.Upper()) Position = Pole.Lower();
1687 }
1688 
1689 //=======================================================================
1690 //function : PrepareInsertKnots
1691 //purpose  : 
1692 //=======================================================================
1693 
1694 Standard_Boolean  BSplCLib::PrepareInsertKnots
1695 (const Standard_Integer         Degree,
1696  const Standard_Boolean         Periodic, 
1697  const TColStd_Array1OfReal&    Knots,
1698  const TColStd_Array1OfInteger& Mults,
1699  const TColStd_Array1OfReal&    AddKnots,
1700  const TColStd_Array1OfInteger* AddMults,
1701  Standard_Integer&              NbPoles,
1702  Standard_Integer&              NbKnots, 
1703  const Standard_Real            Tolerance,
1704  const Standard_Boolean         Add)
1705 {
1706   Standard_Boolean addflat = AddMults == NULL;
1707   
1708   Standard_Integer first,last;
1709   if (Periodic) {
1710     first = Knots.Lower();
1711     last  = Knots.Upper();
1712   }
1713   else {
1714     first = FirstUKnotIndex(Degree,Mults);
1715     last  = LastUKnotIndex(Degree,Mults);
1716   }
1717   Standard_Real adeltaK1 = Knots(first)-AddKnots(AddKnots.Lower());
1718    Standard_Real adeltaK2 = AddKnots(AddKnots.Upper())-Knots(last);
1719   if (adeltaK1 > Tolerance) return Standard_False;
1720   if (adeltaK2  > Tolerance) return Standard_False;
1721   
1722   Standard_Integer sigma = 0, mult, amult;
1723   NbKnots = 0;
1724   Standard_Integer  k  = Knots.Lower() - 1;
1725   Standard_Integer  ak = AddKnots.Lower();
1726 
1727   if(Periodic && AddKnots.Length() > 1)
1728   {
1729     //gka for case when segments was produced on full period only one knot
1730     //was added in the end of curve
1731     if(fabs(adeltaK1) <= gp::Resolution() && 
1732        fabs(adeltaK2) <= gp::Resolution())
1733       ak++;
1734   }
1735   
1736   Standard_Integer aLastKnotMult = Mults (Knots.Upper());
1737   Standard_Real au,oldau = AddKnots(ak),Eps;
1738   
1739   while (ak <= AddKnots.Upper()) {
1740     au = AddKnots(ak);
1741     if (au < oldau) return Standard_False;
1742     oldau = au;
1743 
1744     Eps = Max(Tolerance,Epsilon(au));
1745     
1746     while ((k < Knots.Upper()) && (Knots(k+1)  - au <= Eps)) {
1747       k++;
1748       NbKnots++;
1749       sigma += Mults(k);
1750     }
1751 
1752     if (addflat) amult = 1;
1753     else         amult = Max(0,(*AddMults)(ak));
1754     
1755     while ((ak < AddKnots.Upper()) &&
1756            (Abs(au - AddKnots(ak+1)) <= Eps)) {
1757       ak++;
1758       if (Add) {
1759         if (addflat) amult++;
1760         else         amult += Max(0,(*AddMults)(ak));
1761       }
1762     }
1763     
1764     
1765     if (Abs(au - Knots(k)) <= Eps) {
1766       // identic to existing knot
1767       mult = Mults(k);
1768       if (Add) {
1769         if (mult + amult > Degree)
1770           amult = Max(0,Degree - mult);
1771         sigma += amult;
1772         
1773       }
1774       else if (amult > mult) {
1775         if (amult > Degree) amult = Degree;
1776         if (k == Knots.Upper () && Periodic)
1777         {
1778           aLastKnotMult = Max (amult, mult);
1779           sigma += 2 * (aLastKnotMult - mult);
1780         }
1781         else
1782         {
1783           sigma += amult - mult;
1784         }
1785       }
1786       /*
1787       // on periodic curves if this is the last knot
1788       // the multiplicity is added twice to account for the first knot
1789       if (k == Knots.Upper() && Periodic) {
1790         if (Add)
1791           sigma += amult;
1792         else
1793           sigma += amult - mult;
1794       }
1795       */
1796     }
1797     else {
1798       // not identic to existing knot
1799       if (amult > 0) {
1800         if (amult > Degree) amult = Degree;
1801         NbKnots++;
1802         sigma += amult;
1803       }
1804     }
1805     
1806     ak++;
1807   }
1808   
1809   // count the last knots
1810   while (k < Knots.Upper()) {
1811     k++;
1812     NbKnots++;
1813     sigma += Mults(k);
1814   }
1815 
1816   if (Periodic) {
1817     //for periodic B-Spline the requirement is that multiplicities of the first
1818     //and last knots must be equal (see Geom_BSplineCurve constructor for
1819     //instance);
1820     //respectively AddMults() must meet this requirement if AddKnots() contains
1821     //knot(s) coincident with first or last
1822     NbPoles = sigma - aLastKnotMult;
1823   }
1824   else {
1825     NbPoles = sigma - Degree - 1;
1826   }
1827  
1828   return Standard_True;
1829 }
1830 
1831 //=======================================================================
1832 //function : Copy
1833 //purpose  : copy reals from an array to an other
1834 //        
1835 //   NbValues are copied from OldPoles(OldFirst)
1836 //                 to    NewPoles(NewFirst)
1837 //
1838 //   Periodicity is handled.
1839 //   OldFirst and NewFirst are updated 
1840 //   to the position after the last copied pole.
1841 //
1842 //=======================================================================
1843 
1844 static void Copy(const Standard_Integer      NbPoles,
1845                  Standard_Integer&           OldFirst,
1846                  const TColStd_Array1OfReal& OldPoles,
1847                  Standard_Integer&           NewFirst,
1848                  TColStd_Array1OfReal&       NewPoles)
1849 {
1850   // reset the index in the range for periodicity
1851 
1852   OldFirst = OldPoles.Lower() + 
1853     (OldFirst - OldPoles.Lower()) % (OldPoles.Upper() - OldPoles.Lower() + 1);
1854 
1855   NewFirst = NewPoles.Lower() + 
1856     (NewFirst - NewPoles.Lower()) % (NewPoles.Upper() - NewPoles.Lower() + 1);
1857 
1858   // copy
1859   Standard_Integer i;
1860 
1861   for (i = 1; i <= NbPoles; i++) {
1862     NewPoles(NewFirst) = OldPoles(OldFirst);
1863     OldFirst++;
1864     if (OldFirst > OldPoles.Upper()) OldFirst = OldPoles.Lower();
1865     NewFirst++;
1866     if (NewFirst > NewPoles.Upper()) NewFirst = NewPoles.Lower();
1867   }
1868 }
1869                       
1870 //=======================================================================
1871 //function : InsertKnots
1872 //purpose  : insert an array of knots and multiplicities
1873 //=======================================================================
1874 
1875 void BSplCLib::InsertKnots
1876 (const Standard_Integer         Degree, 
1877  const Standard_Boolean         Periodic,
1878  const Standard_Integer         Dimension, 
1879  const TColStd_Array1OfReal&    Poles,  
1880  const TColStd_Array1OfReal&    Knots,    
1881  const TColStd_Array1OfInteger& Mults, 
1882  const TColStd_Array1OfReal&    AddKnots,    
1883  const TColStd_Array1OfInteger* AddMults, 
1884  TColStd_Array1OfReal&          NewPoles,
1885  TColStd_Array1OfReal&          NewKnots,    
1886  TColStd_Array1OfInteger&       NewMults, 
1887  const Standard_Real            Tolerance,
1888  const Standard_Boolean         Add)
1889 {
1890   Standard_Boolean addflat  = AddMults == NULL;
1891   
1892   Standard_Integer i,k,mult,firstmult;
1893   Standard_Integer index,kn,curnk,curk;
1894   Standard_Integer p,np, curp, curnp, length, depth;
1895   Standard_Real u;
1896   Standard_Integer need;
1897   Standard_Real Eps;
1898 
1899   // -------------------
1900   // create local arrays
1901   // -------------------
1902 
1903   Standard_Real *knots = new Standard_Real[2*Degree];
1904   Standard_Real *poles = new Standard_Real[(2*Degree+1)*Dimension];
1905   
1906   //----------------------------
1907   // loop on the knots to insert
1908   //----------------------------
1909   
1910   curk   = Knots.Lower()-1;          // current position in Knots
1911   curnk  = NewKnots.Lower()-1;       // current position in NewKnots
1912   curp   = Poles.Lower();            // current position in Poles
1913   curnp  = NewPoles.Lower();         // current position in NewPoles
1914 
1915   // NewKnots, NewMults, NewPoles contains the current state of the curve
1916 
1917   // index is the first pole of the current curve for insertion schema
1918 
1919   if (Periodic) index = -Mults(Mults.Lower());
1920   else          index = -Degree-1;
1921 
1922   // on Periodic curves the first knot and the last knot are inserted later
1923   // (they are the same knot)
1924   firstmult = 0;  // multiplicity of the first-last knot for periodic
1925   
1926 
1927   // kn current knot to insert in AddKnots
1928 
1929   for (kn = AddKnots.Lower(); kn <= AddKnots.Upper(); kn++) {
1930     
1931     u = AddKnots(kn);
1932     Eps = Max(Tolerance,Epsilon(u));
1933     
1934     //-----------------------------------
1935     // find the position in the old knots
1936     // and copy to the new knots
1937     //-----------------------------------
1938     
1939     while (curk < Knots.Upper() && Knots(curk+1) - u <= Eps) {
1940       curk++; curnk++;
1941       NewKnots(curnk) = Knots(curk);
1942       index += NewMults(curnk) = Mults(curk);
1943     }
1944     
1945     //-----------------------------------
1946     // Slice the knots and the mults
1947     // to the current size of the new curve
1948     //-----------------------------------
1949 
1950     i = curnk + Knots.Upper() - curk;
1951     TColStd_Array1OfReal    nknots(NewKnots(NewKnots.Lower()),NewKnots.Lower(),i);
1952     TColStd_Array1OfInteger nmults(NewMults(NewMults.Lower()),NewMults.Lower(),i);
1953 
1954     //-----------------------------------
1955     // copy enough knots 
1956     // to compute the insertion schema
1957     //-----------------------------------
1958 
1959     k = curk;
1960     i = curnk;
1961     mult = 0;
1962 
1963     while (mult < Degree && k < Knots.Upper()) {
1964       k++; i++;
1965       nknots(i) = Knots(k);
1966       mult += nmults(i) = Mults(k);
1967     }
1968 
1969     // copy knots at the end for periodic curve
1970     if (Periodic) {
1971       mult = 0;
1972       k = Knots.Upper();
1973       i = nknots.Upper();
1974 
1975       while (mult < Degree && i > curnk) {
1976         nknots(i) = Knots(k);
1977         mult += nmults(i) = Mults(k);
1978         k--;
1979         i--;
1980       }
1981       nmults(nmults.Upper()) = nmults(nmults.Lower());
1982     }
1983 
1984   
1985 
1986     //------------------------------------
1987     // do the boor scheme on the new curve
1988     // to insert the new knot
1989     //------------------------------------
1990     
1991     Standard_Boolean sameknot = (Abs(u-NewKnots(curnk)) <= Eps);
1992     
1993     if (sameknot) length = Max(0,Degree - NewMults(curnk));
1994     else          length = Degree;
1995     
1996     if (addflat) depth = 1;
1997     else         depth = Min(Degree,(*AddMults)(kn));
1998 
1999     if (sameknot) {
2000       if (Add) {
2001         if ((NewMults(curnk) + depth) > Degree)
2002           depth = Degree - NewMults(curnk);
2003       }
2004       else {
2005         depth = Max(0,depth-NewMults(curnk));
2006       }
2007 
2008       if (Periodic) {
2009         // on periodic curve the first and last knot are delayed to the end
2010         if (curk == Knots.Lower() || (curk == Knots.Upper())) {
2011           if (firstmult == 0) // do that only once
2012             firstmult += depth;
2013           depth = 0;
2014         }
2015       }
2016     }
2017     if (depth <= 0) continue;
2018     
2019     BuildKnots(Degree,curnk,Periodic,nknots,&nmults,*knots);
2020 
2021     // copy the poles
2022 
2023     need   = NewPoles.Lower()+(index+length+1)*Dimension - curnp;
2024     need = Min(need,Poles.Upper() - curp + 1);
2025 
2026     p = curp;
2027     np = curnp;
2028     Copy(need,p,Poles,np,NewPoles);
2029     curp  += need;
2030     curnp += need;
2031 
2032     // slice the poles to the current number of poles in case of periodic
2033     TColStd_Array1OfReal npoles(NewPoles(NewPoles.Lower()),NewPoles.Lower(),curnp-1);
2034 
2035     BuildBoor(index,length,Dimension,npoles,*poles);
2036     BoorScheme(u,Degree,*knots,Dimension,*poles,depth,length);
2037     
2038     //-------------------
2039     // copy the new poles
2040     //-------------------
2041 
2042     curnp += depth * Dimension; // number of poles is increased by depth
2043     TColStd_Array1OfReal ThePoles(NewPoles(NewPoles.Lower()),NewPoles.Lower(),curnp-1);
2044     np = NewKnots.Lower()+(index+1)*Dimension;
2045 
2046     for (i = 1; i <= length + depth; i++)
2047       GetPole(i,length,depth,Dimension,*poles,np,ThePoles);
2048     
2049     //-------------------
2050     // insert the knot
2051     //-------------------
2052 
2053     index += depth;
2054     if (sameknot) {
2055       NewMults(curnk) += depth;
2056     }
2057     else {
2058       curnk++;
2059       NewKnots(curnk) = u;
2060       NewMults(curnk) = depth;
2061     }
2062   }
2063   
2064   //------------------------------
2065   // copy the last poles and knots
2066   //------------------------------
2067   
2068   Copy(Poles.Upper() - curp + 1,curp,Poles,curnp,NewPoles);
2069   
2070   while (curk < Knots.Upper()) {
2071     curk++;  curnk++;
2072     NewKnots(curnk) = Knots(curk);
2073     NewMults(curnk) = Mults(curk);
2074   }
2075   
2076   //------------------------------
2077   // process the first-last knot 
2078   // on periodic curves
2079   //------------------------------
2080 
2081   if (firstmult > 0) {
2082     curnk = NewKnots.Lower();
2083     if (NewMults(curnk) + firstmult > Degree) {
2084       firstmult = Degree - NewMults(curnk);
2085     }
2086     if (firstmult > 0) {
2087 
2088       length = Degree - NewMults(curnk);
2089       depth  = firstmult;
2090 
2091       BuildKnots(Degree,curnk,Periodic,NewKnots,&NewMults,*knots);
2092       TColStd_Array1OfReal npoles(NewPoles(NewPoles.Lower()),
2093                                   NewPoles.Lower(),
2094                                   NewPoles.Upper()-depth*Dimension);
2095       BuildBoor(0,length,Dimension,npoles,*poles);
2096       BoorScheme(NewKnots(curnk),Degree,*knots,Dimension,*poles,depth,length);
2097       
2098       //---------------------------
2099       // copy the new poles
2100       // but rotate them with depth
2101       //---------------------------
2102       
2103       np = NewPoles.Lower();
2104 
2105       for (i = depth; i < length + depth; i++)
2106         GetPole(i,length,depth,Dimension,*poles,np,NewPoles);
2107 
2108       np = NewPoles.Upper() - depth*Dimension + 1;
2109 
2110       for (i = 0; i < depth; i++)
2111         GetPole(i,length,depth,Dimension,*poles,np,NewPoles);
2112       
2113       NewMults(NewMults.Lower()) += depth;
2114       NewMults(NewMults.Upper()) += depth;
2115     }
2116   }
2117   // free local arrays
2118   delete [] knots;
2119   delete [] poles;
2120 }
2121 
2122 //=======================================================================
2123 //function : RemoveKnot
2124 //purpose  : 
2125 //=======================================================================
2126 
2127 Standard_Boolean BSplCLib::RemoveKnot 
2128 (const Standard_Integer         Index,       
2129  const Standard_Integer         Mult,        
2130  const Standard_Integer         Degree,  
2131  const Standard_Boolean         Periodic,
2132  const Standard_Integer         Dimension,  
2133  const TColStd_Array1OfReal&    Poles,
2134  const TColStd_Array1OfReal&    Knots,  
2135  const TColStd_Array1OfInteger& Mults,
2136  TColStd_Array1OfReal&          NewPoles,
2137  TColStd_Array1OfReal&          NewKnots,  
2138  TColStd_Array1OfInteger&       NewMults,
2139  const Standard_Real            Tolerance)
2140 {
2141   Standard_Integer index,i,j,k,p,np;
2142 
2143   Standard_Integer TheIndex = Index;
2144 
2145   // protection
2146   Standard_Integer first,last;
2147   if (Periodic) {
2148     first = Knots.Lower();
2149     last  = Knots.Upper();
2150   }
2151   else {
2152     first = BSplCLib::FirstUKnotIndex(Degree,Mults) + 1;
2153     last  = BSplCLib::LastUKnotIndex(Degree,Mults) - 1;
2154   }
2155   if (Index < first) return Standard_False;
2156   if (Index > last)  return Standard_False;
2157 
2158   if ( Periodic && (Index == first)) TheIndex = last;
2159 
2160   Standard_Integer depth  = Mults(TheIndex) - Mult;
2161   Standard_Integer length = Degree - Mult;
2162 
2163   // -------------------
2164   // create local arrays
2165   // -------------------
2166 
2167   Standard_Real *knots = new Standard_Real[4*Degree];
2168   Standard_Real *poles = new Standard_Real[(2*Degree+1)*Dimension];
2169   
2170 
2171   // ------------------------------------
2172   // build the knots for anti Boor Scheme
2173   // ------------------------------------
2174 
2175   // the new sequence of knots
2176   // is obtained from the knots at Index-1 and Index
2177   
2178   BSplCLib::BuildKnots(Degree,TheIndex-1,Periodic,Knots,&Mults,*knots);
2179   index = PoleIndex(Degree,TheIndex-1,Periodic,Mults);
2180   BSplCLib::BuildKnots(Degree,TheIndex,Periodic,Knots,&Mults,knots[2*Degree]);
2181 
2182   index += Mult;
2183 
2184   for (i = 0; i < Degree - Mult; i++)
2185     knots[i] = knots[i+Mult];
2186 
2187   for (i = Degree-Mult; i < 2*Degree; i++)
2188     knots[i] = knots[2*Degree+i];
2189 
2190 
2191   // ------------------------------------
2192   // build the poles for anti Boor Scheme
2193   // ------------------------------------
2194 
2195   p = Poles.Lower()+index * Dimension;
2196 
2197   for (i = 0; i <= length + depth; i++) {
2198     j = Dimension * BoorIndex(i,length,depth);
2199 
2200     for (k = 0; k < Dimension; k++) {
2201       poles[j+k] = Poles(p+k);
2202     }
2203     p += Dimension;
2204     if (p > Poles.Upper()) p = Poles.Lower();
2205   }
2206 
2207 
2208   // ----------------
2209   // Anti Boor Scheme
2210   // ----------------
2211 
2212   Standard_Boolean result = AntiBoorScheme(Knots(TheIndex),Degree,*knots,
2213                                            Dimension,*poles,
2214                                            depth,length,Tolerance);
2215   
2216   // ----------------
2217   // copy the results
2218   // ----------------
2219 
2220   if (result) {
2221 
2222     // poles
2223 
2224     p = Poles.Lower();
2225     np = NewPoles.Lower();
2226     
2227     // unmodified poles before
2228     Copy((index+1)*Dimension,p,Poles,np,NewPoles);
2229     
2230     
2231     // modified
2232 
2233     for (i = 1; i <= length; i++)
2234       BSplCLib::GetPole(i,length,0,Dimension,*poles,np,NewPoles);
2235     p += (length + depth) * Dimension ;
2236     
2237     // unmodified poles after
2238     if (p != Poles.Lower()) {
2239       i = Poles.Upper() - p + 1;
2240       Copy(i,p,Poles,np,NewPoles);
2241     }
2242 
2243     // knots and mults
2244 
2245     if (Mult > 0) {
2246       NewKnots = Knots;
2247       NewMults = Mults;
2248       NewMults(TheIndex) = Mult;
2249       if (Periodic) {
2250         if (TheIndex == first) NewMults(last)  = Mult;
2251         if (TheIndex == last)  NewMults(first) = Mult;
2252       }
2253     }
2254     else {
2255       if (!Periodic || (TheIndex != first && TheIndex != last)) {
2256 
2257         for (i = Knots.Lower(); i < TheIndex; i++) {
2258           NewKnots(i) = Knots(i);
2259           NewMults(i) = Mults(i);
2260         }
2261 
2262         for (i = TheIndex+1; i <= Knots.Upper(); i++) {
2263           NewKnots(i-1) = Knots(i);
2264           NewMults(i-1) = Mults(i);
2265         }
2266       }
2267       else {
2268         // The interesting case of a Periodic curve 
2269         // where the first and last knot is removed.
2270         
2271         for (i = first; i < last-1; i++) {
2272           NewKnots(i) = Knots(i+1);
2273           NewMults(i) = Mults(i+1);
2274         }
2275         NewKnots(last-1) = NewKnots(first) + Knots(last) - Knots(first);
2276         NewMults(last-1) = NewMults(first);
2277       }
2278     }
2279   }
2280 
2281 
2282   // free local arrays
2283   delete [] knots;
2284   delete [] poles;
2285   
2286   return result;
2287 }
2288 
2289 //=======================================================================
2290 //function : IncreaseDegreeCountKnots
2291 //purpose  : 
2292 //=======================================================================
2293 
2294 Standard_Integer  BSplCLib::IncreaseDegreeCountKnots
2295 (const Standard_Integer Degree,
2296  const Standard_Integer NewDegree, 
2297  const Standard_Boolean Periodic, 
2298  const TColStd_Array1OfInteger& Mults)
2299 {
2300   if (Periodic) return Mults.Length();
2301   Standard_Integer f = FirstUKnotIndex(Degree,Mults);
2302   Standard_Integer l = LastUKnotIndex(Degree,Mults);
2303   Standard_Integer m,i,removed = 0, step = NewDegree - Degree;
2304   
2305   i = Mults.Lower();
2306   m = Degree + (f - i + 1) * step + 1;
2307 
2308   while (m > NewDegree+1) {
2309     removed++;
2310     m -= Mults(i) + step;
2311     i++;
2312   }
2313   if (m < NewDegree+1) removed--;
2314 
2315   i = Mults.Upper();
2316   m = Degree + (i - l + 1) * step + 1;
2317 
2318   while (m > NewDegree+1) {
2319     removed++;
2320     m -= Mults(i) + step;
2321     i--;
2322   }
2323   if (m < NewDegree+1) removed--;
2324 
2325   return Mults.Length() - removed;
2326 }
2327 
2328 //=======================================================================
2329 //function : IncreaseDegree
2330 //purpose  : 
2331 //=======================================================================
2332 
2333 void BSplCLib::IncreaseDegree 
2334 (const Standard_Integer         Degree,
2335  const Standard_Integer         NewDegree,
2336  const Standard_Boolean         Periodic,
2337  const Standard_Integer         Dimension,
2338  const TColStd_Array1OfReal&    Poles,
2339  const TColStd_Array1OfReal&    Knots,
2340  const TColStd_Array1OfInteger& Mults,
2341  TColStd_Array1OfReal&          NewPoles,
2342  TColStd_Array1OfReal&          NewKnots,
2343  TColStd_Array1OfInteger&       NewMults)
2344 { 
2345   // Degree elevation of a BSpline Curve
2346 
2347   // This algorithms loops on degree incrementation from Degree to NewDegree.
2348   // The variable curDeg is the current degree to increment.
2349 
2350   // Before degree incrementations a "working curve" is created.
2351   // It has the same knot, poles and multiplicities.
2352 
2353   // If the curve is periodic knots are added on the working curve before
2354   // and after the existing knots to make it a non-periodic curves. 
2355   // The poles are also copied.
2356 
2357   // The first and last multiplicity of the working curve are set to Degree+1,
2358   // null poles are  added if necessary.
2359 
2360   // Then the degree is incremented on the working curve.
2361   // The knots are unchanged but all multiplicities will be incremented.
2362 
2363   // Each degree incrementation is achieved by averaging curDeg+1 curves.
2364 
2365   // See : Degree elevation of B-spline curves
2366   //       Hartmut PRAUTZSCH
2367   //       CAGD 1 (1984)
2368 
2369 
2370   //-------------------------
2371   // create the working curve
2372   //-------------------------
2373 
2374   Standard_Integer i,k,f,l,m,pf,pl,firstknot;
2375 
2376   pf = 0; // number of null poles added at beginning
2377   pl = 0; // number of null poles added at end
2378 
2379   Standard_Integer nbwknots = Knots.Length();
2380   f = FirstUKnotIndex(Degree,Mults);
2381   l = LastUKnotIndex (Degree,Mults);
2382 
2383   if (Periodic) {
2384     // Periodic curves are transformed in non-periodic curves
2385 
2386     nbwknots += f - Mults.Lower();
2387 
2388     pf = -Degree - 1;
2389 
2390     for (i = Mults.Lower(); i <= f; i++)
2391       pf += Mults(i);
2392 
2393     nbwknots += Mults.Upper() - l;
2394 
2395     pl = -Degree - 1;
2396 
2397     for (i = l; i <= Mults.Upper(); i++)
2398       pl += Mults(i);
2399   }
2400 
2401   // copy the knots and multiplicities 
2402   TColStd_Array1OfReal    wknots(1,nbwknots);
2403   TColStd_Array1OfInteger wmults(1,nbwknots);
2404   if (!Periodic) {
2405     wknots  = Knots;
2406     wmults  = Mults;
2407   }
2408   else {
2409     // copy the knots for a periodic curve
2410     Standard_Real period = Knots(Knots.Upper()) - Knots(Knots.Lower());
2411     i = 0;
2412 
2413     for (k = l; k < Knots.Upper(); k++) {
2414       i++; 
2415       wknots(i) = Knots(k) - period;
2416       wmults(i) = Mults(k);
2417     }
2418 
2419     for (k = Knots.Lower(); k <= Knots.Upper(); k++) {
2420       i++; 
2421       wknots(i) = Knots(k);
2422       wmults(i) = Mults(k);
2423     }
2424 
2425     for (k = Knots.Lower()+1; k <= f; k++) {
2426       i++; 
2427       wknots(i) = Knots(k) + period;
2428       wmults(i) = Mults(k);
2429     }
2430   }
2431 
2432   // set the first and last mults to Degree+1
2433   // and add null poles
2434 
2435   pf += Degree + 1 - wmults(1);
2436   wmults(1) = Degree + 1;
2437   pl += Degree + 1 - wmults(nbwknots);
2438   wmults(nbwknots) = Degree + 1;
2439 
2440   //---------------------------
2441   // poles of the working curve
2442   //---------------------------
2443 
2444   Standard_Integer nbwpoles = 0;
2445 
2446   for (i = 1; i <= nbwknots; i++) nbwpoles += wmults(i);
2447   nbwpoles -= Degree + 1;
2448 
2449   // we provide space for degree elevation
2450   TColStd_Array1OfReal 
2451     wpoles(1,(nbwpoles + (nbwknots-1) * (NewDegree - Degree)) * Dimension);
2452 
2453   for (i = 1; i <= pf * Dimension; i++) 
2454     wpoles(i) = 0;
2455 
2456   k = Poles.Lower();
2457 
2458   for (i = pf * Dimension + 1; i <= (nbwpoles - pl) * Dimension; i++) {
2459     wpoles(i) = Poles(k);
2460     k++;
2461     if (k > Poles.Upper()) k = Poles.Lower();
2462   }
2463 
2464   for (i = (nbwpoles-pl)*Dimension+1; i <= nbwpoles*Dimension; i++)
2465     wpoles(i) = 0;
2466   
2467   
2468   //-----------------------------------------------------------
2469   // Declare the temporary arrays used in degree incrementation
2470   //-----------------------------------------------------------
2471 
2472   Standard_Integer nbwp = nbwpoles + (nbwknots-1) * (NewDegree - Degree);
2473   // Arrays for storing the temporary curves
2474   TColStd_Array1OfReal tempc1(1,nbwp * Dimension);
2475   TColStd_Array1OfReal tempc2(1,nbwp * Dimension);
2476 
2477   // Array for storing the knots to insert
2478   TColStd_Array1OfReal iknots(1,nbwknots);
2479 
2480   // Arrays for receiving the knots after insertion
2481   TColStd_Array1OfReal    nknots(1,nbwknots);
2482 
2483 
2484   
2485   //------------------------------
2486   // Loop on degree incrementation
2487   //------------------------------
2488 
2489   Standard_Integer step,curDeg;
2490   Standard_Integer nbp = nbwpoles;
2491   nbwp = nbp;
2492 
2493   for (curDeg = Degree; curDeg < NewDegree; curDeg++) {
2494     
2495     nbp  = nbwp;               // current number of poles
2496     nbwp = nbp + nbwknots - 1; // new number of poles
2497 
2498     // For the averaging
2499     TColStd_Array1OfReal nwpoles(1,nbwp * Dimension);
2500     nwpoles.Init(0.0e0) ;
2501   
2502     
2503     for (step = 0; step <= curDeg; step++) {
2504     
2505       // Compute the bspline of rank step.
2506 
2507       // if not the first time, decrement the multiplicities back
2508       if (step != 0) {
2509         for (i = 1; i <= nbwknots; i++)
2510           wmults(i)--;
2511       }
2512     
2513       // Poles are the current poles 
2514       // but the poles congruent to step are duplicated.
2515       
2516       Standard_Integer offset = 0;
2517 
2518       for (i = 0; i < nbp; i++) {
2519         offset++;
2520 
2521         for (k = 0; k < Dimension; k++) {
2522           tempc1((offset-1)*Dimension+k+1) = 
2523             wpoles(NewPoles.Lower()+i*Dimension+k);
2524         }
2525         if (i % (curDeg+1) == step) {
2526           offset++;
2527 
2528           for (k = 0; k < Dimension; k++) {
2529             tempc1((offset-1)*Dimension+k+1) = 
2530               wpoles(NewPoles.Lower()+i*Dimension+k);
2531           }
2532         }
2533       }
2534         
2535       // Knots multiplicities are increased
2536       // For each knot where the sum of multiplicities is congruent to step
2537       
2538       Standard_Integer stepmult = step+1;
2539       Standard_Integer nbknots = 0;
2540       Standard_Integer smult = 0;
2541 
2542       for (k = 1; k <= nbwknots; k++) {
2543         smult += wmults(k);
2544         if (smult  >= stepmult) {
2545           // this knot is increased
2546           stepmult += curDeg+1;
2547           wmults(k)++;
2548         }
2549         else {
2550           // this knot is inserted
2551           nbknots++;
2552           iknots(nbknots) = wknots(k);
2553         }
2554       }
2555       
2556       // the curve is obtained by inserting the knots
2557       // to raise the multiplicities
2558 
2559       // we build "slices" of the arrays to set the correct size
2560       if (nbknots > 0) {
2561         TColStd_Array1OfReal aknots(iknots(1),1,nbknots);
2562         TColStd_Array1OfReal curve (tempc1(1),1,offset * Dimension);
2563         TColStd_Array1OfReal ncurve(tempc2(1),1,nbwp   * Dimension);
2564 //      InsertKnots(curDeg+1,Standard_False,Dimension,curve,wknots,wmults,
2565 //                  aknots,NoMults(),ncurve,nknots,wmults,Epsilon(1.));
2566 
2567         InsertKnots(curDeg+1,Standard_False,Dimension,curve,wknots,wmults,
2568                     aknots,NoMults(),ncurve,nknots,wmults,0.0);
2569         
2570         // add to the average
2571 
2572         for (i = 1; i <= nbwp * Dimension; i++)
2573           nwpoles(i) += ncurve(i);
2574       }
2575       else {
2576         // add to the average
2577 
2578         for (i = 1; i <= nbwp * Dimension; i++)
2579           nwpoles(i) += tempc1(i);
2580       }
2581     }
2582     
2583     // The result is the average
2584 
2585     for (i = 1; i <= nbwp * Dimension; i++) {
2586       wpoles(i) = nwpoles(i) / (curDeg+1);
2587     }
2588   }
2589   
2590   //-----------------
2591   // Copy the results
2592   //-----------------
2593 
2594   // index in new knots of the first knot of the curve
2595   if (Periodic)
2596     firstknot = Mults.Upper() - l + 1;
2597   else 
2598     firstknot = f;
2599   
2600   // the new curve starts at index firstknot
2601   // so we must remove knots until the sum of multiplicities
2602   // from the first to the start is NewDegree+1
2603 
2604   // m is the current sum of multiplicities
2605   m = 0;
2606 
2607   for (k = 1; k <= firstknot; k++)
2608     m += wmults(k);
2609 
2610   // compute the new first knot (k), pf will be the index of the first pole
2611   k = 1;
2612   pf = 0;
2613 
2614   while (m > NewDegree+1) {
2615     k++;
2616     m  -= wmults(k);
2617     pf += wmults(k);
2618   }
2619   if (m < NewDegree+1) {
2620     k--;
2621     wmults(k) += m - NewDegree - 1;
2622     pf        += m - NewDegree - 1;
2623   }
2624 
2625   // on a periodic curve the knots start with firstknot
2626   if (Periodic)
2627     k = firstknot;
2628 
2629   // copy knots
2630 
2631   for (i = NewKnots.Lower(); i <= NewKnots.Upper(); i++) {
2632     NewKnots(i) = wknots(k);
2633     NewMults(i) = wmults(k);
2634     k++;
2635   }
2636 
2637   // copy poles
2638   pf *= Dimension;
2639 
2640   for (i = NewPoles.Lower(); i <= NewPoles.Upper(); i++) {
2641     pf++;
2642     NewPoles(i) = wpoles(pf);
2643   }
2644 }
2645 
2646 //=======================================================================
2647 //function : PrepareUnperiodize
2648 //purpose  : 
2649 //=======================================================================
2650 
2651 void  BSplCLib::PrepareUnperiodize
2652 (const Standard_Integer         Degree, 
2653  const TColStd_Array1OfInteger& Mults, 
2654  Standard_Integer&        NbKnots, 
2655  Standard_Integer&        NbPoles)
2656 {
2657   Standard_Integer i;
2658   // initialize NbKnots and NbPoles
2659   NbKnots = Mults.Length();
2660   NbPoles = - Degree - 1;
2661 
2662   for (i = Mults.Lower(); i <= Mults.Upper(); i++) 
2663     NbPoles += Mults(i);
2664 
2665   Standard_Integer sigma, k;
2666   // Add knots at the beginning of the curve to raise Multiplicities 
2667   // to Degre + 1;
2668   sigma = Mults(Mults.Lower());
2669   k = Mults.Upper() - 1;
2670 
2671   while ( sigma < Degree + 1) {
2672     sigma   += Mults(k);
2673     NbPoles += Mults(k);
2674     k--;
2675     NbKnots++;
2676   }
2677   // We must add exactly until Degree + 1 -> 
2678   //    Suppress the excedent.
2679   if ( sigma > Degree + 1)
2680     NbPoles -= sigma - Degree - 1;
2681 
2682   // Add knots at the end of the curve to raise Multiplicities 
2683   // to Degre + 1;
2684   sigma = Mults(Mults.Upper());
2685   k = Mults.Lower() + 1;
2686 
2687   while ( sigma < Degree + 1) {
2688     sigma   += Mults(k);
2689     NbPoles += Mults(k);
2690     k++;
2691     NbKnots++;
2692   }
2693   // We must add exactly until Degree + 1 -> 
2694   //    Suppress the excedent.
2695   if ( sigma > Degree + 1)
2696     NbPoles -= sigma - Degree - 1;
2697 }
2698 
2699 //=======================================================================
2700 //function : Unperiodize
2701 //purpose  : 
2702 //=======================================================================
2703 
2704 void  BSplCLib::Unperiodize
2705 (const Standard_Integer         Degree,
2706  const Standard_Integer         , // Dimension,
2707  const TColStd_Array1OfInteger& Mults,
2708  const TColStd_Array1OfReal&    Knots,
2709  const TColStd_Array1OfReal&    Poles,
2710  TColStd_Array1OfInteger& NewMults,
2711  TColStd_Array1OfReal&    NewKnots,
2712  TColStd_Array1OfReal&    NewPoles)
2713 {
2714   Standard_Integer sigma, k, index = 0;
2715   // evaluation of index : number of knots to insert before knot(1) to
2716   // raise sum of multiplicities to <Degree + 1>
2717   sigma = Mults(Mults.Lower());
2718   k = Mults.Upper() - 1;
2719 
2720   while ( sigma < Degree + 1) {
2721     sigma   += Mults(k);
2722     k--;
2723     index++;
2724   }
2725 
2726   Standard_Real    period = Knots(Knots.Upper()) - Knots(Knots.Lower());
2727 
2728   // set the 'interior' knots;
2729 
2730   for ( k = 1; k <= Knots.Length(); k++) {
2731     NewKnots ( k + index ) = Knots( k);
2732     NewMults ( k + index ) = Mults( k);
2733   }
2734   
2735   // set the 'starting' knots;
2736 
2737   for ( k = 1; k <= index; k++) {
2738     NewKnots( k) = NewKnots( k + Knots.Length() - 1) - period;
2739     NewMults( k) = NewMults( k + Knots.Length() - 1);
2740   }
2741   NewMults( 1) -= sigma - Degree -1;
2742   
2743   // set the 'ending' knots;
2744   sigma = NewMults( index + Knots.Length() );
2745 
2746   for ( k = Knots.Length() + index + 1; k <= NewKnots.Length(); k++) {
2747     NewKnots( k) = NewKnots( k - Knots.Length() + 1) + period;
2748     NewMults( k) = NewMults( k - Knots.Length() + 1);
2749     sigma += NewMults( k - Knots.Length() + 1);
2750   }
2751   NewMults(NewMults.Length()) -= sigma - Degree - 1;
2752 
2753   for ( k = 1 ; k <= NewPoles.Length(); k++) {
2754     NewPoles(k ) = Poles( (k-1) % Poles.Length() + 1);
2755   }
2756 }
2757 
2758 //=======================================================================
2759 //function : PrepareTrimming
2760 //purpose  : 
2761 //=======================================================================
2762 
2763 void BSplCLib::PrepareTrimming(const Standard_Integer         Degree,
2764                                const Standard_Boolean         Periodic,
2765                                const TColStd_Array1OfReal&    Knots,
2766                                const TColStd_Array1OfInteger& Mults,
2767                                const Standard_Real            U1,
2768                                const Standard_Real            U2,
2769                                      Standard_Integer&        NbKnots,
2770                                      Standard_Integer&        NbPoles)
2771 {
2772   Standard_Integer i;
2773   Standard_Real NewU1, NewU2;
2774   Standard_Integer index1 = 0, index2 = 0;
2775 
2776   // Eval index1, index2 : position of U1 and U2 in the Array Knots
2777   // such as Knots(index1-1) <= U1 < Knots(index1)
2778   //         Knots(index2-1) <= U2 < Knots(index2)
2779   LocateParameter( Degree, Knots, Mults, U1, Periodic,
2780                    Knots.Lower(), Knots.Upper(), index1, NewU1);
2781   LocateParameter( Degree, Knots, Mults, U2, Periodic,
2782                    Knots.Lower(), Knots.Upper(), index2, NewU2);
2783   index1++;
2784   if ( Abs(Knots(index2) - U2) <= Epsilon( U1))
2785     index2--;
2786 
2787   // eval NbKnots:
2788   NbKnots = index2 - index1 + 3;
2789 
2790   // eval NbPoles:
2791   NbPoles = Degree + 1;
2792 
2793   for ( i = index1; i <= index2; i++) 
2794     NbPoles += Mults(i);
2795 }
2796 
2797 //=======================================================================
2798 //function : Trimming
2799 //purpose  : 
2800 //=======================================================================
2801 
2802 void BSplCLib::Trimming(const Standard_Integer         Degree,
2803                         const Standard_Boolean         Periodic,
2804                         const Standard_Integer         Dimension,
2805                         const TColStd_Array1OfReal&    Knots,
2806                         const TColStd_Array1OfInteger& Mults,
2807                         const TColStd_Array1OfReal&    Poles,
2808                         const Standard_Real            U1,
2809                         const Standard_Real            U2,
2810                               TColStd_Array1OfReal&    NewKnots,
2811                               TColStd_Array1OfInteger& NewMults,
2812                               TColStd_Array1OfReal&    NewPoles)
2813 {
2814   Standard_Integer i, nbpoles=0, nbknots=0;
2815   Standard_Real    kk[2] = { U1, U2 };
2816   Standard_Integer mm[2] = { Degree, Degree };
2817   TColStd_Array1OfReal    K( kk[0], 1, 2 );
2818   TColStd_Array1OfInteger M( mm[0], 1, 2 );
2819   if (!PrepareInsertKnots( Degree, Periodic, Knots, Mults, K, &M, 
2820                           nbpoles, nbknots, Epsilon( U1), 0))
2821   {
2822     throw Standard_OutOfRange();
2823   }
2824 
2825   TColStd_Array1OfReal    TempPoles(1, nbpoles*Dimension);
2826   TColStd_Array1OfReal    TempKnots(1, nbknots);
2827   TColStd_Array1OfInteger TempMults(1, nbknots);
2828 
2829 //
2830 // do not allow the multiplicities to Add : they must be less than Degree
2831 //
2832   InsertKnots(Degree, Periodic, Dimension, Poles, Knots, Mults,
2833               K, &M, TempPoles, TempKnots, TempMults, Epsilon(U1),
2834               Standard_False);
2835 
2836   // find in TempPoles the index of the pole corresponding to U1
2837   Standard_Integer Kindex = 0, Pindex;
2838   Standard_Real NewU1;
2839   LocateParameter( Degree, TempKnots, TempMults, U1, Periodic,
2840                    TempKnots.Lower(), TempKnots.Upper(), Kindex, NewU1);
2841   Pindex = PoleIndex ( Degree, Kindex, Periodic, TempMults);
2842   Pindex *= Dimension;
2843 
2844   for ( i = 1; i <= NewPoles.Length(); i++) NewPoles(i) = TempPoles(Pindex + i);
2845 
2846   for ( i = 1; i <= NewKnots.Length(); i++) {
2847     NewKnots(i) = TempKnots( Kindex+i-1);
2848     NewMults(i) = TempMults( Kindex+i-1);
2849   }
2850   NewMults(1) = Min(Degree, NewMults(1)) + 1 ;
2851   NewMults(NewMults.Length())= Min(Degree, NewMults(NewMults.Length())) + 1 ;
2852 }
2853 
2854 //=======================================================================
2855 //function : Solves a LU factored Matrix 
2856 //purpose  : 
2857 //=======================================================================
2858 
2859 Standard_Integer 
2860 BSplCLib::SolveBandedSystem(const math_Matrix&  Matrix,
2861                             const Standard_Integer UpperBandWidth,
2862                             const Standard_Integer LowerBandWidth,
2863                             const Standard_Integer ArrayDimension,
2864                             Standard_Real&   Array) 
2865 {
2866   Standard_Integer ii,
2867   jj,
2868   kk,
2869   MinIndex,
2870   MaxIndex,
2871   ReturnCode = 0 ;
2872   
2873   Standard_Real   *PolesArray = &Array ;
2874   Standard_Real   Inverse ;
2875   
2876   
2877   if (Matrix.LowerCol() != 1 || 
2878       Matrix.UpperCol() != UpperBandWidth + LowerBandWidth + 1) {
2879     ReturnCode = 1 ;
2880     goto FINISH ;
2881   }
2882   
2883   for (ii = Matrix.LowerRow() + 1; ii <=  Matrix.UpperRow() ; ii++) {
2884     MinIndex = (ii - LowerBandWidth >= Matrix.LowerRow() ?
2885                 ii - LowerBandWidth : Matrix.LowerRow()) ;
2886     
2887     for ( jj = MinIndex  ; jj < ii  ; jj++) {
2888       
2889       for (kk = 0 ; kk < ArrayDimension ; kk++) {
2890         PolesArray[(ii-1) * ArrayDimension + kk] += 
2891           PolesArray[(jj-1) * ArrayDimension + kk] * Matrix(ii, jj - ii + LowerBandWidth + 1) ;
2892       }
2893     }
2894   }
2895   
2896   for (ii = Matrix.UpperRow() ; ii >=  Matrix.LowerRow() ; ii--) {
2897     MaxIndex = (ii + UpperBandWidth <= Matrix.UpperRow() ? 
2898                 ii + UpperBandWidth : Matrix.UpperRow()) ;
2899     
2900     for (jj = MaxIndex  ; jj > ii ; jj--) {
2901       
2902       for (kk = 0 ; kk < ArrayDimension ; kk++) {
2903         PolesArray[(ii-1)  * ArrayDimension + kk] -=
2904           PolesArray[(jj - 1) * ArrayDimension + kk] * 
2905             Matrix(ii, jj - ii + LowerBandWidth + 1) ;
2906       }
2907     }
2908     
2909     //fixing a bug PRO18577 to avoid divizion by zero
2910     
2911     Standard_Real divizor = Matrix(ii,LowerBandWidth + 1) ;
2912     Standard_Real Toler = 1.0e-16;
2913     if ( Abs(divizor) > Toler )
2914       Inverse = 1.0e0 / divizor ;
2915     else {
2916       Inverse = 1.0e0;
2917 //      std::cout << "  BSplCLib::SolveBandedSystem() : zero determinant " << std::endl;
2918       ReturnCode = 1;
2919       goto FINISH;
2920     }
2921         
2922     for (kk = 0 ; kk < ArrayDimension ; kk++) {
2923       PolesArray[(ii-1)  * ArrayDimension + kk] *=  Inverse ; 
2924     }
2925   }
2926   FINISH :
2927     return (ReturnCode) ;
2928 }
2929 
2930 //=======================================================================
2931 //function : Solves a LU factored Matrix 
2932 //purpose  : 
2933 //=======================================================================
2934 
2935 Standard_Integer 
2936 BSplCLib::SolveBandedSystem(const math_Matrix&  Matrix,
2937                             const Standard_Integer UpperBandWidth,
2938                             const Standard_Integer LowerBandWidth,
2939                             const Standard_Boolean HomogeneousFlag,
2940                             const Standard_Integer ArrayDimension,
2941                             Standard_Real&   Poles,
2942                             Standard_Real&   Weights) 
2943 {
2944   Standard_Integer ii,
2945   kk,
2946   ErrorCode = 0,
2947   ReturnCode = 0 ;
2948   
2949   Standard_Real   Inverse,
2950   *PolesArray   = &Poles,
2951   *WeightsArray = &Weights ;
2952   
2953   if (Matrix.LowerCol() != 1 || 
2954       Matrix.UpperCol() != UpperBandWidth + LowerBandWidth + 1) {
2955     ReturnCode = 1 ;
2956     goto FINISH ;
2957   }
2958   if (HomogeneousFlag == Standard_False) {
2959     
2960     for (ii = 0 ; ii <  Matrix.UpperRow() - Matrix.LowerRow() + 1; ii++) {
2961       
2962       for (kk = 0 ; kk < ArrayDimension ; kk++) {
2963         PolesArray[ii * ArrayDimension + kk] *=
2964           WeightsArray[ii] ;
2965       }
2966     }
2967   }
2968   ErrorCode = 
2969     BSplCLib::SolveBandedSystem(Matrix,
2970                                 UpperBandWidth,
2971                                 LowerBandWidth,
2972                                 ArrayDimension,
2973                                 Poles) ;
2974   if (ErrorCode != 0) {
2975     ReturnCode = 2 ;
2976     goto FINISH ;
2977   }
2978   ErrorCode = 
2979     BSplCLib::SolveBandedSystem(Matrix,
2980                                 UpperBandWidth,
2981                                 LowerBandWidth,
2982                                 1,
2983                                 Weights) ;
2984   if (ErrorCode != 0) {
2985     ReturnCode = 3 ;
2986     goto FINISH ;
2987   }
2988   if (HomogeneousFlag == Standard_False) {
2989 
2990     for (ii = 0  ; ii < Matrix.UpperRow() - Matrix.LowerRow() + 1 ; ii++) {
2991       Inverse = 1.0e0 / WeightsArray[ii] ;
2992       
2993       for (kk = 0  ; kk < ArrayDimension ; kk++) {
2994         PolesArray[ii * ArrayDimension + kk] *= Inverse ;
2995       }
2996     }
2997   }
2998   FINISH : return (ReturnCode) ;
2999 }
3000 
3001 //=======================================================================
3002 //function : BuildSchoenbergPoints
3003 //purpose  : 
3004 //=======================================================================
3005 
3006 void  BSplCLib::BuildSchoenbergPoints(const Standard_Integer         Degree,
3007                                       const TColStd_Array1OfReal&    FlatKnots,
3008                                       TColStd_Array1OfReal&          Parameters) 
3009 {
3010   Standard_Integer ii,
3011   jj ;
3012   Standard_Real Inverse ;
3013   Inverse = 1.0e0 / (Standard_Real)Degree ;
3014   
3015   for (ii = Parameters.Lower() ;   ii <= Parameters.Upper() ; ii++) {
3016     Parameters(ii) = 0.0e0 ;
3017     
3018     for (jj = 1 ; jj <= Degree ; jj++) {
3019       Parameters(ii) += FlatKnots(jj + ii) ;
3020     } 
3021     Parameters(ii) *= Inverse ; 
3022   }
3023 }
3024 
3025 //=======================================================================
3026 //function : Interpolate
3027 //purpose  : 
3028 //=======================================================================
3029 
3030 void  BSplCLib::Interpolate(const Standard_Integer         Degree,
3031                             const TColStd_Array1OfReal&    FlatKnots,
3032                             const TColStd_Array1OfReal&    Parameters,
3033                             const TColStd_Array1OfInteger& ContactOrderArray,
3034                             const Standard_Integer         ArrayDimension,
3035                             Standard_Real&                 Poles,
3036                             Standard_Integer&              InversionProblem) 
3037 {
3038   Standard_Integer ErrorCode,
3039   UpperBandWidth,
3040   LowerBandWidth ;
3041 //  Standard_Real *PolesArray = &Poles ;
3042   math_Matrix InterpolationMatrix(1, Parameters.Length(),
3043                                   1, 2 * Degree + 1) ;
3044   ErrorCode =
3045   BSplCLib::BuildBSpMatrix(Parameters,
3046                            ContactOrderArray,
3047                            FlatKnots,
3048                            Degree,
3049                            InterpolationMatrix,
3050                            UpperBandWidth,
3051                            LowerBandWidth) ;
3052   if(ErrorCode)
3053     throw Standard_OutOfRange("BSplCLib::Interpolate");
3054 
3055   ErrorCode =
3056   BSplCLib::FactorBandedMatrix(InterpolationMatrix,
3057                            UpperBandWidth,
3058                            LowerBandWidth,
3059                            InversionProblem) ;
3060   if(ErrorCode)
3061     throw Standard_OutOfRange("BSplCLib::Interpolate");
3062 
3063   ErrorCode  =
3064   BSplCLib::SolveBandedSystem(InterpolationMatrix,
3065                               UpperBandWidth,
3066                               LowerBandWidth,
3067                               ArrayDimension,
3068                               Poles) ;
3069   if(ErrorCode)
3070     throw Standard_OutOfRange("BSplCLib::Interpolate");
3071 }
3072 
3073 //=======================================================================
3074 //function : Interpolate
3075 //purpose  : 
3076 //=======================================================================
3077 
3078 void  BSplCLib::Interpolate(const Standard_Integer         Degree,
3079                             const TColStd_Array1OfReal&    FlatKnots,
3080                             const TColStd_Array1OfReal&    Parameters,
3081                             const TColStd_Array1OfInteger& ContactOrderArray,
3082                             const Standard_Integer         ArrayDimension,
3083                             Standard_Real&                 Poles,
3084                             Standard_Real&                 Weights,
3085                             Standard_Integer&              InversionProblem) 
3086 {
3087   Standard_Integer ErrorCode,
3088   UpperBandWidth,
3089   LowerBandWidth ;
3090 
3091   math_Matrix InterpolationMatrix(1, Parameters.Length(),
3092                                   1, 2 * Degree + 1) ;
3093   ErrorCode =
3094   BSplCLib::BuildBSpMatrix(Parameters,
3095                            ContactOrderArray,
3096                            FlatKnots,
3097                            Degree,
3098                            InterpolationMatrix,
3099                            UpperBandWidth,
3100                            LowerBandWidth) ;
3101   if(ErrorCode)
3102     throw Standard_OutOfRange("BSplCLib::Interpolate");
3103 
3104   ErrorCode =
3105   BSplCLib::FactorBandedMatrix(InterpolationMatrix,
3106                            UpperBandWidth,
3107                            LowerBandWidth,
3108                            InversionProblem) ;
3109   if(ErrorCode)
3110     throw Standard_OutOfRange("BSplCLib::Interpolate");
3111 
3112   ErrorCode  =
3113   BSplCLib::SolveBandedSystem(InterpolationMatrix,
3114                               UpperBandWidth,
3115                               LowerBandWidth,
3116                               Standard_False,
3117                               ArrayDimension,
3118                               Poles,
3119                               Weights) ;
3120   if(ErrorCode)
3121     throw Standard_OutOfRange("BSplCLib::Interpolate");
3122 }
3123 
3124 //=======================================================================
3125 //function : Evaluates a Bspline function : uses the ExtrapMode 
3126 //purpose  : the function is extrapolated using the Taylor expansion
3127 //           of degree ExtrapMode[0] to the left and the Taylor
3128 //           expansion of degree ExtrapMode[1] to the right 
3129 //  this evaluates the numerator by multiplying by the weights
3130 //  and evaluating it but does not call RationalDerivatives after 
3131 //=======================================================================
3132 
3133 void  BSplCLib::Eval
3134 (const Standard_Real                   Parameter,
3135  const Standard_Boolean                PeriodicFlag,
3136  const Standard_Integer                DerivativeRequest,
3137  Standard_Integer&                     ExtrapMode,
3138  const Standard_Integer                Degree,
3139  const  TColStd_Array1OfReal&          FlatKnots, 
3140  const Standard_Integer                ArrayDimension,
3141  Standard_Real&                        Poles,
3142  Standard_Real&                        Weights,
3143  Standard_Real&                        PolesResults,
3144  Standard_Real&                        WeightsResults)
3145 {
3146   Standard_Integer ii,
3147   jj,
3148   kk=0,
3149   Index,
3150   Index1,
3151   Index2,
3152   *ExtrapModeArray,
3153   Modulus,
3154   NewRequest,
3155   ExtrapolatingFlag[2],
3156   ErrorCode,
3157   Order = Degree + 1,
3158   FirstNonZeroBsplineIndex,
3159   LocalRequest = DerivativeRequest ;
3160   Standard_Real  *PResultArray,
3161   *WResultArray,
3162   *PolesArray,
3163   *WeightsArray,
3164   LocalParameter,
3165   Period,
3166   Inverse,
3167   Delta ;
3168   PolesArray = &Poles     ;
3169   WeightsArray = &Weights ;
3170   ExtrapModeArray = &ExtrapMode ;
3171   PResultArray = &PolesResults ;
3172   WResultArray = &WeightsResults ;
3173   LocalParameter = Parameter ;
3174   ExtrapolatingFlag[0] = 
3175     ExtrapolatingFlag[1] = 0 ;
3176   //
3177   // check if we are extrapolating to a degree which is smaller than
3178   // the degree of the Bspline
3179   //
3180   if (PeriodicFlag) {
3181     Period = FlatKnots(FlatKnots.Upper() - 1) - FlatKnots(2) ;
3182 
3183     while (LocalParameter > FlatKnots(FlatKnots.Upper() - 1)) {
3184       LocalParameter -= Period ;
3185     }
3186     
3187     while (LocalParameter < FlatKnots(2)) {
3188       LocalParameter +=  Period ;
3189     }
3190   }
3191   if (Parameter < FlatKnots(2) && 
3192       LocalRequest < ExtrapModeArray[0] &&
3193       ExtrapModeArray[0] < Degree) {
3194     LocalRequest = ExtrapModeArray[0] ;
3195     LocalParameter = FlatKnots(2) ;
3196     ExtrapolatingFlag[0] = 1 ;
3197   }
3198   if (Parameter > FlatKnots(FlatKnots.Upper()-1) &&
3199       LocalRequest < ExtrapModeArray[1]  &&
3200       ExtrapModeArray[1] < Degree) {
3201     LocalRequest = ExtrapModeArray[1] ;
3202     LocalParameter = FlatKnots(FlatKnots.Upper()-1) ;
3203     ExtrapolatingFlag[1] = 1 ;
3204   }
3205   Delta = Parameter - LocalParameter ;
3206   if (LocalRequest >= Order) {
3207     LocalRequest = Degree ;
3208   }
3209   if (PeriodicFlag) {
3210     Modulus = FlatKnots.Length() - Degree -1 ;
3211   }
3212   else {
3213     Modulus = FlatKnots.Length() - Degree ;
3214   }
3215 
3216   BSplCLib_LocalMatrix BsplineBasis (LocalRequest, Order);
3217   ErrorCode =
3218     BSplCLib::EvalBsplineBasis(LocalRequest,
3219                                Order,
3220                                FlatKnots,
3221                                LocalParameter,
3222                                FirstNonZeroBsplineIndex,
3223                                BsplineBasis) ;
3224   if (ErrorCode != 0) {
3225     goto FINISH ;
3226   }
3227   if (ExtrapolatingFlag[0] == 0 && ExtrapolatingFlag[1] == 0) {
3228     Index = 0 ;
3229     Index2 = 0 ;
3230 
3231     for (ii = 1 ; ii <= LocalRequest + 1 ; ii++) {
3232       Index1 = FirstNonZeroBsplineIndex ;
3233 
3234       for (kk = 0 ; kk < ArrayDimension ; kk++) {
3235         PResultArray[Index + kk] = 0.0e0 ;
3236       }
3237       WResultArray[Index] = 0.0e0 ;
3238 
3239       for (jj = 1  ; jj <= Order ; jj++) {
3240         
3241         for (kk = 0 ; kk < ArrayDimension ; kk++) {
3242           PResultArray[Index + kk] += 
3243             PolesArray[(Index1-1) * ArrayDimension + kk] 
3244               * WeightsArray[Index1-1] * BsplineBasis(ii,jj) ;
3245         }
3246         WResultArray[Index2]  += WeightsArray[Index1-1] * BsplineBasis(ii,jj) ;
3247         
3248         Index1 = Index1 % Modulus ;
3249         Index1 += 1 ;
3250       }
3251       Index += ArrayDimension ;
3252       Index2 += 1 ;
3253     }
3254   }
3255   else {
3256     // 
3257     //  store Taylor expansion in LocalRealArray
3258     //
3259     NewRequest = DerivativeRequest ;
3260     if (NewRequest > Degree) {
3261       NewRequest = Degree ;
3262     }
3263     NCollection_LocalArray<Standard_Real> LocalRealArray((LocalRequest + 1)*ArrayDimension);
3264     Index = 0 ;
3265     Inverse = 1.0e0 ;
3266 
3267     for (ii = 1 ; ii <= LocalRequest + 1 ; ii++) {
3268       Index1 = FirstNonZeroBsplineIndex ;
3269       
3270       for (kk = 0 ; kk < ArrayDimension ; kk++) {
3271         LocalRealArray[Index + kk] = 0.0e0 ;
3272       }
3273 
3274       for (jj = 1  ; jj <= Order ; jj++) {
3275 
3276         for (kk = 0 ; kk < ArrayDimension ; kk++) {
3277           LocalRealArray[Index + kk] += 
3278             PolesArray[(Index1-1)*ArrayDimension + kk] * 
3279               WeightsArray[Index1-1] * BsplineBasis(ii,jj) ;
3280         }
3281         Index1 = Index1 % Modulus ;
3282         Index1 += 1 ;
3283       }
3284 
3285       for (kk = 0 ; kk < ArrayDimension ; kk++) {
3286         LocalRealArray[Index + kk] *= Inverse ;
3287       }
3288       Index += ArrayDimension ;
3289       Inverse /= (Standard_Real) ii ;
3290     }
3291     PLib::EvalPolynomial(Delta,
3292                          NewRequest,
3293                          Degree,
3294                          ArrayDimension,
3295                          LocalRealArray[0],
3296                          PolesResults) ;
3297     Index = 0 ;
3298     Inverse = 1.0e0 ;
3299 
3300     for (ii = 1 ; ii <= LocalRequest + 1 ; ii++) {
3301       Index1 = FirstNonZeroBsplineIndex ;
3302       LocalRealArray[Index] = 0.0e0 ;
3303 
3304       for (jj = 1  ; jj <= Order ; jj++) {
3305         LocalRealArray[Index] += 
3306           WeightsArray[Index1-1] * BsplineBasis(ii,jj) ;
3307         Index1 = Index1 % Modulus ;
3308         Index1 += 1 ;
3309       }
3310       LocalRealArray[Index + kk] *= Inverse ;
3311       Index += 1 ;
3312       Inverse /= (Standard_Real) ii ;
3313     }
3314     PLib::EvalPolynomial(Delta,
3315                          NewRequest,
3316                          Degree,
3317                          1,
3318                          LocalRealArray[0],
3319                          WeightsResults) ;
3320   }
3321   FINISH : ;
3322 }
3323 
3324 //=======================================================================
3325 //function : Evaluates a Bspline function : uses the ExtrapMode 
3326 //purpose  : the function is extrapolated using the Taylor expansion
3327 //           of degree ExtrapMode[0] to the left and the Taylor
3328 //           expansion of degree ExtrapMode[1] to the right 
3329 // WARNING : the array Results is supposed to have at least 
3330 // (DerivativeRequest + 1) * ArrayDimension slots and the 
3331 // 
3332 //=======================================================================
3333 
3334 void  BSplCLib::Eval
3335 (const Standard_Real                   Parameter,
3336  const Standard_Boolean                PeriodicFlag,
3337  const Standard_Integer                DerivativeRequest,
3338  Standard_Integer&                     ExtrapMode,
3339  const Standard_Integer                Degree,
3340  const  TColStd_Array1OfReal&          FlatKnots, 
3341  const Standard_Integer                ArrayDimension,
3342  Standard_Real&                        Poles,
3343  Standard_Real&                        Results) 
3344 {
3345   Standard_Integer ii,
3346   jj,
3347   kk,
3348   Index,
3349   Index1,
3350   *ExtrapModeArray,
3351   Modulus,
3352   NewRequest,
3353   ExtrapolatingFlag[2],
3354   ErrorCode,
3355   Order = Degree + 1,
3356   FirstNonZeroBsplineIndex,
3357   LocalRequest = DerivativeRequest ;
3358 
3359   Standard_Real  *ResultArray,
3360   *PolesArray,
3361   LocalParameter,
3362   Period,
3363   Inverse,
3364   Delta ;
3365          
3366   PolesArray = &Poles ;
3367   ExtrapModeArray = &ExtrapMode ;
3368   ResultArray = &Results ;  
3369   LocalParameter = Parameter ;
3370   ExtrapolatingFlag[0] = 
3371     ExtrapolatingFlag[1] = 0 ;
3372   //
3373   // check if we are extrapolating to a degree which is smaller than
3374   // the degree of the Bspline
3375   //
3376   if (PeriodicFlag) {
3377     Period = FlatKnots(FlatKnots.Upper() - 1) - FlatKnots(2) ;
3378 
3379     while (LocalParameter > FlatKnots(FlatKnots.Upper() - 1)) {
3380       LocalParameter -= Period ;
3381     }
3382 
3383     while (LocalParameter < FlatKnots(2)) {
3384       LocalParameter +=  Period ;
3385     }
3386   }
3387   if (Parameter < FlatKnots(2) && 
3388       LocalRequest < ExtrapModeArray[0] &&
3389       ExtrapModeArray[0] < Degree) {
3390     LocalRequest = ExtrapModeArray[0] ;
3391     LocalParameter = FlatKnots(2) ;
3392     ExtrapolatingFlag[0] = 1 ;
3393   }
3394   if (Parameter > FlatKnots(FlatKnots.Upper()-1) &&
3395       LocalRequest < ExtrapModeArray[1]  &&
3396       ExtrapModeArray[1] < Degree) {
3397     LocalRequest = ExtrapModeArray[1] ;
3398     LocalParameter = FlatKnots(FlatKnots.Upper()-1) ;
3399     ExtrapolatingFlag[1] = 1 ;
3400   }
3401   Delta = Parameter - LocalParameter ;
3402   if (LocalRequest >= Order) {
3403     LocalRequest = Degree ;
3404   }
3405   
3406   if (PeriodicFlag) {
3407     Modulus = FlatKnots.Length() - Degree -1 ;
3408   }
3409   else {
3410     Modulus = FlatKnots.Length() - Degree ;
3411   }
3412   
3413   BSplCLib_LocalMatrix BsplineBasis (LocalRequest, Order);
3414   
3415   ErrorCode =
3416     BSplCLib::EvalBsplineBasis(LocalRequest,
3417                                Order,
3418                                FlatKnots,
3419                                LocalParameter,
3420                                FirstNonZeroBsplineIndex,
3421                                BsplineBasis);
3422   if (ErrorCode != 0) {
3423     goto FINISH ;
3424   }
3425   if (ExtrapolatingFlag[0] == 0 && ExtrapolatingFlag[1] == 0) {
3426     Index = 0 ;
3427     
3428     for (ii = 1 ; ii <= LocalRequest + 1 ; ii++) {
3429       Index1 = FirstNonZeroBsplineIndex ;
3430       
3431       for (kk = 0 ; kk < ArrayDimension ; kk++) {
3432         ResultArray[Index + kk] = 0.0e0 ;
3433       }
3434 
3435       for (jj = 1  ; jj <= Order ; jj++) {
3436         
3437         for (kk = 0 ; kk < ArrayDimension ; kk++) {
3438           ResultArray[Index + kk] += 
3439             PolesArray[(Index1-1) * ArrayDimension + kk] * BsplineBasis(ii,jj) ;
3440         }
3441         Index1 = Index1 % Modulus ;
3442         Index1 += 1 ;
3443       }
3444       Index += ArrayDimension ;
3445     }
3446   }
3447   else {
3448     // 
3449     //  store Taylor expansion in LocalRealArray
3450     //
3451     NewRequest = DerivativeRequest ;
3452     if (NewRequest > Degree) {
3453       NewRequest = Degree ;
3454     }
3455     NCollection_LocalArray<Standard_Real> LocalRealArray((LocalRequest + 1)*ArrayDimension);
3456 
3457     Index = 0 ;
3458     Inverse = 1.0e0 ;
3459 
3460     for (ii = 1 ; ii <= LocalRequest + 1 ; ii++) {
3461       Index1 = FirstNonZeroBsplineIndex ;
3462       
3463       for (kk = 0 ; kk < ArrayDimension ; kk++) {
3464         LocalRealArray[Index + kk] = 0.0e0 ;
3465       }
3466 
3467       for (jj = 1  ; jj <= Order ; jj++) {
3468 
3469         for (kk = 0 ; kk < ArrayDimension ; kk++) {
3470           LocalRealArray[Index + kk] += 
3471             PolesArray[(Index1-1)*ArrayDimension + kk] * BsplineBasis(ii,jj) ;
3472         }
3473         Index1 = Index1 % Modulus ;
3474         Index1 += 1 ;
3475       }
3476 
3477       for (kk = 0 ; kk < ArrayDimension ; kk++) {
3478         LocalRealArray[Index + kk] *= Inverse ;
3479       }
3480       Index += ArrayDimension ;
3481       Inverse /= (Standard_Real) ii ;
3482     }
3483     PLib::EvalPolynomial(Delta,
3484                          NewRequest,
3485                          Degree,
3486                          ArrayDimension,
3487                          LocalRealArray[0],
3488                          Results) ;
3489   }
3490   FINISH : ;
3491 }
3492 
3493 //=======================================================================
3494 //function : TangExtendToConstraint 
3495 //purpose  : Extends a Bspline function using the tangency map
3496 // WARNING :  
3497 //  
3498 // 
3499 //=======================================================================
3500 
3501 void  BSplCLib::TangExtendToConstraint
3502 (const  TColStd_Array1OfReal&          FlatKnots, 
3503  const Standard_Real                   C1Coefficient,
3504  const Standard_Integer                NumPoles,
3505  Standard_Real&                        Poles,
3506  const Standard_Integer                CDimension,
3507  const Standard_Integer                CDegree,
3508  const  TColStd_Array1OfReal&          ConstraintPoint, 
3509  const Standard_Integer                Continuity,
3510  const Standard_Boolean                After,
3511  Standard_Integer&                     NbPolesResult,
3512  Standard_Integer&                     NbKnotsResult,
3513  Standard_Real&                        KnotsResult, 
3514  Standard_Real&                        PolesResult) 
3515 {
3516 #ifdef OCCT_DEBUG
3517   if (CDegree<Continuity+1) {
3518     std::cout<<"The BSpline degree must be greater than the order of continuity"<<std::endl;
3519   }
3520 #endif
3521   Standard_Real * Padr = &Poles ;
3522   Standard_Real * KRadr = &KnotsResult ;
3523   Standard_Real * PRadr = &PolesResult ;
3524 
3525 ////////////////////////////////////////////////////////////////////////
3526 //
3527 //    1. calculation of extension nD
3528 //
3529 ////////////////////////////////////////////////////////////////////////
3530 
3531 //  Hermite matrix
3532   Standard_Integer Csize = Continuity + 2;
3533   math_Matrix  MatCoefs(1,Csize, 1,Csize);
3534   if (After) {
3535     PLib::HermiteCoefficients(0, 1,           // Limits 
3536                               Continuity, 0,  // Orders of constraints
3537                               MatCoefs);
3538   }
3539   else {
3540     PLib::HermiteCoefficients(0, 1,           // Limits 
3541                               0, Continuity,  // Orders of constraints
3542                               MatCoefs);    
3543   }
3544 
3545 
3546 //  position at the node of connection
3547   Standard_Real Tbord ;
3548   if (After) {
3549     Tbord = FlatKnots(FlatKnots.Upper()-CDegree);
3550   }
3551   else {
3552     Tbord = FlatKnots(FlatKnots.Lower()+CDegree);
3553   }
3554   Standard_Boolean periodic_flag = Standard_False ;
3555   Standard_Integer ipos, extrap_mode[2], derivative_request = Max(Continuity,1);
3556   extrap_mode[0] = extrap_mode[1] = CDegree;
3557   TColStd_Array1OfReal  EvalBS(1, CDimension * (derivative_request+1)) ; 
3558   Standard_Real * Eadr = (Standard_Real *) &EvalBS(1) ;
3559   BSplCLib::Eval(Tbord,periodic_flag,derivative_request,extrap_mode[0],
3560                   CDegree,FlatKnots,CDimension,Poles,*Eadr);
3561 
3562 //  norm of the tangent at the node of connection
3563   math_Vector Tgte(1,CDimension);
3564 
3565   for (ipos=1;ipos<=CDimension;ipos++) {
3566     Tgte(ipos) = EvalBS(ipos+CDimension);
3567   }
3568   Standard_Real L1=Tgte.Norm();
3569 
3570 
3571 //  matrix of constraints
3572   math_Matrix Contraintes(1,Csize,1,CDimension);
3573   if (After) {
3574 
3575     for (ipos=1;ipos<=CDimension;ipos++) {
3576       Contraintes(1,ipos) = EvalBS(ipos);
3577       Contraintes(2,ipos) = C1Coefficient * EvalBS(ipos+CDimension);
3578       if(Continuity >= 2) Contraintes(3,ipos) = EvalBS(ipos+2*CDimension) * Pow(C1Coefficient,2);
3579       if(Continuity >= 3) Contraintes(4,ipos) = EvalBS(ipos+3*CDimension) * Pow(C1Coefficient,3);
3580       Contraintes(Continuity+2,ipos) = ConstraintPoint(ipos);
3581     }
3582   }
3583   else {
3584 
3585     for (ipos=1;ipos<=CDimension;ipos++) {
3586       Contraintes(1,ipos) = ConstraintPoint(ipos);
3587       Contraintes(2,ipos) = EvalBS(ipos);
3588       if(Continuity >= 1) Contraintes(3,ipos) = C1Coefficient * EvalBS(ipos+CDimension);
3589       if(Continuity >= 2) Contraintes(4,ipos) = EvalBS(ipos+2*CDimension) * Pow(C1Coefficient,2);
3590       if(Continuity >= 3) Contraintes(5,ipos) = EvalBS(ipos+3*CDimension) * Pow(C1Coefficient,3);
3591     }
3592   }
3593 
3594 //  calculate the coefficients of extension
3595   Standard_Integer ii, jj, kk;
3596   TColStd_Array1OfReal ExtraCoeffs(1,Csize*CDimension);
3597   ExtraCoeffs.Init(0.);
3598 
3599   for (ii=1; ii<=Csize; ii++) {
3600 
3601     for (jj=1; jj<=Csize; jj++) {
3602 
3603       for (kk=1; kk<=CDimension; kk++) {
3604         ExtraCoeffs(kk+(jj-1)*CDimension) += MatCoefs(ii,jj)*Contraintes(ii,kk);
3605       }
3606     }
3607   }
3608 
3609 //  calculate the poles of extension
3610   TColStd_Array1OfReal ExtrapPoles(1,Csize*CDimension);
3611   Standard_Real * EPadr = &ExtrapPoles(1) ;
3612   PLib::CoefficientsPoles(CDimension,
3613                           ExtraCoeffs, PLib::NoWeights(),
3614                           ExtrapPoles, PLib::NoWeights());
3615 
3616 //  calculate the nodes of extension with multiplicities
3617   TColStd_Array1OfReal ExtrapNoeuds(1,2);
3618   ExtrapNoeuds(1) = 0.;
3619   ExtrapNoeuds(2) = 1.;
3620   TColStd_Array1OfInteger ExtrapMults(1,2);
3621   ExtrapMults(1) = Csize;
3622   ExtrapMults(2) = Csize;
3623 
3624 // flat nodes of extension
3625   TColStd_Array1OfReal FK2(1, Csize*2);
3626   BSplCLib::KnotSequence(ExtrapNoeuds,ExtrapMults,FK2);
3627 
3628 //  norm of the tangent at the connection point 
3629   if (After) {
3630     BSplCLib::Eval(0.,periodic_flag,1,extrap_mode[0],
3631                   Csize-1,FK2,CDimension,*EPadr,*Eadr);
3632   }
3633   else {
3634     BSplCLib::Eval(1.,periodic_flag,1,extrap_mode[0],
3635                   Csize-1,FK2,CDimension,*EPadr,*Eadr);
3636   }
3637 
3638   for (ipos=1;ipos<=CDimension;ipos++) {
3639     Tgte(ipos) = EvalBS(ipos+CDimension);
3640   }
3641   Standard_Real L2 = Tgte.Norm();
3642 
3643 //  harmonisation of degrees
3644   TColStd_Array1OfReal NewP2(1, (CDegree+1)*CDimension);
3645   TColStd_Array1OfReal NewK2(1, 2);
3646   TColStd_Array1OfInteger NewM2(1, 2);
3647   if (Csize-1<CDegree) {
3648     BSplCLib::IncreaseDegree(Csize-1,CDegree,Standard_False,CDimension,
3649                              ExtrapPoles,ExtrapNoeuds,ExtrapMults,
3650                              NewP2,NewK2,NewM2);
3651   }
3652   else {
3653     NewP2 = ExtrapPoles;
3654     NewK2 = ExtrapNoeuds;
3655     NewM2 = ExtrapMults;
3656   }
3657 
3658 //  flat nodes of extension after harmonization of degrees
3659   TColStd_Array1OfReal NewFK2(1, (CDegree+1)*2);
3660   BSplCLib::KnotSequence(NewK2,NewM2,NewFK2);
3661 
3662 
3663 ////////////////////////////////////////////////////////////////////////
3664 //
3665 //    2.  concatenation C0
3666 //
3667 ////////////////////////////////////////////////////////////////////////
3668 
3669 //  ratio of reparametrization
3670   Standard_Real Ratio=1, Delta;
3671   if ( (L1 > Precision::Confusion()) && (L2 > Precision::Confusion()) ) {
3672     Ratio = L2 / L1;
3673   }
3674   if ( (Ratio < 1.e-5) || (Ratio > 1.e5) ) Ratio = 1;
3675 
3676   if (After) {
3677 //    do not touch the first BSpline
3678     Delta = Ratio*NewFK2(NewFK2.Lower()) - FlatKnots(FlatKnots.Upper());
3679   }
3680   else {
3681 //    do not touch the second BSpline
3682     Delta = Ratio*NewFK2(NewFK2.Upper()) - FlatKnots(FlatKnots.Lower());
3683   }
3684 
3685 //  result of the concatenation
3686   Standard_Integer NbP1 = NumPoles, NbP2 = CDegree+1;
3687   Standard_Integer NbK1 = FlatKnots.Length(), NbK2 = 2*(CDegree+1);
3688   TColStd_Array1OfReal NewPoles (1, (NbP1+ NbP2-1)*CDimension);
3689   TColStd_Array1OfReal NewFlats (1, NbK1+NbK2-CDegree-2);
3690 
3691 //  poles
3692   Standard_Integer indNP, indP, indEP;
3693   if (After) {
3694 
3695     for (ii=1;  ii<=NbP1+NbP2-1; ii++) {
3696 
3697       for (jj=1;  jj<=CDimension; jj++) {
3698         indNP = (ii-1)*CDimension+jj;
3699         indP = (ii-1)*CDimension+jj-1;
3700         indEP = (ii-NbP1)*CDimension+jj;
3701         if (ii<NbP1) NewPoles(indNP) =  Padr[indP];
3702         else NewPoles(indNP) = NewP2(indEP);
3703       }
3704     }
3705   }
3706   else {
3707 
3708     for (ii=1;  ii<=NbP1+NbP2-1; ii++) {
3709 
3710       for (jj=1;  jj<=CDimension; jj++) {
3711         indNP = (ii-1)*CDimension+jj;
3712         indEP = (ii-1)*CDimension+jj;
3713         indP = (ii-NbP2)*CDimension+jj-1;
3714         if (ii<NbP2) NewPoles(indNP) =  NewP2(indEP);
3715         else NewPoles(indNP) = Padr[indP];
3716       }
3717     }
3718   }
3719 
3720 //  flat nodes 
3721   if (After) {
3722 //    start with the nodes of the initial surface
3723 
3724     for (ii=1; ii<NbK1; ii++) {
3725       NewFlats(ii) = FlatKnots(FlatKnots.Lower()+ii-1);
3726     }
3727 //    continue with the reparameterized nodes of the extension
3728 
3729     for (ii=1; ii<=NbK2-CDegree-1; ii++) {
3730       NewFlats(NbK1+ii-1) = Ratio*NewFK2(NewFK2.Lower()+ii+CDegree) - Delta;
3731     }
3732   }
3733   else {
3734 //    start with the reparameterized nodes of the extension
3735 
3736     for (ii=1; ii<NbK2-CDegree; ii++) {
3737       NewFlats(ii) = Ratio*NewFK2(NewFK2.Lower()+ii-1) - Delta;
3738     }
3739 //    continue with the nodes of the initial surface
3740 
3741     for (ii=2; ii<=NbK1; ii++) {
3742       NewFlats(NbK2+ii-CDegree-2) = FlatKnots(FlatKnots.Lower()+ii-1);
3743     }
3744   }
3745 
3746 
3747 ////////////////////////////////////////////////////////////////////////
3748 //
3749 //    3.  reduction of multiplicite at the node of connection
3750 //
3751 ////////////////////////////////////////////////////////////////////////
3752 
3753 //  number of separate nodes
3754   Standard_Integer KLength = 1;
3755 
3756   for (ii=2; ii<=NbK1+NbK2-CDegree-2;ii++) {
3757     if (NewFlats(ii) != NewFlats(ii-1)) KLength++;
3758   }
3759 
3760 //  flat nodes --> nodes + multiplicities
3761   TColStd_Array1OfReal NewKnots (1, KLength);
3762   TColStd_Array1OfInteger NewMults (1, KLength);
3763   NewMults.Init(1);
3764   jj = 1;
3765   NewKnots(jj) = NewFlats(1);
3766 
3767   for (ii=2; ii<=NbK1+NbK2-CDegree-2;ii++) {
3768     if (NewFlats(ii) == NewFlats(ii-1)) NewMults(jj)++;
3769     else {
3770       jj++;
3771       NewKnots(jj) = NewFlats(ii);
3772     }
3773   }
3774 
3775 //  reduction of multiplicity at the second or the last but one node
3776   Standard_Integer Index = 2, M = CDegree;
3777   if (After) Index = KLength-1;
3778   TColStd_Array1OfReal ResultPoles (1, (NbP1+ NbP2-1)*CDimension);
3779   TColStd_Array1OfReal ResultKnots (1, KLength);
3780   TColStd_Array1OfInteger ResultMults (1, KLength);
3781   Standard_Real Tol = 1.e-6;
3782   Standard_Boolean Ok = Standard_True;
3783 
3784   while ( (M>CDegree-Continuity) && Ok) {
3785     Ok = RemoveKnot(Index, M-1, CDegree, Standard_False, CDimension,
3786                     NewPoles, NewKnots, NewMults,
3787                     ResultPoles, ResultKnots, ResultMults, Tol);
3788     if (Ok) M--;
3789   }
3790 
3791   if (M == CDegree) {
3792 //    number of poles of the concatenation
3793     NbPolesResult = NbP1 + NbP2 - 1;
3794 //    the poles of the concatenation
3795     Standard_Integer PLength = NbPolesResult*CDimension;
3796 
3797     for (jj=1; jj<=PLength; jj++) {
3798       PRadr[jj-1] = NewPoles(jj);
3799     }
3800   
3801 //    flat nodes of the concatenation
3802     Standard_Integer ideb = 0;
3803 
3804     for (jj=0; jj<NewKnots.Length(); jj++) {
3805       for (ii=0; ii<NewMults(jj+1); ii++) {
3806         KRadr[ideb+ii] = NewKnots(jj+1);
3807       }
3808       ideb += NewMults(jj+1);
3809     }
3810     NbKnotsResult = ideb;
3811   }
3812 
3813   else {
3814 //    number of poles of the result
3815     NbPolesResult = NbP1 + NbP2 - 1 - CDegree + M;
3816 //    the poles of the result
3817     Standard_Integer PLength = NbPolesResult*CDimension;
3818 
3819     for (jj=0; jj<PLength; jj++) {
3820       PRadr[jj] = ResultPoles(jj+1);
3821     }
3822   
3823 //    flat nodes of the result
3824     Standard_Integer ideb = 0;
3825 
3826     for (jj=0; jj<ResultKnots.Length(); jj++) {
3827       for (ii=0; ii<ResultMults(jj+1); ii++) {
3828         KRadr[ideb+ii] = ResultKnots(jj+1);
3829       }
3830       ideb += ResultMults(jj+1);
3831     }
3832     NbKnotsResult = ideb;
3833   }
3834 }
3835 
3836 //=======================================================================
3837 //function : Resolution
3838 //purpose  : 
3839 //                           d
3840 //  Let C(t) = SUM      Ci Bi(t)  a Bspline curve of degree d  
3841 //            i = 1,n      
3842 //  with nodes tj for j = 1,n+d+1 
3843 //
3844 //
3845 //         '                    C1 - Ci-1   d-1
3846 //  Then C (t) = SUM     d *  ---------  Bi (t) 
3847 //                i = 2,n      ti+d - ti
3848 //
3849 //                          d-1
3850 //  for the base of BSpline  Bi  (t) of degree d-1.
3851 //
3852 //  Consequently the upper bound of the norm of the derivative from C is :
3853 //
3854 //
3855 //                        |  Ci - Ci-1  |
3856 //          d *   Max     |  ---------  |
3857 //              i = 2,n |  ti+d - ti  |
3858 //     
3859 //                                      N(t) 
3860 //  In the rational case set    C(t) = -----
3861 //                                      D(t) 
3862 //
3863 //  
3864 //  D(t) =  SUM    Di Bi(t) 
3865 //        i=1,n
3866 //
3867 //  N(t) =  SUM   Di * Ci Bi(t) 
3868 //          i =1,n
3869 //
3870 //          N'(t)  -    D'(t) C(t) 
3871 //   C'(t) = -----------------------
3872 //                   D(t)
3873 //
3874 //                                   
3875 //   N'(t) - D'(t) C(t) = 
3876 //      
3877 //                     Di * (Ci - C(t)) - Di-1 * (Ci-1 - C(t))    d-1
3878 //      SUM   d *   ---------------------------------------- * Bi  (t)  =
3879 //        i=2,n                   ti+d   - ti
3880 //
3881 //    
3882 //                   Di * (Ci - Cj) - Di-1 * (Ci-1 - Cj)                d-1
3883 // SUM   SUM     d * -----------------------------------  * Betaj(t) * Bi  (t) 
3884 //i=2,n j=1,n               ti+d  - ti  
3885 //  
3886 //
3887 //
3888 //                 Dj Bj(t) 
3889 //    Betaj(t) =   --------
3890 //                 D(t) 
3891 //
3892 //  Betaj(t) form a partition >= 0 of the entity with support
3893 //  tj, tj+d+1. Consequently if Rj = {j-d, ....,  j+d+d+1} 
3894 //  obtain an upper bound of the derivative of C by taking :
3895 //
3896 //
3897 //
3898 //
3899 //
3900 //    
3901 //                         Di * (Ci - Cj) - Di-1 * (Ci-1 - Cj) 
3902 //   Max   Max       d  *  -----------------------------------  
3903 // j=1,n  i dans Rj                   ti+d  - ti  
3904 //
3905 //  --------------------------------------------------------
3906 //
3907 //               Min    Di
3908 //              i =1,n
3909 //  
3910 //
3911 //=======================================================================
3912 
3913 void BSplCLib::Resolution(      Standard_Real&        Poles,
3914                           const Standard_Integer      ArrayDimension,
3915                           const Standard_Integer      NumPoles,
3916                           const TColStd_Array1OfReal* Weights,
3917                           const TColStd_Array1OfReal& FlatKnots,
3918                           const Standard_Integer      Degree,
3919                           const Standard_Real         Tolerance3D,
3920                           Standard_Real&              UTolerance) 
3921 {
3922   Standard_Integer ii,num_poles,ii_index,jj_index,ii_inDim;
3923   Standard_Integer lower,upper,ii_minus,jj,ii_miDim;
3924   Standard_Integer Deg1 = Degree + 1;
3925   Standard_Integer Deg2 = (Degree << 1) + 1;
3926   Standard_Real value,factor,W,min_weights,inverse;
3927   Standard_Real pa_ii_inDim_0, pa_ii_inDim_1, pa_ii_inDim_2, pa_ii_inDim_3;
3928   Standard_Real pa_ii_miDim_0, pa_ii_miDim_1, pa_ii_miDim_2, pa_ii_miDim_3;
3929   Standard_Real wg_ii_index, wg_ii_minus;
3930   Standard_Real *PA,max_derivative;
3931   const Standard_Real * FK = &FlatKnots(FlatKnots.Lower());
3932   PA = &Poles;
3933   max_derivative = 0.0e0;
3934   num_poles = FlatKnots.Length() - Deg1;
3935   switch (ArrayDimension) {
3936   case 2 : {
3937     if (Weights != NULL) {
3938       const Standard_Real * WG = &(*Weights)(Weights->Lower());
3939       min_weights = WG[0];
3940       
3941       for (ii = 1 ; ii < NumPoles ; ii++) {
3942         W = WG[ii];
3943         if (W < min_weights) min_weights = W;
3944       }
3945       
3946       for (ii = 1 ; ii < num_poles ; ii++) {
3947         ii_index = ii % NumPoles;
3948         ii_inDim = ii_index << 1;
3949         ii_minus = (ii - 1) % NumPoles;
3950         ii_miDim = ii_minus << 1;
3951         pa_ii_inDim_0 = PA[ii_inDim]; ii_inDim++;
3952         pa_ii_inDim_1 = PA[ii_inDim];
3953         pa_ii_miDim_0 = PA[ii_miDim]; ii_miDim++;
3954         pa_ii_miDim_1 = PA[ii_miDim];
3955         wg_ii_index   = WG[ii_index];
3956         wg_ii_minus   = WG[ii_minus];
3957         inverse = FK[ii + Degree] - FK[ii];
3958         inverse = 1.0e0 / inverse;
3959         lower = ii - Deg1;
3960         if (lower < 0) lower = 0;
3961         upper = Deg2 + ii;
3962         if (upper > num_poles) upper = num_poles;
3963         
3964         for (jj = lower ; jj < upper ; jj++) {
3965           jj_index = jj % NumPoles;
3966           jj_index = jj_index << 1;
3967           value = 0.0e0;
3968           factor  = (((PA[jj_index] - pa_ii_inDim_0) * wg_ii_index) -
3969                      ((PA[jj_index] - pa_ii_miDim_0) * wg_ii_minus));
3970           if (factor < 0) factor = - factor;
3971           value += factor; jj_index++;
3972           factor  = (((PA[jj_index] - pa_ii_inDim_1) * wg_ii_index) -
3973                      ((PA[jj_index] - pa_ii_miDim_1) * wg_ii_minus));
3974           if (factor < 0) factor = - factor;
3975           value += factor;
3976           value *= inverse;
3977           if (max_derivative < value) max_derivative = value;
3978         }
3979       }
3980       max_derivative /= min_weights;
3981     }
3982     else {
3983       
3984       for (ii = 1 ; ii < num_poles ; ii++) {
3985         ii_index = ii % NumPoles;
3986         ii_index = ii_index << 1;
3987         ii_minus = (ii - 1) % NumPoles;
3988         ii_minus = ii_minus << 1;
3989         inverse = FK[ii + Degree] - FK[ii];
3990         inverse = 1.0e0 / inverse;
3991         value = 0.0e0;
3992         factor = PA[ii_index] - PA[ii_minus];
3993         if (factor < 0) factor = - factor;
3994         value += factor; ii_index++; ii_minus++;
3995         factor = PA[ii_index] - PA[ii_minus];
3996         if (factor < 0) factor = - factor;
3997         value += factor;
3998         value *= inverse;
3999         if (max_derivative < value) max_derivative = value;
4000       }
4001     }
4002     break;
4003   }
4004   case 3 : {
4005     if (Weights != NULL) {
4006       const Standard_Real * WG = &(*Weights)(Weights->Lower());
4007       min_weights = WG[0];
4008       
4009       for (ii = 1 ; ii < NumPoles ; ii++) {
4010         W = WG[ii];
4011         if (W < min_weights) min_weights = W;
4012       }
4013       
4014       for (ii = 1 ; ii < num_poles ; ii++) {
4015         ii_index = ii % NumPoles;
4016         ii_inDim = (ii_index << 1) + ii_index;
4017         ii_minus = (ii - 1) % NumPoles;
4018         ii_miDim = (ii_minus << 1) + ii_minus;
4019         pa_ii_inDim_0 = PA[ii_inDim]; ii_inDim++;
4020         pa_ii_inDim_1 = PA[ii_inDim]; ii_inDim++;
4021         pa_ii_inDim_2 = PA[ii_inDim];
4022         pa_ii_miDim_0 = PA[ii_miDim]; ii_miDim++;
4023         pa_ii_miDim_1 = PA[ii_miDim]; ii_miDim++;
4024         pa_ii_miDim_2 = PA[ii_miDim];
4025         wg_ii_index   = WG[ii_index];
4026         wg_ii_minus   = WG[ii_minus];
4027         inverse = FK[ii + Degree] - FK[ii];
4028         inverse = 1.0e0 / inverse;
4029         lower = ii - Deg1;
4030         if (lower < 0) lower = 0;
4031         upper = Deg2 + ii;
4032         if (upper > num_poles) upper = num_poles;
4033         
4034         for (jj = lower ; jj < upper ; jj++) {
4035           jj_index = jj % NumPoles;
4036           jj_index = (jj_index << 1) + jj_index;
4037           value = 0.0e0;
4038           factor  = (((PA[jj_index] - pa_ii_inDim_0) * wg_ii_index) -
4039                      ((PA[jj_index] - pa_ii_miDim_0) * wg_ii_minus));
4040           if (factor < 0) factor = - factor;
4041           value += factor; jj_index++;
4042           factor  = (((PA[jj_index] - pa_ii_inDim_1) * wg_ii_index) -
4043                      ((PA[jj_index] - pa_ii_miDim_1) * wg_ii_minus));
4044           if (factor < 0) factor = - factor;
4045           value += factor; jj_index++;
4046           factor  = (((PA[jj_index] - pa_ii_inDim_2) * wg_ii_index) -
4047                      ((PA[jj_index] - pa_ii_miDim_2) * wg_ii_minus));
4048           if (factor < 0) factor = - factor;
4049           value += factor;
4050           value *= inverse;
4051           if (max_derivative < value) max_derivative = value;
4052         }
4053       }
4054       max_derivative /= min_weights;
4055     }
4056     else {
4057       
4058       for (ii = 1 ; ii < num_poles ; ii++) {
4059         ii_index = ii % NumPoles;
4060         ii_index = (ii_index << 1) + ii_index;
4061         ii_minus = (ii - 1) % NumPoles;
4062         ii_minus = (ii_minus << 1) + ii_minus;
4063         inverse = FK[ii + Degree] - FK[ii];
4064         inverse = 1.0e0 / inverse;
4065         value = 0.0e0;
4066         factor = PA[ii_index] - PA[ii_minus];
4067         if (factor < 0) factor = - factor;
4068         value += factor; ii_index++; ii_minus++;
4069         factor = PA[ii_index] - PA[ii_minus];
4070         if (factor < 0) factor = - factor;
4071         value += factor; ii_index++; ii_minus++;
4072         factor = PA[ii_index] - PA[ii_minus];
4073         if (factor < 0) factor = - factor;
4074         value += factor;
4075         value *= inverse;
4076         if (max_derivative < value) max_derivative = value;
4077       }
4078     }
4079     break;
4080   }
4081   case 4 : {
4082     if (Weights != NULL) {
4083       const Standard_Real * WG = &(*Weights)(Weights->Lower());
4084       min_weights = WG[0];
4085       
4086       for (ii = 1 ; ii < NumPoles ; ii++) {
4087         W = WG[ii];
4088         if (W < min_weights) min_weights = W;
4089       }
4090       
4091       for (ii = 1 ; ii < num_poles ; ii++) {
4092         ii_index = ii % NumPoles;
4093         ii_inDim = ii_index << 2;
4094         ii_minus = (ii - 1) % NumPoles;
4095         ii_miDim = ii_minus << 2;
4096         pa_ii_inDim_0 = PA[ii_inDim]; ii_inDim++;
4097         pa_ii_inDim_1 = PA[ii_inDim]; ii_inDim++;
4098         pa_ii_inDim_2 = PA[ii_inDim]; ii_inDim++;
4099         pa_ii_inDim_3 = PA[ii_inDim];
4100         pa_ii_miDim_0 = PA[ii_miDim]; ii_miDim++;
4101         pa_ii_miDim_1 = PA[ii_miDim]; ii_miDim++;
4102         pa_ii_miDim_2 = PA[ii_miDim]; ii_miDim++;
4103         pa_ii_miDim_3 = PA[ii_miDim];
4104         wg_ii_index   = WG[ii_index];
4105         wg_ii_minus   = WG[ii_minus];
4106         inverse = FK[ii + Degree] - FK[ii];
4107         inverse = 1.0e0 / inverse;
4108         lower = ii - Deg1;
4109         if (lower < 0) lower = 0;
4110         upper = Deg2 + ii;
4111         if (upper > num_poles) upper = num_poles;
4112         
4113         for (jj = lower ; jj < upper ; jj++) {
4114           jj_index = jj % NumPoles;
4115           jj_index = jj_index << 2;
4116           value = 0.0e0;
4117           factor  = (((PA[jj_index] - pa_ii_inDim_0) * wg_ii_index) -
4118                      ((PA[jj_index] - pa_ii_miDim_0) * wg_ii_minus));
4119           if (factor < 0) factor = - factor;
4120           value += factor; jj_index++;
4121           factor  = (((PA[jj_index] - pa_ii_inDim_1) * wg_ii_index) -
4122                      ((PA[jj_index] - pa_ii_miDim_1) * wg_ii_minus));
4123           if (factor < 0) factor = - factor;
4124           value += factor; jj_index++;
4125           factor  = (((PA[jj_index] - pa_ii_inDim_2) * wg_ii_index) -
4126                      ((PA[jj_index] - pa_ii_miDim_2) * wg_ii_minus));
4127           if (factor < 0) factor = - factor;
4128           value += factor; jj_index++;
4129           factor  = (((PA[jj_index] - pa_ii_inDim_3) * wg_ii_index) -
4130                      ((PA[jj_index] - pa_ii_miDim_3) * wg_ii_minus));
4131           if (factor < 0) factor = - factor;
4132           value += factor;
4133           value *= inverse;
4134           if (max_derivative < value) max_derivative = value;
4135         }
4136       }
4137       max_derivative /= min_weights;
4138     }
4139     else {
4140       
4141       for (ii = 1 ; ii < num_poles ; ii++) {
4142         ii_index = ii % NumPoles;
4143         ii_index = ii_index << 2;
4144         ii_minus = (ii - 1) % NumPoles;
4145         ii_minus = ii_minus << 2;
4146         inverse = FK[ii + Degree] - FK[ii];
4147         inverse = 1.0e0 / inverse;
4148         value = 0.0e0;
4149         factor = PA[ii_index] - PA[ii_minus];
4150         if (factor < 0) factor = - factor;
4151         value += factor; ii_index++; ii_minus++;
4152         factor = PA[ii_index] - PA[ii_minus];
4153         if (factor < 0) factor = - factor;
4154         value += factor; ii_index++; ii_minus++;
4155         factor = PA[ii_index] - PA[ii_minus];
4156         if (factor < 0) factor = - factor;
4157         value += factor; ii_index++; ii_minus++;
4158         factor = PA[ii_index] - PA[ii_minus];
4159         if (factor < 0) factor = - factor;
4160         value += factor;
4161         value *= inverse;
4162         if (max_derivative < value) max_derivative = value;
4163       }
4164     }
4165     break;
4166   }
4167     default : {
4168       Standard_Integer kk;
4169       if (Weights != NULL) {
4170         const Standard_Real * WG = &(*Weights)(Weights->Lower());
4171         min_weights = WG[0];
4172         
4173         for (ii = 1 ; ii < NumPoles ; ii++) {
4174           W = WG[ii];
4175           if (W < min_weights) min_weights = W;
4176         }
4177         
4178         for (ii = 1 ; ii < num_poles ; ii++) {
4179           ii_index  = ii % NumPoles;
4180           ii_inDim  = ii_index * ArrayDimension;
4181           ii_minus  = (ii - 1) % NumPoles;
4182           ii_miDim  = ii_minus * ArrayDimension;
4183           wg_ii_index   = WG[ii_index];
4184           wg_ii_minus   = WG[ii_minus];
4185           inverse = FK[ii + Degree] - FK[ii];
4186           inverse = 1.0e0 / inverse;
4187           lower = ii - Deg1;
4188           if (lower < 0) lower = 0;
4189           upper = Deg2 + ii;
4190           if (upper > num_poles) upper = num_poles;
4191           
4192           for (jj = lower ; jj < upper ; jj++) {
4193             jj_index = jj % NumPoles;
4194             jj_index *= ArrayDimension;
4195             value = 0.0e0;
4196             
4197             for (kk = 0 ; kk < ArrayDimension ; kk++) {
4198               factor  = (((PA[jj_index + kk] - PA[ii_inDim + kk]) * wg_ii_index) -
4199                          ((PA[jj_index + kk] - PA[ii_miDim + kk]) * wg_ii_minus));
4200               if (factor < 0) factor = - factor;
4201               value += factor;
4202             }
4203             value *= inverse;
4204             if (max_derivative < value) max_derivative = value;
4205           }
4206         }
4207         max_derivative /= min_weights;
4208       }
4209       else {
4210         
4211         for (ii = 1 ; ii < num_poles ; ii++) {
4212           ii_index  = ii % NumPoles;
4213           ii_index *= ArrayDimension;
4214           ii_minus  = (ii - 1) % NumPoles;
4215           ii_minus *= ArrayDimension;
4216           inverse = FK[ii + Degree] - FK[ii];
4217           inverse = 1.0e0 / inverse;
4218           value = 0.0e0;
4219           
4220           for (kk = 0 ; kk < ArrayDimension ; kk++) {
4221             factor = PA[ii_index + kk] - PA[ii_minus + kk];
4222             if (factor < 0) factor = - factor;
4223             value += factor;
4224           }
4225           value *= inverse;
4226           if (max_derivative < value) max_derivative = value;
4227         }
4228       }
4229     }
4230   }
4231   max_derivative *= Degree;
4232   if (max_derivative > RealSmall())
4233     UTolerance = Tolerance3D / max_derivative; 
4234   else
4235     UTolerance = Tolerance3D / RealSmall();
4236 }
4237 
4238 //=======================================================================
4239 // function : Intervals 
4240 // purpose  : 
4241 //=======================================================================
4242 Standard_Integer BSplCLib::Intervals (const TColStd_Array1OfReal& theKnots,
4243                                       const TColStd_Array1OfInteger& theMults,
4244                                       Standard_Integer theDegree,
4245                                       Standard_Boolean isPeriodic,
4246                                       Standard_Integer theContinuity,
4247                                       Standard_Real theFirst,
4248                                       Standard_Real theLast,
4249                                       Standard_Real theTolerance,
4250                                       TColStd_Array1OfReal* theIntervals) 
4251 {
4252   // remove all knots with multiplicity less or equal than (degree - continuity) except first and last
4253   Standard_Integer aFirstIndex = isPeriodic ? 1 : FirstUKnotIndex (theDegree, theMults);
4254   Standard_Integer aLastIndex = isPeriodic ? theKnots.Size() : LastUKnotIndex (theDegree, theMults);
4255   TColStd_Array1OfReal aNewKnots (1, aLastIndex - aFirstIndex + 1);
4256   Standard_Integer aNbNewKnots = 0;
4257   for (Standard_Integer anIndex = aFirstIndex; anIndex <= aLastIndex; anIndex++)
4258   {
4259     if (theMults(anIndex) > (theDegree - theContinuity) ||
4260     anIndex == aFirstIndex ||
4261     anIndex == aLastIndex)
4262     {
4263       aNbNewKnots++;
4264       aNewKnots(aNbNewKnots) = theKnots[anIndex];
4265     }
4266   }
4267   aNewKnots.Resize (1, aNbNewKnots, Standard_True);
4268 
4269   // the range boundaries 
4270   Standard_Real aCurFirst = theFirst;
4271   Standard_Real aCurLast = theLast;
4272   Standard_Real aPeriod = 0.0;
4273   Standard_Integer aFirstPeriod = 0;
4274   Standard_Integer aLastPeriod = 0;
4275   // move boundaries into period
4276   if (isPeriodic)
4277   {
4278     Standard_Real aLower = theKnots.First();
4279     Standard_Real anUpper = theKnots.Last();
4280     aPeriod = anUpper - aLower;
4281     
4282     while (aCurFirst < aLower)
4283     {
4284       aCurFirst += aPeriod;
4285       aFirstPeriod--;
4286     }
4287     while (aCurLast < aLower)
4288     {
4289       aCurLast += aPeriod;
4290       aLastPeriod--;
4291     }
4292     while (aCurFirst >= anUpper)
4293     {
4294       aCurFirst -= aPeriod;
4295       aFirstPeriod += 1;
4296     }
4297     while (aCurLast >= anUpper)
4298     {
4299       aCurLast -= aPeriod;
4300       aLastPeriod += 1;
4301     }
4302   }
4303   // locate the left and nearest knot for boundaries
4304   Standard_Integer anIndex1 = 0;
4305   Standard_Integer anIndex2 = 0;
4306   Standard_Real aDummyDouble;
4307   // we use version of LocateParameter that doesn't need multiplicities
4308   LocateParameter(theDegree, aNewKnots, TColStd_Array1OfInteger(), aCurFirst, Standard_False, 1, aNbNewKnots, anIndex1, aDummyDouble);
4309   LocateParameter(theDegree, aNewKnots, TColStd_Array1OfInteger(), aCurLast, Standard_False, 1, aNbNewKnots, anIndex2, aDummyDouble);
4310   // the case when the beginning of the range coincides with the next knot
4311   if (anIndex1 < aNbNewKnots && Abs(aNewKnots[anIndex1 + 1] - aCurFirst) < theTolerance)
4312   {
4313     anIndex1 += 1;
4314   }
4315   // the case when the ending of the range coincides with the current knot
4316   if (aNbNewKnots && Abs(aNewKnots[anIndex2] - aCurLast) < theTolerance)
4317   {
4318     anIndex2 -= 1;
4319   }
4320   Standard_Integer aNbIntervals = anIndex2 - anIndex1 + 1 + (aLastPeriod - aFirstPeriod) * (aNbNewKnots - 1);
4321   
4322   // fill the interval array
4323   if (theIntervals)
4324   {
4325     theIntervals->Resize (1, aNbIntervals + 1, Standard_False);
4326     if (isPeriodic && aLastPeriod != aFirstPeriod)
4327     {
4328       Standard_Integer anIndex = 1;
4329       // part from the begging of range to the end of the first period
4330       for (Standard_Integer i = anIndex1; i < aNewKnots.Size(); i++, anIndex++)
4331       {
4332         theIntervals->ChangeValue(anIndex) = aNewKnots[i] + aFirstPeriod * aPeriod;
4333       }
4334       // full periods
4335       for (Standard_Integer aPeriodNum = aFirstPeriod + 1; aPeriodNum < aLastPeriod; aPeriodNum++)
4336       {
4337         for (Standard_Integer i = 1; i < aNewKnots.Size(); i++, anIndex++)
4338         {
4339           theIntervals->ChangeValue(anIndex) = aNewKnots[i] + aPeriodNum * aPeriod;
4340         }
4341       }
4342       // part from the begging of the last period to the end of range
4343       for (Standard_Integer i = 1; i <= anIndex2; i++, anIndex++)
4344       {
4345         theIntervals->ChangeValue(anIndex) = aNewKnots[i] + aLastPeriod * aPeriod;
4346       }
4347     }
4348     else
4349     {
4350       Standard_Integer anIndex = 1;
4351       for (Standard_Integer i = anIndex1; i <= anIndex2; i++, anIndex++)
4352       {
4353         theIntervals->ChangeValue(anIndex) = aNewKnots[i] + aFirstPeriod * aPeriod;
4354       }
4355     }
4356     // update the first position (the begging of range doesn't coincide with the knot at anIndex1 in general)
4357     theIntervals->ChangeValue(1) = theFirst;
4358     // write the ending of the range (we didn't write it at all)
4359     theIntervals->ChangeValue(aNbIntervals + 1) = theLast;
4360   }
4361   
4362   return aNbIntervals;
4363 }
4364 
4365 //=======================================================================
4366 // function: FlatBezierKnots
4367 // purpose :
4368 //=======================================================================
4369 
4370 // array of flat knots for bezier curve of maximum 25 degree
4371 static const Standard_Real knots[52] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
4372                                          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
4373 const Standard_Real& BSplCLib::FlatBezierKnots (const Standard_Integer Degree)
4374 {
4375   Standard_OutOfRange_Raise_if (Degree < 1 || Degree > MaxDegree() || MaxDegree() != 25,
4376     "Bezier curve degree greater than maximal supported");
4377 
4378   return knots[25-Degree];
4379 }
Open CASCADE Technology repositoryRSSAtom
