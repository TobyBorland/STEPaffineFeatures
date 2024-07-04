gitprojects / occt.git / blob

commit
 ? search: 
 re
summary | shortlog | log | commit | commitdiff | tree
blame | history | raw | HEAD
0033261: Data Exchange, Step Import - Empty shape after reading process
[occt.git] / src / BSplSLib / BSplSLib.cxx
   1 // Created on: 1991-08-26
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
  17 // Modified RLE Aug 93 - Complete rewrite
  18 // xab  21-Mar-95  implemented cache mechanism
  19 // pmn  25-09-96   Interpolation
  20 // jct  25-09-96 : Correction de l'alloc de LocalArray dans RationalDerivative.
  21 // pmn  07-10-96 : Correction de DN dans le cas rationnal.
  22 // pmn  06-02-97 : Correction des poids dans RationalDerivative. (PRO700)
  23 
  24 #include <BSplCLib.hxx>
  25 #include <BSplSLib.hxx>
  26 #include <gp_Pnt.hxx>
  27 #include <gp_Vec.hxx>
  28 #include <math_Matrix.hxx>
  29 #include <NCollection_LocalArray.hxx>
  30 #include <PLib.hxx>
  31 #include <Standard_ConstructionError.hxx>
  32 #include <Standard_NotImplemented.hxx>
  33 #include <TColgp_Array1OfXYZ.hxx>
  34 #include <TColgp_Array2OfXYZ.hxx>
  35 #include <TColStd_HArray1OfInteger.hxx>
  36 
  37 // for null derivatives
  38 static Standard_Real BSplSLib_zero[3] = {0.0, 0.0, 0.0};
  39 
  40 //=======================================================================
  41 //struct : BSplCLib_DataContainer 
  42 //purpose: Auxiliary structure providing buffers for poles and knots used in
  43 //         evaluation of bspline (allocated in the stack)
  44 //=======================================================================
  45 
  46 struct BSplSLib_DataContainer
  47 {
  48   BSplSLib_DataContainer (Standard_Integer UDegree, Standard_Integer VDegree)
  49   {
  50     (void)UDegree; (void)VDegree; // just to avoid compiler warning in Release mode
  51     Standard_OutOfRange_Raise_if (UDegree > BSplCLib::MaxDegree() ||
  52         VDegree > BSplCLib::MaxDegree() || BSplCLib::MaxDegree() > 25,
  53         "BSplSLib: bspline degree is greater than maximum supported");
  54   }
  55   Standard_Real poles[4*(25+1)*(25+1)];
  56   Standard_Real knots1[2*25];
  57   Standard_Real knots2[2*25];
  58   Standard_Real ders[48];
  59 };
  60 
  61 //**************************************************************************
  62 //                     Evaluation methods
  63 //**************************************************************************
  64 
  65 //=======================================================================
  66 //function : RationalDerivative
  67 //purpose  : computes the rational derivatives when whe have the
  68 //           the derivatives of the homogeneous numerator and the
  69 //           the derivatives of the denominator 
  70 //=======================================================================
  71 
  72 void  BSplSLib::RationalDerivative(const Standard_Integer UDeg,
  73                                    const Standard_Integer VDeg,
  74                                    const Standard_Integer N, 
  75                                    const Standard_Integer M, 
  76                                    Standard_Real& HDerivatives,
  77                                    Standard_Real& RDerivatives,
  78                                    const Standard_Boolean All)
  79 {
  80   //
  81   //  if All is True all derivatives are computed. if Not only
  82   //  the requested N, M is computed 
  83   //
  84   //                                           Numerator(u,v) 
  85   //   let f(u,v) be a rational function  = ------------------
  86   //                                         Denominator(u,v)
  87   //
  88   //
  89   //   Let (N,M) the order of the derivatives we want : then since
  90   //   we have :
  91   //
  92   //         Numerator = f * Denominator 
  93   //
  94   //   we derive :
  95   //
  96   //   (N,M)         1           (         (N M)                  (p q)            (N -p M-q) )
  97   //  f       =  ------------    (  Numerator   -  SUM SUM  a   * f    * Denominator          )
  98   //                       (0,0) (                 p<N q<M   p q                              )
  99   //             Denominator
 100   //
 101   //   with :
 102   //
 103   //              ( N )  ( M )
 104   //    a      =  (   )  (   )
 105   //     p q      ( p )  ( q )
 106   //
 107   //  
 108   //   HDerivatives is an array where derivatives are stored in the following form
 109   //   Numerator is assumee to have 3 functions that is a vector of dimension
 110   //   3
 111   // 
 112   //             (0,0)           (0,0)                       (0, DegV)           (0, DegV) 
 113   //     Numerator     Denominator      ...         Numerator          Denominator
 114   //
 115   //             (1,0)           (1,0)                       (1, DegV)           (1, DegV) 
 116   //     Numerator     Denominator      ...         Numerator          Denominator
 117   //
 118   //         ...........................................................
 119   //
 120   //
 121   //            (DegU,0)        (DegU,0)                  (DegU, DegV)           (DegU, DegV) 
 122   //     Numerator     Denominator      ...         Numerator          Denominator
 123   //
 124   //
 125   Standard_Integer ii,jj,pp,qq,index,index1,index2;
 126   Standard_Integer M1,M3,M4,N1,iiM1,iiM3,jjM1,ppM1,ppM3;
 127   Standard_Integer MinN,MinN1,MinM,MinM1;
 128   Standard_Integer index_u,index_u1,index_v,index_v1,index_w;
 129 
 130   M1 = M + 1;
 131   N1 = N + 1;
 132   ii = N1 * M1;
 133   M3 = (M1 << 1) + M1;
 134   M4 = (VDeg + 1) << 2;
 135   
 136   NCollection_LocalArray<Standard_Real> StoreDerivatives (All ? 0 : ii * 3);
 137   Standard_Real *RArray = (All ? &RDerivatives : (Standard_Real*)StoreDerivatives);
 138   NCollection_LocalArray<Standard_Real> StoreW (ii);  
 139   Standard_Real *HomogeneousArray = &HDerivatives; 
 140   Standard_Real denominator,Pii,Pip,Pjq;
 141   
 142   denominator = 1.0e0 / HomogeneousArray[3];
 143   index_u  = 0;
 144   index_u1 = 0;
 145   if (UDeg < N) MinN = UDeg;
 146   else          MinN = N;
 147   if (VDeg < M) MinM = VDeg;
 148   else          MinM = M;
 149   MinN1 = MinN + 1;
 150   MinM1 = MinM + 1;
 151   iiM1 = - M1;
 152 
 153   for (ii = 0 ; ii < MinN1 ; ii++) {
 154     iiM1    += M1;
 155     index_v  = index_u;
 156     index_v1 = index_u1;
 157     index_w  = iiM1;
 158     
 159     for (jj = 0 ; jj < MinM1 ; jj++) {
 160       RArray[index_v++] = HomogeneousArray[index_v1++];
 161       RArray[index_v++] = HomogeneousArray[index_v1++];
 162       RArray[index_v++] = HomogeneousArray[index_v1++];
 163       StoreW[index_w++] = HomogeneousArray[index_v1++];
 164     }
 165 
 166     for (jj = MinM1 ; jj < M1 ; jj++) {
 167       RArray[index_v++] = 0.;
 168       RArray[index_v++] = 0.;
 169       RArray[index_v++] = 0.;
 170       StoreW[index_w++] = 0.;
 171     }
 172     index_u1 += M4;
 173     index_u  += M3;
 174   }
 175   index_v = MinN1 * M3;
 176   index_w = MinN1 * M1;
 177   
 178   for (ii = MinN1 ; ii < N1 ; ii++) {
 179     
 180     for (jj = 0 ; jj < M1 ; jj++) {  
 181       RArray[index_v++] = 0.0e0;
 182       RArray[index_v++] = 0.0e0;
 183       RArray[index_v++] = 0.0e0;
 184       StoreW[index_w++] = 0.0e0;
 185     }
 186   } 
 187 
 188   // ---------------  Calculation ----------------
 189 
 190   iiM1 = - M1;
 191   iiM3 = - M3;
 192 
 193   for (ii = 0 ; ii <= N  ; ii++) {
 194     iiM1  += M1;
 195     iiM3  += M3;
 196     index1 = iiM3 - 3;
 197     jjM1   = iiM1;
 198     
 199     for (jj = 0 ; jj <= M ; jj++) {
 200       jjM1 ++;
 201       ppM1    = - M1;
 202       ppM3    = - M3;
 203       index1 += 3;
 204 
 205       for (pp = 0 ; pp < ii ; pp++) {
 206         ppM1  += M1;
 207         ppM3  += M3;
 208         index  = ppM3;
 209         index2 = jjM1 - ppM1;
 210         Pip    = PLib::Bin(ii,pp);
 211         
 212         for (qq = 0 ; qq <= jj ; qq++) {
 213           index2--;
 214           Pjq    = Pip * PLib::Bin(jj,qq) * StoreW[index2]; 
 215           RArray[index1] -= Pjq * RArray[index]; index++; index1++;
 216           RArray[index1] -= Pjq * RArray[index]; index++; index1++;
 217           RArray[index1] -= Pjq * RArray[index]; index++;
 218           index1 -= 2;
 219         }
 220       }
 221       index  = iiM3;
 222       index2 = jj + 1;
 223       Pii = PLib::Bin(ii,ii);
 224 
 225       for (qq = 0 ; qq < jj ; qq++) {
 226         index2--;
 227         Pjq = Pii * PLib::Bin(jj,qq) * StoreW[index2]; 
 228         RArray[index1] -= Pjq * RArray[index]; index++; index1++;
 229         RArray[index1] -= Pjq * RArray[index]; index++; index1++;
 230         RArray[index1] -= Pjq * RArray[index]; index++;
 231         index1 -= 2;
 232       }
 233       RArray[index1] *= denominator; index1++;
 234       RArray[index1] *= denominator; index1++;
 235       RArray[index1] *= denominator;
 236       index1 -= 2;
 237     }
 238   }
 239   if (!All) {
 240     RArray = &RDerivatives;
 241     index = N * M1 + M;
 242     index = (index << 1) + index;
 243     RArray[0] = StoreDerivatives[index]; index++;
 244     RArray[1] = StoreDerivatives[index]; index++;
 245     RArray[2] = StoreDerivatives[index];
 246   }
 247 }
 248 
 249 //=======================================================================
 250 //function : PrepareEval
 251 //purpose  : 
 252 //=======================================================================
 253 
 254 //
 255 // PrepareEval :
 256 //
 257 // Prepare all data for computing points :
 258 //  local arrays of knots
 259 //  local array  of poles (multiplied by the weights if rational)
 260 //
 261 //  The first direction to compute (smaller degree) is returned 
 262 //  and the poles are stored according to this direction.
 263 
 264 static Standard_Boolean  PrepareEval (const Standard_Real            U,
 265                                       const Standard_Real            V,
 266                                       const Standard_Integer         Uindex,
 267                                       const Standard_Integer         Vindex,
 268                                       const Standard_Integer         UDegree,
 269                                       const Standard_Integer         VDegree,
 270                                       const Standard_Boolean         URat,
 271                                       const Standard_Boolean         VRat,
 272                                       const Standard_Boolean         UPer,
 273                                       const Standard_Boolean         VPer,
 274                                       const TColgp_Array2OfPnt&      Poles,
 275                                       const TColStd_Array2OfReal*    Weights,
 276                                       const TColStd_Array1OfReal&    UKnots,
 277                                       const TColStd_Array1OfReal&    VKnots,
 278                                       const TColStd_Array1OfInteger* UMults,
 279                                       const TColStd_Array1OfInteger* VMults,
 280                                       Standard_Real& u1,     // first  parameter to use
 281                                       Standard_Real& u2,     // second parameter to use
 282                                       Standard_Integer& d1,  // first degree
 283                                       Standard_Integer& d2,  // second degree
 284                                       Standard_Boolean& rational,
 285                                       BSplSLib_DataContainer& dc)
 286   {
 287   rational = URat || VRat;
 288   Standard_Integer uindex = Uindex;
 289   Standard_Integer vindex = Vindex;
 290   Standard_Integer UKLower = UKnots.Lower();
 291   Standard_Integer UKUpper = UKnots.Upper();
 292   Standard_Integer VKLower = VKnots.Lower();
 293   Standard_Integer VKUpper = VKnots.Upper();
 294 
 295   if (UDegree <= VDegree)
 296     {
 297     // compute the indices
 298     if (uindex < UKLower || uindex > UKUpper)
 299       BSplCLib::LocateParameter(UDegree,UKnots,UMults,U,UPer,uindex,u1);
 300     else
 301       u1 = U;
 302 
 303     if (vindex < VKLower || vindex > VKUpper)
 304       BSplCLib::LocateParameter(VDegree,VKnots,VMults,V,VPer,vindex,u2);
 305     else
 306       u2 = V;
 307 
 308     // get the knots
 309     d1 = UDegree;
 310     d2 = VDegree;
 311     BSplCLib::BuildKnots(UDegree,uindex,UPer,UKnots,UMults,*dc.knots1);
 312     BSplCLib::BuildKnots(VDegree,vindex,VPer,VKnots,VMults,*dc.knots2);
 313     
 314     if (UMults == NULL)
 315       uindex -= UKLower + UDegree;
 316     else
 317       uindex  = BSplCLib::PoleIndex(UDegree,uindex,UPer,*UMults);
 318 
 319     if (VMults == NULL)
 320       vindex -= VKLower + VDegree;
 321     else
 322       vindex  = BSplCLib::PoleIndex(VDegree,vindex,VPer,*VMults);
 323 
 324     // get the poles
 325     Standard_Integer i,j,ip,jp;
 326     Standard_Real w, *pole = dc.poles;
 327     d1 = UDegree;
 328     d2 = VDegree;
 329     Standard_Integer PLowerRow = Poles.LowerRow();
 330     Standard_Integer PUpperRow = Poles.UpperRow();
 331     Standard_Integer PLowerCol = Poles.LowerCol();
 332     Standard_Integer PUpperCol = Poles.UpperCol();
 333 
 334     // verify if locally non rational
 335     if (rational) 
 336       {
 337       rational = Standard_False;
 338       ip = PLowerRow + uindex;
 339       jp = PLowerCol + vindex;
 340       
 341       if(ip < PLowerRow)        ip = PUpperRow;
 342       if(jp < PLowerCol)        jp = PUpperCol;
 343       
 344       w  = Weights->Value(ip,jp);
 345       Standard_Real eps = Epsilon(w);
 346       Standard_Real dw;
 347 
 348       for (i = 0; i <= UDegree && !rational; i++)
 349         {
 350         jp = PLowerCol + vindex;
 351 
 352         if(jp < PLowerCol)
 353           jp = PUpperCol;
 354 
 355         for (j = 0; j <= VDegree && !rational; j++)
 356           {
 357           dw = Weights->Value(ip,jp) - w;
 358           if (dw < 0)
 359             dw = - dw;
 360 
 361           rational = (dw > eps);
 362 
 363           jp++;
 364 
 365           if (jp > PUpperCol)
 366             jp = PLowerCol;
 367           }
 368 
 369         ip++;
 370 
 371         if (ip > PUpperRow)
 372           ip = PLowerRow;
 373 
 374         }
 375       }
 376 
 377     // copy the poles
 378     ip = PLowerRow + uindex;
 379 
 380     if(ip < PLowerRow)
 381       ip = PUpperRow;
 382 
 383     if (rational)
 384       {
 385       for (i = 0; i <= d1; i++)
 386         {
 387         jp = PLowerCol + vindex;
 388 
 389         if(jp < PLowerCol)
 390           jp = PUpperCol;
 391 
 392         for (j = 0; j <= d2; j++)
 393           {
 394           const gp_Pnt& P = Poles  .Value(ip,jp);
 395           pole[3] = w     = Weights->Value(ip,jp);
 396           pole[0] = P.X() * w;
 397           pole[1] = P.Y() * w;
 398           pole[2] = P.Z() * w;
 399           pole   += 4;
 400           jp++;
 401 
 402           if (jp > PUpperCol)
 403             jp = PLowerCol;
 404           }
 405 
 406         ip++;
 407 
 408         if (ip > PUpperRow)
 409           ip = PLowerRow;
 410 
 411         }
 412       }
 413     else
 414       {
 415       for (i = 0; i <= d1; i++)
 416         {
 417         jp = PLowerCol + vindex;
 418         
 419         if(jp < PLowerCol)
 420           jp = PUpperCol;
 421 
 422         for (j = 0; j <= d2; j++)
 423           {
 424           const gp_Pnt& P = Poles.Value(ip,jp);
 425           pole[0] = P.X();
 426           pole[1] = P.Y();
 427           pole[2] = P.Z();
 428           pole   += 3;
 429           jp++;
 430 
 431           if (jp > PUpperCol)
 432             jp = PLowerCol;
 433           }
 434 
 435         ip++;
 436         
 437         if (ip > PUpperRow)
 438           ip = PLowerRow;
 439         }
 440       }
 441 
 442     return Standard_True;
 443     }
 444   else
 445     {
 446     // compute the indices
 447     if (uindex < UKLower || uindex > UKUpper)
 448       BSplCLib::LocateParameter(UDegree,UKnots,UMults,U,UPer,uindex,u2);
 449     else
 450       u2 = U;
 451 
 452     if (vindex < VKLower || vindex > VKUpper)
 453       BSplCLib::LocateParameter(VDegree,VKnots,VMults,V,VPer,vindex,u1);
 454     else
 455       u1 = V;
 456 
 457     // get the knots
 458 
 459     d2 = UDegree;
 460     d1 = VDegree;
 461 
 462     BSplCLib::BuildKnots(UDegree,uindex,UPer,UKnots,UMults,*dc.knots2);
 463     BSplCLib::BuildKnots(VDegree,vindex,VPer,VKnots,VMults,*dc.knots1);
 464 
 465     if (UMults == NULL)
 466       uindex -= UKLower + UDegree;
 467     else
 468       uindex  = BSplCLib::PoleIndex(UDegree,uindex,UPer,*UMults);
 469 
 470     if (VMults == NULL)
 471       vindex -= VKLower + VDegree;
 472     else
 473       vindex  = BSplCLib::PoleIndex(VDegree,vindex,VPer,*VMults);
 474 
 475     // get the poles
 476     Standard_Integer i,j,ip,jp;
 477     Standard_Real w, *pole = dc.poles;
 478     d1 = VDegree;
 479     d2 = UDegree;
 480     Standard_Integer PLowerRow = Poles.LowerRow();
 481     Standard_Integer PUpperRow = Poles.UpperRow();
 482     Standard_Integer PLowerCol = Poles.LowerCol();
 483     Standard_Integer PUpperCol = Poles.UpperCol();
 484 
 485     // verify if locally non rational
 486     if (rational)
 487       { 
 488       rational = Standard_False;
 489       ip = PLowerRow + uindex;
 490       jp = PLowerCol + vindex;
 491       
 492       if(ip < PLowerRow)
 493         ip = PUpperRow;
 494 
 495       if(jp < PLowerCol)
 496         jp = PUpperCol;
 497 
 498       w  = Weights->Value(ip,jp);
 499       Standard_Real eps = Epsilon(w);
 500       Standard_Real dw;
 501 
 502       for (i = 0; i <= UDegree && !rational; i++)
 503         {
 504         jp = PLowerCol + vindex;
 505 
 506         if(jp < PLowerCol)
 507           jp = PUpperCol;
 508 
 509         for (j = 0; j <= VDegree && !rational; j++)
 510           {
 511           dw = Weights->Value(ip,jp) - w;
 512           if (dw < 0) dw = - dw;
 513           rational = dw > eps;
 514 
 515           jp++;
 516 
 517           if (jp > PUpperCol)
 518             jp = PLowerCol;
 519           }
 520 
 521         ip++;
 522 
 523         if (ip > PUpperRow)
 524           ip = PLowerRow;
 525 
 526         }
 527       }
 528 
 529     // copy the poles
 530     jp = PLowerCol + vindex;
 531 
 532     if(jp < PLowerCol)
 533       jp = PUpperCol;
 534 
 535     if (rational)
 536       {
 537       for (i = 0; i <= d1; i++)
 538         {
 539         ip = PLowerRow + uindex;
 540 
 541         if(ip < PLowerRow)
 542           ip = PUpperRow;
 543 
 544         for (j = 0; j <= d2; j++)
 545           {
 546           const gp_Pnt& P = Poles.Value(ip,jp);
 547           pole[3] = w     = Weights->Value(ip,jp);
 548           pole[0] = P.X() * w;
 549           pole[1] = P.Y() * w;
 550           pole[2] = P.Z() * w;
 551           pole   += 4;
 552           ip++;
 553 
 554           if (ip > PUpperRow)
 555             ip = PLowerRow;
 556 
 557           }
 558 
 559         jp++;
 560 
 561         if (jp > PUpperCol)
 562           jp = PLowerCol;
 563 
 564         }
 565       }
 566     else
 567       {
 568       for (i = 0; i <= d1; i++)
 569         {
 570         ip = PLowerRow + uindex;
 571 
 572         if(ip < PLowerRow)
 573           ip = PUpperRow;
 574 
 575         if(ip > PUpperRow)
 576           ip = PLowerRow;
 577 
 578         for (j = 0; j <= d2; j++)
 579           {
 580           const gp_Pnt& P = Poles.Value(ip,jp);
 581           pole[0] = P.X();
 582           pole[1] = P.Y();
 583           pole[2] = P.Z();
 584           pole   += 3;
 585           ip++;
 586 
 587           if (ip > PUpperRow)
 588             ip = PLowerRow;
 589           }
 590 
 591         jp++;
 592 
 593         if (jp > PUpperCol)
 594           jp = PLowerCol;
 595 
 596         }
 597       }
 598 
 599     return Standard_False;
 600     }
 601   }
 602 
 603 //=======================================================================
 604 //function : D0
 605 //purpose  : 
 606 //=======================================================================
 607 
 608 void  BSplSLib::D0
 609 (const Standard_Real U, 
 610  const Standard_Real V, 
 611  const Standard_Integer UIndex, 
 612  const Standard_Integer VIndex,
 613  const TColgp_Array2OfPnt& Poles,
 614  const TColStd_Array2OfReal* Weights,
 615  const TColStd_Array1OfReal& UKnots,
 616  const TColStd_Array1OfReal& VKnots,
 617  const TColStd_Array1OfInteger* UMults,
 618  const TColStd_Array1OfInteger* VMults,
 619  const Standard_Integer UDegree,
 620  const Standard_Integer VDegree,
 621  const Standard_Boolean URat,
 622  const Standard_Boolean VRat,
 623  const Standard_Boolean UPer,
 624  const Standard_Boolean VPer,
 625  gp_Pnt& P)
 626 {
 627 //  Standard_Integer k ;
 628   Standard_Real W ;
 629   HomogeneousD0(U,
 630                 V, 
 631                 UIndex, 
 632                 VIndex,
 633                 Poles,
 634                 Weights,
 635                 UKnots,
 636                 VKnots,
 637                 UMults,
 638                 VMults,
 639                 UDegree,
 640                 VDegree,
 641                 URat,
 642                 VRat,
 643                 UPer,
 644                 VPer,
 645                 W,
 646                 P) ;
 647   P.SetX(P.X() / W);
 648   P.SetY(P.Y() / W);
 649   P.SetZ(P.Z() / W);
 650 }
 651 
 652 //=======================================================================
 653 //function : D0
 654 //purpose  : 
 655 //=======================================================================
 656 
 657 void  BSplSLib::HomogeneousD0
 658 (const Standard_Real U, 
 659  const Standard_Real V, 
 660  const Standard_Integer UIndex, 
 661  const Standard_Integer VIndex,
 662  const TColgp_Array2OfPnt& Poles,
 663  const TColStd_Array2OfReal* Weights,
 664  const TColStd_Array1OfReal& UKnots,
 665  const TColStd_Array1OfReal& VKnots,
 666  const TColStd_Array1OfInteger* UMults,
 667  const TColStd_Array1OfInteger* VMults,
 668  const Standard_Integer UDegree,
 669  const Standard_Integer VDegree,
 670  const Standard_Boolean URat,
 671  const Standard_Boolean VRat,
 672  const Standard_Boolean UPer,
 673  const Standard_Boolean VPer,
 674  Standard_Real & W,
 675  gp_Pnt& P)
 676 {
 677   Standard_Boolean rational;
 678 //  Standard_Integer k,dim;
 679   Standard_Integer dim;
 680   Standard_Real u1,u2;
 681   Standard_Integer d1,d2;
 682   W = 1.0e0 ;
 683   
 684   BSplSLib_DataContainer dc (UDegree, VDegree);
 685   PrepareEval(U,V,UIndex,VIndex,UDegree,VDegree,URat,VRat,UPer,VPer,
 686               Poles,Weights,UKnots,VKnots,UMults,VMults,
 687               u1,u2,d1,d2,rational,dc);
 688   if (rational) {
 689     dim = 4;
 690     BSplCLib::Eval(u1,d1,*dc.knots1,dim * (d2 + 1),*dc.poles);
 691     BSplCLib::Eval(u2,d2,*dc.knots2,dim,*dc.poles);
 692     W = dc.poles[3];
 693     P.SetX(dc.poles[0]);
 694     P.SetY(dc.poles[1]);
 695     P.SetZ(dc.poles[2]);
 696   }
 697   else {
 698     dim = 3;
 699     BSplCLib::Eval(u1,d1,*dc.knots1,dim * (d2 + 1),*dc.poles);
 700     BSplCLib::Eval(u2,d2,*dc.knots2,dim,*dc.poles);
 701     P.SetX(dc.poles[0]);
 702     P.SetY(dc.poles[1]);
 703     P.SetZ(dc.poles[2]);
 704   }
 705 }
 706 
 707 //=======================================================================
 708 //function : D1
 709 //purpose  : 
 710 //=======================================================================
 711 
 712 void  BSplSLib::D1
 713 (const Standard_Real U,
 714  const Standard_Real V,
 715  const Standard_Integer UIndex,
 716  const Standard_Integer VIndex,
 717  const TColgp_Array2OfPnt& Poles,
 718  const TColStd_Array2OfReal* Weights,
 719  const TColStd_Array1OfReal& UKnots,
 720  const TColStd_Array1OfReal& VKnots,
 721  const TColStd_Array1OfInteger* UMults,
 722  const TColStd_Array1OfInteger* VMults,
 723  const Standard_Integer UDegree,
 724  const Standard_Integer VDegree,
 725  const Standard_Boolean URat,
 726  const Standard_Boolean VRat,
 727  const Standard_Boolean UPer,
 728  const Standard_Boolean VPer,
 729  gp_Pnt& P,
 730  gp_Vec& Vu,
 731  gp_Vec& Vv)
 732 {
 733   Standard_Boolean rational;
 734 //  Standard_Integer k,dim,dim2;
 735   Standard_Integer dim,dim2;
 736   Standard_Real u1,u2;
 737   Standard_Integer d1,d2;
 738   Standard_Real *result, *resVu, *resVv;
 739   BSplSLib_DataContainer dc (UDegree, VDegree);
 740   if (PrepareEval
 741     (U,V,UIndex,VIndex,UDegree,VDegree,URat,VRat,UPer,VPer,
 742      Poles,Weights,UKnots,VKnots,UMults,VMults,
 743      u1,u2,d1,d2,rational,dc)) {
 744     if (rational) {
 745       dim  = 4;
 746       dim2 = (d2 + 1) << 2;
 747       BSplCLib::Bohm(u1,d1,1,*dc.knots1,dim2,*dc.poles);
 748       BSplCLib::Bohm(u2,d2,1,*dc.knots2,dim ,*dc.poles);
 749       BSplCLib::Eval(u2,d2,  *dc.knots2,dim ,*(dc.poles + dim2));
 750       BSplSLib::RationalDerivative(d1,d2,1,1,*dc.poles,*dc.ders);
 751       result = dc.ders;
 752       resVu  = result + 6;
 753       resVv  = result + 3;
 754     }
 755     else {
 756       dim  = 3;
 757       dim2 = d2 + 1;
 758       dim2 = (dim2 << 1) + dim2;
 759       BSplCLib::Bohm(u1,d1,1,*dc.knots1,dim2,*dc.poles);
 760       BSplCLib::Bohm(u2,d2,1,*dc.knots2,dim ,*dc.poles);
 761       BSplCLib::Eval(u2,d2,  *dc.knots2,dim ,*(dc.poles + dim2));
 762       result = dc.poles;
 763       resVu  = result + dim2;
 764       resVv  = result + 3;
 765     }
 766   }
 767   else {
 768     if (rational) {
 769       dim  = 4;
 770       dim2 = (d2 + 1) << 2;
 771       BSplCLib::Bohm(u1,d1,1,*dc.knots1,dim2,*dc.poles);
 772       BSplCLib::Bohm(u2,d2,1,*dc.knots2,dim ,*dc.poles);
 773       BSplCLib::Eval(u2,d2,  *dc.knots2,dim ,*(dc.poles + dim2));
 774       BSplSLib::RationalDerivative(d1,d2,1,1,*dc.poles,*dc.ders);
 775       result = dc.ders;
 776       resVu  = result + 3;
 777       resVv  = result + 6;
 778     }
 779     else {
 780       dim  = 3;
 781       dim2 = d2 + 1;
 782       dim2 = (dim2 << 1) + dim2;
 783       BSplCLib::Bohm(u1,d1,1,*dc.knots1,dim2,*dc.poles);
 784       BSplCLib::Bohm(u2,d2,1,*dc.knots2,dim ,*dc.poles);
 785       BSplCLib::Eval(u2,d2  ,*dc.knots2,dim ,*(dc.poles + dim2));
 786       result = dc.poles;
 787       resVu  = result + 3;
 788       resVv  = result + dim2;
 789     }
 790   }
 791   
 792   P .SetX(result[0]);
 793   Vu.SetX(resVu [0]);
 794   Vv.SetX(resVv [0]);
 795   
 796   P .SetY(result[1]);
 797   Vu.SetY(resVu [1]);
 798   Vv.SetY(resVv [1]);
 799   
 800   P .SetZ(result[2]);
 801   Vu.SetZ(resVu [2]);
 802   Vv.SetZ(resVv [2]);
 803 }
 804 
 805 //=======================================================================
 806 //function : D1
 807 //purpose  : 
 808 //=======================================================================
 809 
 810 void  BSplSLib::HomogeneousD1
 811 (const Standard_Real U,
 812  const Standard_Real V,
 813  const Standard_Integer UIndex,
 814  const Standard_Integer VIndex,
 815  const TColgp_Array2OfPnt& Poles,
 816  const TColStd_Array2OfReal* Weights,
 817  const TColStd_Array1OfReal& UKnots,
 818  const TColStd_Array1OfReal& VKnots,
 819  const TColStd_Array1OfInteger* UMults,
 820  const TColStd_Array1OfInteger* VMults,
 821  const Standard_Integer UDegree,
 822  const Standard_Integer VDegree,
 823  const Standard_Boolean URat,
 824  const Standard_Boolean VRat,
 825  const Standard_Boolean UPer,
 826  const Standard_Boolean VPer,
 827  gp_Pnt& N,
 828  gp_Vec& Nu,
 829  gp_Vec& Nv,
 830  Standard_Real& D,
 831  Standard_Real& Du,
 832  Standard_Real& Dv)
 833 {
 834   Standard_Boolean rational;
 835 //  Standard_Integer k,dim;
 836   Standard_Integer dim;
 837   Standard_Real u1,u2;
 838   Standard_Integer d1,d2;
 839   
 840   D = 1.0e0 ;
 841   Du = 0.0e0 ;
 842   Dv = 0.0e0 ;
 843   BSplSLib_DataContainer dc (UDegree, VDegree);
 844   Standard_Boolean ufirst = PrepareEval
 845     (U,V,UIndex,VIndex,UDegree,VDegree,URat,VRat,UPer,VPer,
 846      Poles,Weights,UKnots,VKnots,UMults,VMults,
 847      u1,u2,d1,d2,rational,dc);
 848   dim  = rational ? 4 : 3;
 849   
 850   BSplCLib::Bohm(u1,d1,1,*dc.knots1,dim * (d2 + 1),*dc.poles);
 851   BSplCLib::Bohm(u2,d2,1,*dc.knots2,dim,*dc.poles);
 852   BSplCLib::Eval(u2,d2,*dc.knots2,dim,*(dc.poles+dim*(d2+1)));
 853   
 854   Standard_Real *result, *resVu, *resVv;
 855   result = dc.poles;
 856   resVu  = result + (ufirst ? dim*(d2+1) : dim);
 857   resVv  = result + (ufirst ? dim : dim*(d2+1));
 858   
 859   N .SetX(result[0]);
 860   Nu.SetX(resVu [0]);
 861   Nv.SetX(resVv [0]);
 862   
 863   N .SetY(result[1]);
 864   Nu.SetY(resVu [1]);
 865   Nv.SetY(resVv [1]);
 866   
 867   N .SetZ(result[2]);
 868   Nu.SetZ(resVu [2]);
 869   Nv.SetZ(resVv [2]);
 870   
 871   if (rational) {
 872     D  = result[3];
 873     Du = resVu [3];
 874     Dv = resVv [3];
 875   }
 876 }
 877 
 878 //=======================================================================
 879 //function : D2
 880 //purpose  : 
 881 //=======================================================================
 882 
 883 void  BSplSLib::D2
 884 (const Standard_Real U,
 885  const Standard_Real V,
 886  const Standard_Integer UIndex,
 887  const Standard_Integer VIndex,
 888  const TColgp_Array2OfPnt& Poles,
 889  const TColStd_Array2OfReal* Weights,
 890  const TColStd_Array1OfReal& UKnots,
 891  const TColStd_Array1OfReal& VKnots,
 892  const TColStd_Array1OfInteger* UMults,
 893  const TColStd_Array1OfInteger* VMults,
 894  const Standard_Integer UDegree,
 895  const Standard_Integer VDegree,
 896  const Standard_Boolean URat,
 897  const Standard_Boolean VRat,
 898  const Standard_Boolean UPer,
 899  const Standard_Boolean VPer,
 900  gp_Pnt& P,
 901  gp_Vec& Vu,
 902  gp_Vec& Vv,
 903  gp_Vec& Vuu,
 904  gp_Vec& Vvv,
 905  gp_Vec& Vuv)
 906 {
 907   Standard_Boolean rational;
 908 //  Standard_Integer k,dim,dim2;
 909   Standard_Integer dim,dim2;
 910   Standard_Real u1,u2;
 911   Standard_Integer d1,d2;
 912   Standard_Real *result, *resVu, *resVv, *resVuu, *resVvv, *resVuv;
 913   BSplSLib_DataContainer dc (UDegree, VDegree);
 914   if (PrepareEval
 915       (U,V,UIndex,VIndex,UDegree,VDegree,URat,VRat,UPer,VPer,
 916        Poles,Weights,UKnots,VKnots,UMults,VMults,
 917        u1,u2,d1,d2,rational,dc)) {
 918     if (rational) {
 919       dim = 4;
 920       dim2 = (d2 + 1) << 2;
 921       BSplCLib::Bohm(u1,d1,2,*dc.knots1,dim2,*dc.poles);
 922       BSplCLib::Bohm(u2,d2,2,*dc.knots2,dim ,*dc.poles);
 923       BSplCLib::Bohm(u2,d2,1,*dc.knots2,dim ,*(dc.poles + dim2));
 924       if (d1 > 1)
 925         BSplCLib::Eval(u2,d2,*dc.knots2,dim ,*(dc.poles + (dim2 << 1)));
 926       BSplSLib::RationalDerivative(d1,d2,2,2,*dc.poles,*dc.ders);
 927       result = dc.ders;
 928       resVu  = result + 9;
 929       resVv  = result + 3;
 930       resVuu = result + 18;
 931       resVvv = result + 6;
 932       resVuv = result + 12;
 933     }
 934     else {
 935       dim = 3;
 936       dim2 = d2 + 1;
 937       dim2 = (dim2 << 1) + dim2;
 938       BSplCLib::Bohm(u1,d1,2,*dc.knots1,dim2,*dc.poles);
 939       BSplCLib::Bohm(u2,d2,2,*dc.knots2,dim ,*dc.poles);
 940       BSplCLib::Bohm(u2,d2,1,*dc.knots2,dim ,*(dc.poles + dim2));
 941       if (d1 > 1)
 942         BSplCLib::Eval(u2,d2,*dc.knots2,dim ,*(dc.poles + (dim2 << 1)));
 943       result = dc.poles;
 944       resVu  = result + dim2;
 945       resVv  = result + 3;
 946       if (UDegree <= 1) resVuu = BSplSLib_zero;
 947       else              resVuu = result + (dim2 << 1);
 948       if (VDegree <= 1) resVvv = BSplSLib_zero;
 949       else              resVvv = result + 6;
 950       resVuv = result + (d2 << 1) + d2 + 6; 
 951     }
 952   }
 953   else {
 954     if (rational) {
 955       dim = 4;
 956       dim2 = (d2 + 1) << 2;
 957       BSplCLib::Bohm(u1,d1,2,*dc.knots1,dim2,*dc.poles);
 958       BSplCLib::Bohm(u2,d2,2,*dc.knots2,dim ,*dc.poles);
 959       BSplCLib::Bohm(u2,d2,1,*dc.knots2,dim ,*(dc.poles + dim2));
 960       if (d1 > 1)
 961         BSplCLib::Eval(u2,d2,*dc.knots2,dim ,*(dc.poles + (dim2 << 1)));
 962       BSplSLib::RationalDerivative(d1,d2,2,2,*dc.poles,*dc.ders);
 963       result = dc.ders;
 964       resVu  = result + 3;
 965       resVv  = result + 9;
 966       resVuu = result + 6;
 967       resVvv = result + 18;
 968       resVuv = result + 12;
 969     }
 970     else {
 971       dim = 3;
 972       dim2 = d2 + 1;
 973       dim2 = (dim2 << 1) + dim2;
 974       BSplCLib::Bohm(u1,d1,2,*dc.knots1,dim2,*dc.poles);
 975       BSplCLib::Bohm(u2,d2,2,*dc.knots2,dim ,*dc.poles);
 976       BSplCLib::Bohm(u2,d2,1,*dc.knots2,dim ,*(dc.poles + dim2));
 977       if (d1 > 1)
 978         BSplCLib::Eval(u2,d2,*dc.knots2,dim ,*(dc.poles + (dim2 << 1)));
 979       result = dc.poles;
 980       resVu  = result + 3;
 981       resVv  = result + dim2;
 982       if (UDegree <= 1) resVuu = BSplSLib_zero;
 983       else              resVuu = result + 6;
 984       if (VDegree <= 1) resVvv = BSplSLib_zero;
 985       else              resVvv = result + (dim2 << 1);
 986       resVuv = result + (d2 << 1) + d2 + 6; 
 987     }
 988   }
 989 
 990   P  .SetX(result[0]);
 991   Vu .SetX(resVu [0]);
 992   Vv .SetX(resVv [0]);
 993   Vuu.SetX(resVuu[0]);
 994   Vvv.SetX(resVvv[0]);
 995   Vuv.SetX(resVuv[0]);
 996   
 997   P  .SetY(result[1]);
 998   Vu .SetY(resVu [1]);
 999   Vv .SetY(resVv [1]);
1000   Vuu.SetY(resVuu[1]);
1001   Vvv.SetY(resVvv[1]);
1002   Vuv.SetY(resVuv[1]);
1003   
1004   P  .SetZ(result[2]);
1005   Vu .SetZ(resVu [2]);
1006   Vv .SetZ(resVv [2]);
1007   Vuu.SetZ(resVuu[2]);
1008   Vvv.SetZ(resVvv[2]);
1009   Vuv.SetZ(resVuv[2]);
1010 }
1011 
1012 //=======================================================================
1013 //function : D3
1014 //purpose  : 
1015 //=======================================================================
1016 
1017 void  BSplSLib::D3
1018 (const Standard_Real U,
1019  const Standard_Real V,
1020  const Standard_Integer UIndex,
1021  const Standard_Integer VIndex,
1022  const TColgp_Array2OfPnt& Poles,
1023  const TColStd_Array2OfReal* Weights,
1024  const TColStd_Array1OfReal& UKnots,
1025  const TColStd_Array1OfReal& VKnots,
1026  const TColStd_Array1OfInteger* UMults,
1027  const TColStd_Array1OfInteger* VMults,
1028  const Standard_Integer UDegree,
1029  const Standard_Integer VDegree,
1030  const Standard_Boolean URat,
1031  const Standard_Boolean VRat,
1032  const Standard_Boolean UPer,
1033  const Standard_Boolean VPer,
1034  gp_Pnt& P,
1035  gp_Vec& Vu,
1036  gp_Vec& Vv,
1037  gp_Vec& Vuu,
1038  gp_Vec& Vvv,
1039  gp_Vec& Vuv,
1040  gp_Vec& Vuuu,
1041  gp_Vec& Vvvv,
1042  gp_Vec& Vuuv,
1043  gp_Vec& Vuvv)
1044 {
1045   Standard_Boolean rational;
1046 //  Standard_Integer k,dim,dim2;
1047   Standard_Integer dim,dim2;
1048   Standard_Real u1,u2;
1049   Standard_Integer d1,d2;
1050   Standard_Real *result, *resVu, *resVv, *resVuu, *resVvv, *resVuv,
1051   *resVuuu, *resVvvv, *resVuuv, *resVuvv;
1052   BSplSLib_DataContainer dc (UDegree, VDegree);
1053   if (PrepareEval
1054     (U,V,UIndex,VIndex,UDegree,VDegree,URat,VRat,UPer,VPer,
1055      Poles,Weights,UKnots,VKnots,UMults,VMults,
1056      u1,u2,d1,d2,rational,dc)) {
1057     if (rational) {
1058       dim = 4;
1059       dim2 = (d2 + 1) << 2;
1060       BSplCLib::Bohm  (u1,d1,3,*dc.knots1,dim2,*dc.poles);
1061       BSplCLib::Bohm  (u2,d2,3,*dc.knots2,dim ,*dc.poles);
1062       BSplCLib::Bohm  (u2,d2,2,*dc.knots2,dim ,*(dc.poles + dim2));
1063       if (d1 > 1)
1064         BSplCLib::Bohm(u2,d2,1,*dc.knots2,dim ,*(dc.poles + (dim2 << 1)));
1065       if (d1 > 2)
1066         BSplCLib::Eval(u2,d2  ,*dc.knots2,dim ,*(dc.poles + (dim2 << 1) + dim2));
1067       BSplSLib::RationalDerivative(d1,d2,3,3,*dc.poles,*dc.ders);
1068       result  = dc.ders;
1069       resVu   = result + 12;
1070       resVv   = result + 3;
1071       resVuu  = result + 24;
1072       resVvv  = result + 6;
1073       resVuv  = result + 15;
1074       resVuuu = result + 36;
1075       resVvvv = result + 9;
1076       resVuuv = result + 27;
1077       resVuvv = result + 18;
1078     }
1079     else {
1080       dim = 3;
1081       dim2 = (d2 + 1);
1082       dim2 = (dim2 << 1) + dim2;
1083       BSplCLib::Bohm  (u1,d1,3,*dc.knots1,dim2,*dc.poles);
1084       BSplCLib::Bohm  (u2,d2,3,*dc.knots2,dim ,*dc.poles);
1085       BSplCLib::Bohm  (u2,d2,2,*dc.knots2,dim ,*(dc.poles + dim2));
1086       if (d1 > 1)
1087         BSplCLib::Bohm(u2,d2,1,*dc.knots2,dim ,*(dc.poles + (dim2 << 1)));
1088       if (d1 > 2)
1089         BSplCLib::Eval(u2,d2  ,*dc.knots2,dim ,*(dc.poles + (dim2 << 1) + dim2));
1090       result = dc.poles;
1091       resVu  = result + dim2;
1092       resVv  = result + 3;
1093       if (UDegree <= 1) {
1094         resVuu  = BSplSLib_zero;
1095         resVuuv = BSplSLib_zero;
1096       }
1097       else {
1098         resVuu  = result + (dim2 << 1);
1099         resVuuv = result + (dim2 << 1) + 3;
1100       }
1101       if (VDegree <= 1) {
1102         resVvv  = BSplSLib_zero;
1103         resVuvv = BSplSLib_zero;
1104       }
1105       else {
1106         resVvv  = result + 6;
1107         resVuvv = result + dim2 + 6;
1108       }
1109       resVuv = result + (d2 << 1) + d2 + 6;
1110       if (UDegree <= 2) resVuuu = BSplSLib_zero;
1111       else              resVuuu = result + (dim2 << 1) + dim2;
1112       if (VDegree <= 2) resVvvv = BSplSLib_zero;
1113       else              resVvvv = result + 9;
1114     }
1115   }
1116   else {
1117     if (rational) {
1118       dim = 4;
1119       dim2 = (d2 + 1) << 2;
1120       BSplCLib::Bohm  (u1,d1,3,*dc.knots1,dim2,*dc.poles);
1121       BSplCLib::Bohm  (u2,d2,3,*dc.knots2,dim ,*dc.poles);
1122       BSplCLib::Bohm  (u2,d2,2,*dc.knots2,dim ,*(dc.poles + dim2));
1123       if (d1 > 1)
1124         BSplCLib::Bohm(u2,d2,1,*dc.knots2,dim ,*(dc.poles + (dim2 << 1)));
1125       if (d1 > 2)
1126         BSplCLib::Eval(u2,d2  ,*dc.knots2,dim ,*(dc.poles + (dim2 << 1) + dim2));
1127       BSplSLib::RationalDerivative(d1,d2,3,3,*dc.poles,*dc.ders);
1128       result  = dc.ders;
1129       resVu   = result + 3;
1130       resVv   = result + 12;
1131       resVuu  = result + 6;
1132       resVvv  = result + 24;
1133       resVuv  = result + 15;
1134       resVuuu = result + 9;
1135       resVvvv = result + 36;
1136       resVuuv = result + 18;
1137       resVuvv = result + 27;
1138     }
1139     else {
1140       dim = 3;
1141       dim2 = (d2 + 1);
1142       dim2 = (dim2 << 1) + dim2;
1143       BSplCLib::Bohm  (u1,d1,3,*dc.knots1,dim2,*dc.poles);
1144       BSplCLib::Bohm  (u2,d2,3,*dc.knots2,dim ,*dc.poles);
1145       BSplCLib::Bohm  (u2,d2,2,*dc.knots2,dim ,*(dc.poles + dim2));
1146       if (d1 > 1)
1147         BSplCLib::Bohm(u2,d2,1,*dc.knots2,dim ,*(dc.poles + (dim2 << 1)));
1148       if (d1 > 2)
1149         BSplCLib::Eval(u2,d2  ,*dc.knots2,dim ,*(dc.poles + (dim2 << 1) + dim2));
1150       result = dc.poles;
1151       resVu  = result + 3;
1152       resVv  = result + dim2;
1153       if (UDegree <= 1) {
1154         resVuu  = BSplSLib_zero;
1155         resVuuv = BSplSLib_zero;
1156       }
1157       else {
1158         resVuu  = result + 6;
1159         resVuuv = result + dim2 + 6;
1160       }
1161       if (VDegree <= 1) {
1162         resVvv  = BSplSLib_zero;
1163         resVuvv = BSplSLib_zero;
1164       }
1165       else {
1166         resVvv  = result + (dim2 << 1);
1167         resVuvv = result + (dim2 << 1) + 3;
1168       }
1169       resVuv = result + (d2 << 1) + d2 + 6;
1170       if (UDegree <= 2) resVuuu = BSplSLib_zero;
1171       else              resVuuu = result + 9;
1172       if (VDegree <= 2) resVvvv = BSplSLib_zero;
1173       else              resVvvv = result + (dim2 << 1) + dim2;
1174     }
1175   }
1176   
1177   P   .SetX(result [0]);
1178   Vu  .SetX(resVu  [0]);
1179   Vv  .SetX(resVv  [0]);
1180   Vuu .SetX(resVuu [0]);
1181   Vvv .SetX(resVvv [0]);
1182   Vuv .SetX(resVuv [0]);
1183   Vuuu.SetX(resVuuu[0]);
1184   Vvvv.SetX(resVvvv[0]);
1185   Vuuv.SetX(resVuuv[0]);
1186   Vuvv.SetX(resVuvv[0]);
1187   
1188   P   .SetY(result [1]);
1189   Vu  .SetY(resVu  [1]);
1190   Vv  .SetY(resVv  [1]);
1191   Vuu .SetY(resVuu [1]);
1192   Vvv .SetY(resVvv [1]);
1193   Vuv .SetY(resVuv [1]);
1194   Vuuu.SetY(resVuuu[1]);
1195   Vvvv.SetY(resVvvv[1]);
1196   Vuuv.SetY(resVuuv[1]);
1197   Vuvv.SetY(resVuvv[1]);
1198   
1199   P   .SetZ(result [2]);
1200   Vu  .SetZ(resVu  [2]);
1201   Vv  .SetZ(resVv  [2]);
1202   Vuu .SetZ(resVuu [2]);
1203   Vvv .SetZ(resVvv [2]);
1204   Vuv .SetZ(resVuv [2]);
1205   Vuuu.SetZ(resVuuu[2]);
1206   Vvvv.SetZ(resVvvv[2]);
1207   Vuuv.SetZ(resVuuv[2]);
1208   Vuvv.SetZ(resVuvv[2]);
1209 }
1210 
1211 //=======================================================================
1212 //function : DN
1213 //purpose  : 
1214 //=======================================================================
1215 
1216 void  BSplSLib::DN
1217 (const Standard_Real U,
1218  const Standard_Real V,
1219  const Standard_Integer Nu,
1220  const Standard_Integer Nv,
1221  const Standard_Integer UIndex,
1222  const Standard_Integer VIndex,
1223  const TColgp_Array2OfPnt& Poles,
1224  const TColStd_Array2OfReal* Weights,
1225  const TColStd_Array1OfReal& UKnots,
1226  const TColStd_Array1OfReal& VKnots,
1227  const TColStd_Array1OfInteger* UMults,
1228  const TColStd_Array1OfInteger* VMults,
1229  const Standard_Integer UDegree,
1230  const Standard_Integer VDegree,
1231  const Standard_Boolean URat,
1232  const Standard_Boolean VRat,
1233  const Standard_Boolean UPer,
1234  const Standard_Boolean VPer,
1235  gp_Vec& Vn)
1236 {
1237   Standard_Boolean rational;
1238   Standard_Integer k,dim;
1239   Standard_Real u1,u2;
1240   Standard_Integer d1,d2;
1241   
1242   BSplSLib_DataContainer dc (UDegree, VDegree);
1243   Standard_Boolean ufirst = PrepareEval
1244     (U,V,UIndex,VIndex,UDegree,VDegree,URat,VRat,UPer,VPer,
1245      Poles,Weights,UKnots,VKnots,UMults,VMults,
1246      u1,u2,d1,d2,rational,dc);
1247   dim  = rational ? 4 : 3;
1248   
1249   if (!rational) {
1250     if ((Nu > UDegree) || (Nv > VDegree)) {
1251       Vn.SetX(0.);
1252       Vn.SetY(0.);
1253       Vn.SetZ(0.);
1254       return;
1255     }
1256   }
1257   
1258   Standard_Integer n1 = ufirst ? Nu : Nv;
1259   Standard_Integer n2 = ufirst ? Nv : Nu;
1260   
1261   BSplCLib::Bohm(u1,d1,n1,*dc.knots1,dim * (d2 + 1),*dc.poles);
1262 
1263   for (k = 0; k <= Min(n1,d1); k++) 
1264     BSplCLib::Bohm(u2,d2,n2,*dc.knots2,dim,*(dc.poles+k*dim*(d2+1)));
1265   
1266   Standard_Real *result;
1267   if (rational) {
1268     BSplSLib::RationalDerivative(d1,d2,n1,n2,*dc.poles,*dc.ders,Standard_False);
1269     result = dc.ders; // because Standard_False ci-dessus.
1270     
1271   }
1272   else {
1273     result = dc.poles + (n1 * (d2+1) + n2) * dim ;
1274   }
1275   
1276   Vn.SetX(result[0]);
1277   Vn.SetY(result[1]);
1278   Vn.SetZ(result[2]);
1279 }
1280 
1281 //
1282 // Surface modifications
1283 // 
1284 // a surface is processed as a curve of curves.
1285 // i.e : a curve of parameter u where the current point is the set of poles 
1286 //       of the iso.
1287 //
1288 
1289 //=======================================================================
1290 //function : Iso
1291 //purpose  : 
1292 //=======================================================================
1293 
1294 void  BSplSLib::Iso(const Standard_Real            Param, 
1295                     const Standard_Boolean         IsU,
1296                     const TColgp_Array2OfPnt&      Poles,
1297                     const TColStd_Array2OfReal*    Weights,
1298                     const TColStd_Array1OfReal&    Knots,
1299                     const TColStd_Array1OfInteger* Mults,
1300                     const Standard_Integer         Degree,
1301                     const Standard_Boolean         Periodic,
1302                     TColgp_Array1OfPnt&            CPoles,
1303                     TColStd_Array1OfReal*          CWeights)
1304 {
1305   Standard_Integer index = 0;
1306   Standard_Real    u = Param;
1307   Standard_Boolean rational = Weights != NULL;
1308   Standard_Integer dim = rational ? 4 : 3;
1309   
1310   // compute local knots
1311   
1312   NCollection_LocalArray<Standard_Real> locknots1 (2*Degree);  
1313   BSplCLib::LocateParameter(Degree,Knots,Mults,u,Periodic,index,u);
1314   BSplCLib::BuildKnots(Degree,index,Periodic,Knots,Mults,*locknots1);
1315   if (Mults == NULL)
1316     index -= Knots.Lower() + Degree;
1317   else
1318     index = BSplCLib::PoleIndex(Degree,index,Periodic,*Mults);
1319   
1320   
1321   // copy the local poles
1322   
1323 //  Standard_Integer f1,l1,f2,l2,i,j,k;
1324   Standard_Integer f1,l1,f2,l2,i,j;
1325   
1326   if (IsU) {
1327     f1 = Poles.LowerRow();
1328     l1 = Poles.UpperRow();
1329     f2 = Poles.LowerCol();
1330     l2 = Poles.UpperCol();
1331   }
1332   else {
1333     f1 = Poles.LowerCol();
1334     l1 = Poles.UpperCol();
1335     f2 = Poles.LowerRow();
1336     l2 = Poles.UpperRow();
1337   }
1338   
1339   NCollection_LocalArray<Standard_Real> locpoles ((Degree+1) * (l2-f2+1) * dim);
1340   
1341   Standard_Real w, *pole = locpoles;
1342   index += f1;
1343 
1344   for (i = 0; i <= Degree; i++) {
1345 
1346     for (j = f2; j <= l2; j++) {
1347       
1348       const gp_Pnt& P  = IsU ? Poles(index,j)   : Poles(j,index);
1349       if (rational) { 
1350         pole[3] = w      = IsU ? (*Weights)(index,j) : (*Weights)(j,index);
1351         pole[0] = P.X() * w;
1352         pole[1] = P.Y() * w;
1353         pole[2] = P.Z() * w;
1354       }
1355       else {
1356         pole[0] = P.X();
1357         pole[1] = P.Y();
1358         pole[2] = P.Z();
1359       }
1360       pole += dim;
1361     }
1362     index++;
1363     if (index > l1) index = f1;
1364   }
1365   
1366   // compute the iso
1367   BSplCLib::Eval(u,Degree,*locknots1,(l2-f2+1)*dim,*locpoles);
1368   
1369   // get the result
1370   pole = locpoles;
1371 
1372   for (i = CPoles.Lower(); i <= CPoles.Upper(); i++) {
1373     gp_Pnt& P = CPoles(i);
1374     if (rational) {
1375       (*CWeights)(i) = w = pole[3];
1376       P.SetX( pole[0] / w);
1377       P.SetY( pole[1] / w);
1378       P.SetZ( pole[2] / w);
1379     }
1380     else {
1381       P.SetX( pole[0]);
1382       P.SetY( pole[1]);
1383       P.SetZ( pole[2]);
1384     }
1385     pole += dim;
1386   }
1387   
1388   // if the input is not rational but weights are wanted
1389   if (!rational && (CWeights != NULL)) {
1390 
1391     for (i = CWeights->Lower(); i <= CWeights->Upper(); i++)
1392       (*CWeights)(i) = 1.;
1393   }
1394 }
1395 
1396 //=======================================================================
1397 //function : Reverse
1398 //purpose  : 
1399 //=======================================================================
1400 
1401 void  BSplSLib::Reverse(      TColgp_Array2OfPnt& Poles,
1402                         const Standard_Integer    Last, 
1403                         const Standard_Boolean    UDirection)
1404 {
1405   Standard_Integer i,j, l = Last;
1406   if ( UDirection) {
1407     l = Poles.LowerRow() + 
1408       (l - Poles.LowerRow())%(Poles.ColLength());
1409     TColgp_Array2OfPnt temp(0, Poles.ColLength()-1,
1410                             Poles.LowerCol(), Poles.UpperCol());
1411     
1412     for (i = Poles.LowerRow(); i <= l; i++) {
1413 
1414       for (j = Poles.LowerCol(); j <= Poles.UpperCol(); j++) {
1415         temp(l-i,j) = Poles(i,j);
1416       }
1417     }
1418 
1419     for (i = l+1; i <= Poles.UpperRow(); i++) {
1420 
1421       for (j = Poles.LowerCol(); j <= Poles.UpperCol(); j++) {
1422         temp(l+Poles.ColLength()-i,j) = Poles(i,j);
1423       }
1424     }
1425 
1426     for (i = Poles.LowerRow(); i <= Poles.UpperRow(); i++) {
1427 
1428       for (j = Poles.LowerCol(); j <= Poles.UpperCol(); j++) {
1429         Poles(i,j) = temp (i-Poles.LowerRow(),j);
1430       }
1431     }
1432   }
1433   else {
1434     l = Poles.LowerCol() + 
1435       (l - Poles.LowerCol())%(Poles.RowLength());
1436     TColgp_Array2OfPnt temp(Poles.LowerRow(), Poles.UpperRow(),
1437                             0, Poles.RowLength()-1);
1438 
1439     for (j = Poles.LowerCol(); j <= l; j++) {
1440 
1441       for (i = Poles.LowerRow(); i <= Poles.UpperRow(); i++) {
1442         temp(i,l-j) = Poles(i,j);
1443       }
1444     }
1445 
1446     for (j = l+1; j <= Poles.UpperCol(); j++) {
1447 
1448       for (i = Poles.LowerRow(); i <= Poles.UpperRow(); i++) {
1449         temp(i,l+Poles.RowLength()-j) = Poles(i,j);
1450       }
1451     }
1452 
1453     for (i = Poles.LowerRow(); i <= Poles.UpperRow(); i++) {
1454 
1455       for (j = Poles.LowerCol(); j <= Poles.UpperCol(); j++) {
1456         Poles(i,j) = temp (i,j-Poles.LowerCol());
1457       }
1458     }
1459   }
1460 }
1461 
1462 //=======================================================================
1463 //function : Reverse
1464 //purpose  : 
1465 //=======================================================================
1466 
1467 void  BSplSLib::Reverse(      TColStd_Array2OfReal& Weights, 
1468                         const Standard_Integer      Last, 
1469                         const Standard_Boolean      UDirection)
1470 {
1471   Standard_Integer i,j, l = Last;
1472   if ( UDirection) {
1473     l = Weights.LowerRow() + 
1474       (l - Weights.LowerRow())%(Weights.ColLength());
1475     TColStd_Array2OfReal temp(0, Weights.ColLength()-1,
1476                               Weights.LowerCol(), Weights.UpperCol());
1477 
1478     for (i = Weights.LowerRow(); i <= l; i++) {
1479 
1480       for (j = Weights.LowerCol(); j <= Weights.UpperCol(); j++) {
1481         temp(l-i,j) = Weights(i,j);
1482       }
1483     }
1484 
1485     for (i = l+1; i <= Weights.UpperRow(); i++) {
1486 
1487       for (j = Weights.LowerCol(); j <= Weights.UpperCol(); j++) {
1488         temp(l+Weights.ColLength()-i,j) = Weights(i,j);
1489       }
1490     }
1491 
1492     for (i = Weights.LowerRow(); i <= Weights.UpperRow(); i++) {
1493 
1494       for (j = Weights.LowerCol(); j <= Weights.UpperCol(); j++) {
1495         Weights(i,j) = temp (i-Weights.LowerRow(),j);
1496       }
1497     }
1498   }
1499   else {
1500     l = Weights.LowerCol() + 
1501       (l - Weights.LowerCol())%(Weights.RowLength());
1502     TColStd_Array2OfReal temp(Weights.LowerRow(), Weights.UpperRow(),
1503                               0, Weights.RowLength()-1);
1504 
1505     for (j = Weights.LowerCol(); j <= l; j++) {
1506 
1507       for (i = Weights.LowerRow(); i <= Weights.UpperRow(); i++) {
1508         temp(i,l-j) = Weights(i,j);
1509       }
1510     }
1511 
1512     for (j = l+1; j <= Weights.UpperCol(); j++) {
1513 
1514       for (i = Weights.LowerRow(); i <= Weights.UpperRow(); i++) {
1515         temp(i,l+Weights.RowLength()-j) = Weights(i,j);
1516       }
1517     }
1518 
1519     for (i = Weights.LowerRow(); i <= Weights.UpperRow(); i++) {
1520 
1521       for (j = Weights.LowerCol(); j <= Weights.UpperCol(); j++) {
1522         Weights(i,j) = temp (i,j-Weights.LowerCol());
1523       }
1524     }
1525   }
1526 }
1527 
1528 //=======================================================================
1529 //function : IsRational
1530 //purpose  : 
1531 //=======================================================================
1532 
1533 Standard_Boolean BSplSLib::IsRational
1534 (const TColStd_Array2OfReal& Weights, 
1535  const Standard_Integer      I1, 
1536  const Standard_Integer      I2, 
1537  const Standard_Integer      J1, 
1538  const Standard_Integer      J2, 
1539  const Standard_Real         Epsi)
1540 {
1541   Standard_Real eps = (Epsi > 0.0) ? Epsi : Epsilon(Weights(I1,I2));
1542   Standard_Integer i,j;
1543   Standard_Integer fi = Weights.LowerRow(), li = Weights.ColLength();
1544   Standard_Integer fj = Weights.LowerCol(), lj = Weights.RowLength();
1545 
1546   for (i = I1 - fi; i < I2 - fi; i++) {
1547 
1548     for (j = J1 - fj; j < J2 - fj; j++) {
1549       if (Abs(Weights(fi+i%li,fj+j%lj)-Weights(fi+(i+1)%li,fj+j%lj))>eps)
1550         return Standard_True;
1551     }
1552   }
1553   return Standard_False;
1554 }
1555 
1556 //=======================================================================
1557 //function : SetPoles
1558 //purpose  : 
1559 //=======================================================================
1560 
1561 void  BSplSLib::SetPoles(const TColgp_Array2OfPnt&   Poles, 
1562                          TColStd_Array1OfReal& FP,
1563                          const Standard_Boolean      UDirection)
1564 {
1565   Standard_Integer i, j, l = FP.Lower();
1566   Standard_Integer PLowerRow = Poles.LowerRow();
1567   Standard_Integer PUpperRow = Poles.UpperRow();
1568   Standard_Integer PLowerCol = Poles.LowerCol();
1569   Standard_Integer PUpperCol = Poles.UpperCol();
1570   if (UDirection) {
1571 
1572     for ( i = PLowerRow; i <= PUpperRow; i++) {
1573 
1574       for ( j = PLowerCol; j <= PUpperCol; j++) {
1575         const gp_Pnt& P = Poles.Value(i,j);
1576         FP(l) = P.X(); l++;
1577         FP(l) = P.Y(); l++;
1578         FP(l) = P.Z(); l++;
1579       }
1580     }
1581   }
1582   else {
1583 
1584     for ( j = PLowerCol; j <= PUpperCol; j++) {
1585 
1586       for ( i = PLowerRow; i <= PUpperRow; i++) {
1587         const gp_Pnt& P = Poles.Value(i,j);
1588         FP(l) = P.X(); l++;
1589         FP(l) = P.Y(); l++;
1590         FP(l) = P.Z(); l++;
1591       }
1592     }
1593   }
1594 }
1595 
1596 //=======================================================================
1597 //function : SetPoles
1598 //purpose  : 
1599 //=======================================================================
1600 
1601 void  BSplSLib::SetPoles(const TColgp_Array2OfPnt&   Poles, 
1602                          const TColStd_Array2OfReal& Weights, 
1603                          TColStd_Array1OfReal& FP,
1604                          const Standard_Boolean      UDirection)
1605 {
1606   Standard_Integer i, j, l = FP.Lower();
1607   Standard_Integer PLowerRow = Poles.LowerRow();
1608   Standard_Integer PUpperRow = Poles.UpperRow();
1609   Standard_Integer PLowerCol = Poles.LowerCol();
1610   Standard_Integer PUpperCol = Poles.UpperCol();
1611   if (UDirection) {
1612 
1613     for ( i = PLowerRow; i <= PUpperRow; i++) {
1614 
1615       for ( j = PLowerCol; j <= PUpperCol; j++) {
1616         const gp_Pnt& P = Poles  .Value(i,j);
1617         Standard_Real w = Weights.Value(i,j);
1618         FP(l) = P.X() * w; l++;
1619         FP(l) = P.Y() * w; l++;
1620         FP(l) = P.Z() * w; l++;
1621         FP(l) = w; l++;
1622       }
1623     }
1624   }
1625   else {
1626 
1627     for ( j = PLowerCol; j <= PUpperCol; j++) {
1628 
1629       for ( i = PLowerRow; i <= PUpperRow; i++) {
1630         const gp_Pnt& P = Poles  .Value(i,j);
1631         Standard_Real w = Weights.Value(i,j);
1632         FP(l) = P.X() * w; l++;
1633         FP(l) = P.Y() * w; l++;
1634         FP(l) = P.Z() * w; l++;
1635         FP(l) = w; l++;
1636       }
1637     }
1638   }
1639 }
1640 
1641 //=======================================================================
1642 //function : GetPoles
1643 //purpose  : 
1644 //=======================================================================
1645 
1646 void  BSplSLib::GetPoles(const TColStd_Array1OfReal& FP, 
1647                          TColgp_Array2OfPnt&   Poles,
1648                          const Standard_Boolean      UDirection)
1649 {
1650   Standard_Integer i, j, l = FP.Lower();
1651   Standard_Integer PLowerRow = Poles.LowerRow();
1652   Standard_Integer PUpperRow = Poles.UpperRow();
1653   Standard_Integer PLowerCol = Poles.LowerCol();
1654   Standard_Integer PUpperCol = Poles.UpperCol();
1655   if (UDirection) {
1656 
1657     for ( i = PLowerRow; i <= PUpperRow; i++) {
1658 
1659       for ( j = PLowerCol; j <= PUpperCol; j++) {
1660         gp_Pnt& P = Poles.ChangeValue(i,j);
1661         P.SetX(FP(l)); l++;
1662         P.SetY(FP(l)); l++;
1663         P.SetZ(FP(l)); l++;
1664       }
1665     }
1666   }
1667   else {
1668 
1669     for ( j = PLowerCol; j <= PUpperCol; j++) {
1670 
1671       for ( i = PLowerRow; i <= PUpperRow; i++) {
1672         gp_Pnt& P = Poles.ChangeValue(i,j);
1673         P.SetX(FP(l)); l++;
1674         P.SetY(FP(l)); l++;
1675         P.SetZ(FP(l)); l++;
1676       }
1677     }
1678   }
1679 }
1680 
1681 //=======================================================================
1682 //function : GetPoles
1683 //purpose  : 
1684 //=======================================================================
1685 
1686 void  BSplSLib::GetPoles(const TColStd_Array1OfReal& FP, 
1687                          TColgp_Array2OfPnt&   Poles, 
1688                          TColStd_Array2OfReal& Weights,
1689                          const Standard_Boolean      UDirection)
1690 {
1691   Standard_Integer i, j, l = FP.Lower();
1692   Standard_Integer PLowerRow = Poles.LowerRow();
1693   Standard_Integer PUpperRow = Poles.UpperRow();
1694   Standard_Integer PLowerCol = Poles.LowerCol();
1695   Standard_Integer PUpperCol = Poles.UpperCol();
1696   if (UDirection) {
1697 
1698     for ( i = PLowerRow; i <= PUpperRow; i++) {
1699 
1700       for ( j = PLowerCol; j <= PUpperCol; j++) {
1701         Standard_Real w = FP( l + 3);
1702         Weights(i,j) = w;
1703         gp_Pnt& P = Poles.ChangeValue(i,j);
1704         P.SetX(FP(l) / w); l++;
1705         P.SetY(FP(l) / w); l++;
1706         P.SetZ(FP(l) / w); l++;
1707         l++;
1708       }
1709     }
1710   }
1711   else {
1712 
1713     for ( j = PLowerCol; j <= PUpperCol; j++) {
1714 
1715       for ( i = PLowerRow; i <= PUpperRow; i++) {
1716         Standard_Real w = FP( l + 3);
1717         Weights(i,j) = w;
1718         gp_Pnt& P = Poles.ChangeValue(i,j);
1719         P.SetX(FP(l) / w); l++;
1720         P.SetY(FP(l) / w); l++;
1721         P.SetZ(FP(l) / w); l++;
1722         l++;
1723       }
1724     }
1725   }
1726 }
1727 
1728 //=======================================================================
1729 //function : InsertKnots
1730 //purpose  : 
1731 //=======================================================================
1732 
1733 void  BSplSLib::InsertKnots(const Standard_Boolean         UDirection, 
1734                             const Standard_Integer         Degree, 
1735                             const Standard_Boolean         Periodic, 
1736                             const TColgp_Array2OfPnt&      Poles, 
1737                             const TColStd_Array2OfReal*    Weights, 
1738                             const TColStd_Array1OfReal&    Knots, 
1739                             const TColStd_Array1OfInteger& Mults, 
1740                             const TColStd_Array1OfReal&    AddKnots, 
1741                             const TColStd_Array1OfInteger* AddMults, 
1742                             TColgp_Array2OfPnt&      NewPoles,
1743                             TColStd_Array2OfReal*    NewWeights,
1744                             TColStd_Array1OfReal&    NewKnots,
1745                             TColStd_Array1OfInteger& NewMults, 
1746                             const Standard_Real            Epsilon, 
1747                             const Standard_Boolean         Add )
1748 {
1749   Standard_Boolean rational = Weights != NULL;
1750   Standard_Integer dim = 3;
1751   if (rational) dim++;
1752   
1753   TColStd_Array1OfReal poles( 1, dim*Poles.RowLength()*Poles.ColLength());
1754   TColStd_Array1OfReal 
1755     newpoles( 1, dim*NewPoles.RowLength()*NewPoles.ColLength());
1756   
1757   if (rational) SetPoles(Poles,*Weights,poles,UDirection);
1758   else          SetPoles(Poles,poles,UDirection);
1759   
1760   if (UDirection) {
1761     dim *= Poles.RowLength();
1762   }
1763   else {
1764     dim *= Poles.ColLength();
1765   }
1766   BSplCLib::InsertKnots(Degree,Periodic,dim,poles,Knots,Mults,
1767                         AddKnots,AddMults,newpoles,NewKnots,NewMults,
1768                         Epsilon,Add);
1769   
1770   if (rational) GetPoles(newpoles,NewPoles,*NewWeights,UDirection);
1771   else          GetPoles(newpoles,NewPoles,UDirection);
1772 }
1773 
1774 //=======================================================================
1775 //function : RemoveKnot
1776 //purpose  : 
1777 //=======================================================================
1778 
1779 Standard_Boolean  BSplSLib::RemoveKnot
1780 (const Standard_Boolean         UDirection, 
1781  const Standard_Integer         Index,
1782  const Standard_Integer         Mult,
1783  const Standard_Integer         Degree, 
1784  const Standard_Boolean         Periodic, 
1785  const TColgp_Array2OfPnt&      Poles, 
1786  const TColStd_Array2OfReal*    Weights, 
1787  const TColStd_Array1OfReal&    Knots, 
1788  const TColStd_Array1OfInteger& Mults,
1789  TColgp_Array2OfPnt&      NewPoles,
1790  TColStd_Array2OfReal*    NewWeights,
1791  TColStd_Array1OfReal&    NewKnots,
1792  TColStd_Array1OfInteger& NewMults,
1793  const Standard_Real            Tolerance)
1794 {
1795   Standard_Boolean rational = Weights != NULL;
1796   Standard_Integer dim = 3;
1797   if (rational) dim++;
1798   
1799   TColStd_Array1OfReal poles( 1, dim*Poles.RowLength()*Poles.ColLength());
1800   TColStd_Array1OfReal 
1801     newpoles( 1, dim*NewPoles.RowLength()*NewPoles.ColLength());
1802   
1803   if (rational) SetPoles(Poles,*Weights,poles,UDirection);
1804   else          SetPoles(Poles,poles,UDirection);
1805   
1806   if (UDirection) {
1807     dim *= Poles.RowLength();
1808   }
1809   else {
1810     dim *= Poles.ColLength();
1811   }
1812   
1813   if ( !BSplCLib::RemoveKnot(Index,Mult,Degree,Periodic,dim,
1814                              poles,Knots,Mults,newpoles,NewKnots,NewMults,
1815                              Tolerance))
1816     return Standard_False;
1817   
1818   if (rational) GetPoles(newpoles,NewPoles,*NewWeights,UDirection);
1819   else          GetPoles(newpoles,NewPoles,UDirection);
1820   return Standard_True;
1821 }
1822 
1823 //=======================================================================
1824 //function : IncreaseDegree
1825 //purpose  : 
1826 //=======================================================================
1827 
1828 void  BSplSLib::IncreaseDegree
1829 (const Standard_Boolean         UDirection,
1830  const Standard_Integer         Degree, 
1831  const Standard_Integer         NewDegree, 
1832  const Standard_Boolean         Periodic, 
1833  const TColgp_Array2OfPnt&      Poles, 
1834  const TColStd_Array2OfReal*    Weights,
1835  const TColStd_Array1OfReal&    Knots,
1836  const TColStd_Array1OfInteger& Mults, 
1837  TColgp_Array2OfPnt&      NewPoles, 
1838  TColStd_Array2OfReal*    NewWeights, 
1839  TColStd_Array1OfReal&    NewKnots, 
1840  TColStd_Array1OfInteger& NewMults)
1841 {
1842   Standard_Boolean rational = Weights != NULL;
1843   Standard_Integer dim = 3;
1844   if (rational) dim++;
1845   
1846   TColStd_Array1OfReal poles( 1, dim*Poles.RowLength()*Poles.ColLength());
1847   TColStd_Array1OfReal 
1848     newpoles( 1, dim*NewPoles.RowLength()*NewPoles.ColLength());
1849   
1850   if (rational) SetPoles(Poles,*Weights,poles,UDirection);
1851   else          SetPoles(Poles,poles,UDirection);
1852   
1853   if (UDirection) {
1854     dim *= Poles.RowLength();
1855   }
1856   else {
1857     dim *= Poles.ColLength();
1858   }
1859   
1860   BSplCLib::IncreaseDegree(Degree,NewDegree,Periodic,dim,poles,Knots,Mults,
1861                            newpoles,NewKnots,NewMults);
1862   
1863   if (rational) GetPoles(newpoles,NewPoles,*NewWeights,UDirection);
1864   else          GetPoles(newpoles,NewPoles,UDirection);
1865 }
1866 
1867 //=======================================================================
1868 //function : Unperiodize
1869 //purpose  : 
1870 //=======================================================================
1871 
1872 void  BSplSLib::Unperiodize
1873 (const Standard_Boolean         UDirection,
1874  const Standard_Integer         Degree,
1875  const TColStd_Array1OfInteger& Mults,
1876  const TColStd_Array1OfReal&    Knots,
1877  const TColgp_Array2OfPnt&      Poles, 
1878  const TColStd_Array2OfReal*    Weights,
1879  TColStd_Array1OfInteger& NewMults,
1880  TColStd_Array1OfReal&    NewKnots,
1881  TColgp_Array2OfPnt&      NewPoles,
1882  TColStd_Array2OfReal*    NewWeights)
1883 {
1884   Standard_Boolean rational = Weights != NULL;
1885   Standard_Integer dim = 3;
1886   if (rational) dim++;
1887   
1888   TColStd_Array1OfReal poles( 1, dim*Poles.RowLength()*Poles.ColLength());
1889   TColStd_Array1OfReal 
1890     newpoles( 1, dim*NewPoles.RowLength()*NewPoles.ColLength());
1891   
1892   if (rational) SetPoles(Poles,*Weights,poles,UDirection);
1893   else          SetPoles(Poles,poles,UDirection);
1894   
1895   if (UDirection) {
1896     dim *= Poles.RowLength();
1897   }
1898   else {
1899     dim *= Poles.ColLength();
1900   }
1901   
1902   BSplCLib::Unperiodize(Degree, dim, Mults, Knots, poles, 
1903                         NewMults, NewKnots, newpoles);
1904   
1905   if (rational) GetPoles(newpoles,NewPoles,*NewWeights,UDirection);
1906   else          GetPoles(newpoles,NewPoles,UDirection);
1907 }
1908 
1909 //=======================================================================
1910 //function : BuildCache
1911 //purpose  : Stores theTaylor expansion normalized between 0,1 in the
1912 //           Cache : in case of  a rational function the Poles are
1913 //           stored in homogeneous form 
1914 //=======================================================================
1915 
1916 void BSplSLib::BuildCache
1917 (const Standard_Real            U,   
1918  const Standard_Real            V,
1919  const Standard_Real            USpanDomain,
1920  const Standard_Real            VSpanDomain,  
1921  const Standard_Boolean         UPeriodic,
1922  const Standard_Boolean         VPeriodic,
1923  const Standard_Integer         UDegree,
1924  const Standard_Integer         VDegree,
1925  const Standard_Integer         UIndex,
1926  const Standard_Integer         VIndex,
1927  const TColStd_Array1OfReal&    UFlatKnots,   
1928  const TColStd_Array1OfReal&    VFlatKnots,   
1929  const TColgp_Array2OfPnt&      Poles,  
1930  const TColStd_Array2OfReal*    Weights,
1931  TColgp_Array2OfPnt&            CachePoles,
1932  TColStd_Array2OfReal*          CacheWeights)
1933 {  
1934   Standard_Boolean rational,rational_u,rational_v,flag_u_or_v;                  
1935   Standard_Integer kk,d1,d1p1,d2,d2p1,ii,jj,iii,jjj,Index;
1936   Standard_Real u1,min_degree_domain,max_degree_domain,f,factor[2],u2;
1937   if (Weights != NULL) 
1938     rational_u = rational_v = Standard_True;
1939   else
1940     rational_u = rational_v = Standard_False;
1941   BSplSLib_DataContainer dc (UDegree, VDegree);
1942   flag_u_or_v =
1943     PrepareEval  (U,
1944                   V,
1945                   UIndex,
1946                   VIndex,
1947                   UDegree,
1948                   VDegree,
1949                   rational_u,
1950                   rational_v,
1951                   UPeriodic,
1952                   VPeriodic,
1953                   Poles,
1954                   Weights,
1955                   UFlatKnots,
1956                   VFlatKnots,
1957                   (BSplCLib::NoMults()),
1958                   (BSplCLib::NoMults()),
1959                   u1,
1960                   u2,
1961                   d1,
1962                   d2,
1963                   rational,
1964           dc);
1965   d1p1 = d1 + 1;
1966   d2p1 = d2 + 1;
1967   if (rational) {
1968     BSplCLib::Bohm(u1,d1,d1,*dc.knots1,4 * d2p1,*dc.poles);
1969     
1970     for (kk = 0; kk <= d1 ; kk++) 
1971       BSplCLib::Bohm(u2,d2,d2,*dc.knots2,4,*(dc.poles + kk * 4 * d2p1));
1972     if (flag_u_or_v) {
1973       min_degree_domain = USpanDomain ;
1974       max_degree_domain = VSpanDomain ;
1975     }
1976     else {
1977       min_degree_domain = VSpanDomain ;
1978       max_degree_domain = USpanDomain ;
1979     }
1980     factor[0] = 1.0e0 ;
1981     
1982     for (ii = 0 ; ii <= d2 ; ii++) {
1983       iii = ii + 1;
1984       factor[1] = 1.0e0 ;
1985       
1986       for (jj = 0 ; jj <= d1 ; jj++) {
1987         jjj = jj + 1;
1988         Index = jj * d2p1 + ii ;
1989         Index = Index << 2;
1990         gp_Pnt& P = CachePoles(iii,jjj);
1991         f = factor[0] * factor[1];
1992         P.SetX( f * dc.poles[Index]); Index++;
1993         P.SetY( f * dc.poles[Index]); Index++;
1994         P.SetZ( f * dc.poles[Index]); Index++;
1995         (*CacheWeights)(iii ,jjj) = f * dc.poles[Index] ;
1996         factor[1] *= min_degree_domain / (Standard_Real) (jjj) ;
1997       }
1998       factor[0] *= max_degree_domain / (Standard_Real) (iii) ;
1999     }
2000   }
2001   else {
2002     BSplCLib::Bohm(u1,d1,d1,*dc.knots1,3 * d2p1,*dc.poles);
2003     
2004     for (kk = 0; kk <= d1 ; kk++) 
2005       BSplCLib::Bohm(u2,d2,d2,*dc.knots2,3,*(dc.poles + kk * 3 * d2p1));
2006     if (flag_u_or_v) {
2007       min_degree_domain = USpanDomain ;
2008       max_degree_domain = VSpanDomain ;
2009     }
2010     else {
2011       min_degree_domain = VSpanDomain ;
2012       max_degree_domain = USpanDomain ;
2013     }
2014     factor[0] = 1.0e0 ;
2015     
2016     for (ii = 0 ; ii <= d2 ; ii++) {
2017       iii = ii + 1;
2018       factor[1] = 1.0e0 ;
2019       
2020       for (jj = 0 ; jj <= d1 ; jj++) {
2021         jjj = jj + 1;
2022         Index = jj * d2p1 + ii ;
2023         Index = (Index << 1) + Index;
2024         gp_Pnt& P = CachePoles(iii,jjj);
2025         f = factor[0] * factor[1];
2026         P.SetX( f * dc.poles[Index]); Index++;
2027         P.SetY( f * dc.poles[Index]); Index++;
2028         P.SetZ( f * dc.poles[Index]);
2029         factor[1] *= min_degree_domain / (Standard_Real) (jjj) ;
2030       }
2031       factor[0] *= max_degree_domain / (Standard_Real) (iii) ;
2032     }
2033     if (Weights != NULL) {
2034       //
2035       // means that PrepareEval did found out that the surface was 
2036       // locally polynomial but since the surface is constructed
2037       // with some weights we need to set the weight polynomial to constant
2038       // 
2039       
2040       for (ii = 1 ; ii <= d2p1 ; ii++) {
2041         
2042         for (jj = 1 ; jj <= d1p1 ; jj++) {
2043           (*CacheWeights)(ii,jj) = 0.0e0 ;
2044         }
2045       }
2046       (*CacheWeights)(1,1) = 1.0e0 ;
2047     }
2048   }
2049 }
2050 
2051 void BSplSLib::BuildCache(const Standard_Real          theU,
2052                           const Standard_Real          theV,
2053                           const Standard_Real          theUSpanDomain,
2054                           const Standard_Real          theVSpanDomain,
2055                           const Standard_Boolean       theUPeriodicFlag,
2056                           const Standard_Boolean       theVPeriodicFlag,
2057                           const Standard_Integer       theUDegree,
2058                           const Standard_Integer       theVDegree,
2059                           const Standard_Integer       theUIndex,
2060                           const Standard_Integer       theVIndex,
2061                           const TColStd_Array1OfReal&  theUFlatKnots,
2062                           const TColStd_Array1OfReal&  theVFlatKnots,
2063                           const TColgp_Array2OfPnt&    thePoles,
2064                           const TColStd_Array2OfReal*  theWeights,
2065                                 TColStd_Array2OfReal&  theCacheArray)
2066 {
2067   Standard_Boolean flag_u_or_v;
2068   Standard_Integer d1, d2;
2069   Standard_Real    u1, u2;
2070   Standard_Boolean isRationalOnParam = (theWeights != NULL);
2071   Standard_Boolean isRational;
2072 
2073   BSplSLib_DataContainer dc(theUDegree, theVDegree);
2074   flag_u_or_v =
2075     PrepareEval(theU, theV, theUIndex, theVIndex, theUDegree, theVDegree,
2076                 isRationalOnParam, isRationalOnParam,
2077                 theUPeriodicFlag, theVPeriodicFlag,
2078                 thePoles, theWeights,
2079                 theUFlatKnots, theVFlatKnots,
2080                 (BSplCLib::NoMults()), (BSplCLib::NoMults()),
2081                 u1, u2, d1, d2, isRational, dc);
2082 
2083   Standard_Integer d2p1 = d2 + 1;
2084   Standard_Integer aDimension = isRational ? 4 : 3;
2085   Standard_Integer aCacheShift = // helps to store weights when PrepareEval did not found that the surface is locally polynomial
2086     (isRationalOnParam && !isRational) ? aDimension + 1 : aDimension;
2087 
2088   Standard_Real aDomains[2];
2089   // aDomains[0] corresponds to variable with minimal degree
2090   // aDomains[1] corresponds to variable with maximal degree
2091   if (flag_u_or_v)
2092   {
2093     aDomains[0] = theUSpanDomain;
2094     aDomains[1] = theVSpanDomain;
2095   }
2096   else
2097   {
2098     aDomains[0] = theVSpanDomain;
2099     aDomains[1] = theUSpanDomain;
2100   }
2101 
2102   BSplCLib::Bohm(u1, d1, d1, *dc.knots1, aDimension * d2p1, *dc.poles);
2103   for (Standard_Integer kk = 0; kk <= d1 ; kk++) 
2104     BSplCLib::Bohm(u2, d2, d2, *dc.knots2, aDimension, *(dc.poles + kk * aDimension * d2p1));
2105 
2106   Standard_Real* aCache = (Standard_Real *) &(theCacheArray(theCacheArray.LowerRow(), theCacheArray.LowerCol()));
2107   Standard_Real* aPolyCoeffs = dc.poles;
2108 
2109   Standard_Real aFactors[2];
2110   // aFactors[0] corresponds to variable with minimal degree
2111   // aFactors[1] corresponds to variable with maximal degree
2112   aFactors[1] = 1.0;
2113   Standard_Integer aRow, aCol, i;
2114   Standard_Real aCoeff;
2115   for (aRow = 0; aRow <= d2; aRow++)
2116   {
2117     aFactors[0] = 1.0;
2118     for (aCol = 0; aCol <= d1; aCol++)
2119     {
2120       aPolyCoeffs = dc.poles + (aCol * d2p1 + aRow) * aDimension;
2121       aCoeff = aFactors[0] * aFactors[1];
2122       for (i = 0; i < aDimension; i++)
2123         aCache[i] = aPolyCoeffs[i] * aCoeff;
2124       aCache += aCacheShift;
2125       aFactors[0] *= aDomains[0] / (aCol + 1);
2126     }
2127     aFactors[1] *= aDomains[1] / (aRow + 1);
2128   }
2129 
2130   // Fill the weights for the surface which is not locally polynomial
2131   if (aCacheShift > aDimension)
2132   {
2133     aCache = (Standard_Real *) &(theCacheArray(theCacheArray.LowerRow(), theCacheArray.LowerCol()));
2134     aCache += aCacheShift - 1;
2135     for (aRow = 0; aRow <= d2; aRow++)
2136       for (aCol = 0; aCol <= d1; aCol++)
2137       {
2138         *aCache = 0.0;
2139         aCache += aCacheShift;
2140       }
2141     theCacheArray.SetValue(theCacheArray.LowerRow(), theCacheArray.LowerCol() + aCacheShift - 1, 1.0);
2142   }
2143 }
2144 
2145 
2146 //=======================================================================
2147 //function : CacheD0
2148 //purpose  : Evaluates the polynomial cache of the Bspline Curve
2149 //           
2150 //=======================================================================
2151 
2152 void  BSplSLib::CacheD0(const Standard_Real                  UParameter,
2153                         const Standard_Real                  VParameter,
2154                         const  Standard_Integer              UDegree,
2155                         const  Standard_Integer              VDegree,
2156                         const  Standard_Real                 UCacheParameter,
2157                         const  Standard_Real                 VCacheParameter,
2158                         const  Standard_Real                 USpanLenght,
2159                         const  Standard_Real                 VSpanLenght,
2160                         const  TColgp_Array2OfPnt&           PolesArray,
2161                         const  TColStd_Array2OfReal*         WeightsArray,
2162                         gp_Pnt&                              aPoint)
2163 {
2164   //
2165   // the CacheParameter is where the cache polynomial was evaluated in homogeneous
2166   // form
2167   // the SpanLenght     is the normalizing factor so that the polynomial is between
2168   // 0 and 1 
2169   Standard_Integer 
2170 //    ii,
2171   dimension,
2172   min_degree,
2173   max_degree  ;
2174   
2175   Standard_Real 
2176     new_parameter[2] ,
2177   inverse ;
2178   
2179   Standard_Real * 
2180     PArray = (Standard_Real *) 
2181       &(PolesArray(PolesArray.LowerCol(), PolesArray.LowerRow())) ;
2182   Standard_Real *
2183     myPoint = (Standard_Real *) &aPoint  ;
2184   if (UDegree <= VDegree) {
2185     min_degree = UDegree ;
2186     max_degree = VDegree ;
2187     new_parameter[1] = (UParameter - UCacheParameter) / USpanLenght ;
2188     new_parameter[0] = (VParameter - VCacheParameter) / VSpanLenght ; 
2189     dimension = 3 * (UDegree + 1) ;
2190   }
2191   else {
2192     min_degree = VDegree ;
2193     max_degree = UDegree ;
2194     new_parameter[0] = (UParameter - UCacheParameter) / USpanLenght ;
2195     new_parameter[1] = (VParameter - VCacheParameter) / VSpanLenght ; 
2196     dimension = 3 * (VDegree + 1) ;
2197   }
2198   NCollection_LocalArray<Standard_Real> locpoles(dimension);
2199   
2200   PLib::NoDerivativeEvalPolynomial(new_parameter[0],
2201                        max_degree,
2202                        dimension,
2203                        max_degree*dimension,
2204                        PArray[0],
2205                        locpoles[0]) ;
2206   
2207   PLib::NoDerivativeEvalPolynomial(new_parameter[1],
2208                        min_degree,
2209                        3,
2210                        (min_degree << 1) + min_degree,
2211                        locpoles[0],
2212                        myPoint[0]) ;
2213   if (WeightsArray != NULL) {
2214     dimension = min_degree + 1 ;
2215     const TColStd_Array2OfReal& refWeights = *WeightsArray;
2216     Standard_Real *
2217       WArray = (Standard_Real *) 
2218         &refWeights(refWeights.LowerCol(),refWeights.LowerRow()) ;
2219     PLib::NoDerivativeEvalPolynomial(new_parameter[0],
2220                          max_degree,
2221                          dimension,
2222                          max_degree*dimension,
2223                          WArray[0],
2224                          locpoles[0]) ;
2225     
2226     PLib::NoDerivativeEvalPolynomial(new_parameter[1],
2227                          min_degree,
2228                          1,
2229                          min_degree,
2230                          locpoles[0],
2231                          inverse) ;
2232     inverse = 1.0e0/ inverse ;
2233     
2234     myPoint[0] *= inverse ;
2235     myPoint[1] *= inverse ;
2236     myPoint[2] *= inverse ;
2237   }
2238 }
2239 
2240 //=======================================================================
2241 //function : CacheD1
2242 //purpose  : Evaluates the polynomial cache of the Bspline Curve
2243 //           
2244 //=======================================================================
2245 
2246 void  BSplSLib::CacheD1(const Standard_Real                  UParameter,
2247                         const Standard_Real                  VParameter,
2248                         const  Standard_Integer              UDegree,
2249                         const  Standard_Integer              VDegree,
2250                         const  Standard_Real                 UCacheParameter,
2251                         const  Standard_Real                 VCacheParameter,
2252                         const  Standard_Real                 USpanLenght,
2253                         const  Standard_Real                 VSpanLenght,
2254                         const  TColgp_Array2OfPnt&           PolesArray,
2255                         const  TColStd_Array2OfReal*         WeightsArray,
2256                         gp_Pnt&                              aPoint,
2257                         gp_Vec&                              aVecU,
2258                         gp_Vec&                              aVecV)
2259 {
2260   //
2261   // the CacheParameter is where the cache polynomial was evaluated in homogeneous
2262   // form
2263   // the SpanLenght     is the normalizing factor so that the polynomial is between
2264   // 0 and 1 
2265   Standard_Integer 
2266 //    ii,
2267 //  jj,
2268 //  kk,
2269   dimension,
2270   min_degree,
2271   max_degree  ;
2272   
2273   Standard_Real
2274     inverse_min,
2275   inverse_max, 
2276   new_parameter[2]  ;
2277   
2278   Standard_Real * 
2279     PArray = (Standard_Real *) 
2280       &(PolesArray(PolesArray.LowerCol(), PolesArray.LowerRow())) ;
2281   Standard_Real local_poles_array[2][2][3],
2282   local_poles_and_weights_array[2][2][4],
2283   local_weights_array[2][2]  ;
2284   Standard_Real * my_vec_min,
2285   * my_vec_max,
2286   * my_point ;
2287   my_point = (Standard_Real *) &aPoint  ;
2288   //
2289   // initialize in case of rational evaluation
2290   // because RationalDerivative will use all
2291   // the coefficients
2292   //
2293   //
2294   if (WeightsArray != NULL) {
2295 
2296     local_poles_array            [0][0][0] = 0.0e0 ;
2297     local_poles_array            [0][0][1] = 0.0e0 ;
2298     local_poles_array            [0][0][2] = 0.0e0 ;
2299     local_weights_array          [0][0]    = 0.0e0 ;
2300     local_poles_and_weights_array[0][0][0] = 0.0e0 ;
2301     local_poles_and_weights_array[0][0][1] = 0.0e0 ;
2302     local_poles_and_weights_array[0][0][2] = 0.0e0 ;
2303     local_poles_and_weights_array[0][0][3] = 0.0e0 ;
2304     
2305     local_poles_array            [0][1][0] = 0.0e0 ;
2306     local_poles_array            [0][1][1] = 0.0e0 ;
2307     local_poles_array            [0][1][2] = 0.0e0 ;
2308     local_weights_array          [0][1]    = 0.0e0 ;
2309     local_poles_and_weights_array[0][1][0] = 0.0e0 ;
2310     local_poles_and_weights_array[0][1][1] = 0.0e0 ;
2311     local_poles_and_weights_array[0][1][2] = 0.0e0 ;
2312     local_poles_and_weights_array[0][1][3] = 0.0e0 ;
2313 
2314     local_poles_array            [1][0][0] = 0.0e0 ;
2315     local_poles_array            [1][0][1] = 0.0e0 ;
2316     local_poles_array            [1][0][2] = 0.0e0 ;
2317     local_weights_array          [1][0]    = 0.0e0 ;
2318     local_poles_and_weights_array[1][0][0] = 0.0e0 ;
2319     local_poles_and_weights_array[1][0][1] = 0.0e0 ;
2320     local_poles_and_weights_array[1][0][2] = 0.0e0 ;
2321     local_poles_and_weights_array[1][0][3] = 0.0e0 ;
2322     
2323     local_poles_array            [1][1][0] = 0.0e0 ;
2324     local_poles_array            [1][1][1] = 0.0e0 ;
2325     local_poles_array            [1][1][2] = 0.0e0 ;
2326     local_weights_array          [1][1]    = 0.0e0 ;
2327     local_poles_and_weights_array[1][1][0] = 0.0e0 ;
2328     local_poles_and_weights_array[1][1][1] = 0.0e0 ;
2329     local_poles_and_weights_array[1][1][2] = 0.0e0 ;
2330     local_poles_and_weights_array[1][1][3] = 0.0e0 ;
2331   }
2332 
2333   if (UDegree <= VDegree) {
2334     min_degree = UDegree ;
2335     max_degree = VDegree ;
2336     inverse_min = 1.0e0/USpanLenght ;
2337     inverse_max = 1.0e0/VSpanLenght ;
2338     new_parameter[0] = (VParameter - VCacheParameter) * inverse_max ; 
2339     new_parameter[1] = (UParameter - UCacheParameter) * inverse_min ;
2340     
2341     dimension = 3 * (UDegree + 1) ;
2342     my_vec_min = (Standard_Real *) &aVecU ;
2343     my_vec_max = (Standard_Real *) &aVecV ;
2344   }
2345   else {
2346     min_degree = VDegree ;
2347     max_degree = UDegree ;
2348     inverse_min = 1.0e0/VSpanLenght ;
2349     inverse_max = 1.0e0/USpanLenght ;
2350     new_parameter[0] = (UParameter - UCacheParameter) * inverse_max ;
2351     new_parameter[1] = (VParameter - VCacheParameter) * inverse_min ; 
2352     dimension = 3 * (VDegree + 1) ;
2353     my_vec_min = (Standard_Real *) &aVecV ;
2354     my_vec_max = (Standard_Real *) &aVecU ;
2355   }
2356 
2357   NCollection_LocalArray<Standard_Real> locpoles (2 * dimension);
2358   
2359   PLib::EvalPolynomial(new_parameter[0],
2360                        1,
2361                        max_degree,
2362                        dimension,
2363                        PArray[0],
2364                        locpoles[0]) ;
2365   
2366   PLib::EvalPolynomial(new_parameter[1],
2367                        1,
2368                        min_degree,
2369                        3,
2370                        locpoles[0],
2371                        local_poles_array[0][0][0]) ;
2372   PLib::NoDerivativeEvalPolynomial(new_parameter[1],
2373                        min_degree,
2374                        3,
2375                        (min_degree << 1) + min_degree,
2376                        locpoles[dimension],
2377                        local_poles_array[1][0][0]) ;
2378   
2379   if (WeightsArray != NULL) {
2380     dimension = min_degree + 1 ;
2381     const TColStd_Array2OfReal& refWeights = *WeightsArray;
2382     Standard_Real *
2383       WArray = (Standard_Real *) 
2384         &refWeights(refWeights.LowerCol(),refWeights.LowerRow()) ;
2385     PLib::EvalPolynomial(new_parameter[0],
2386                          1,
2387                          max_degree,
2388                          dimension,
2389                          WArray[0],
2390                          locpoles[0]) ;
2391     
2392     PLib::EvalPolynomial(new_parameter[1],
2393                          1,
2394                          min_degree,
2395                          1,
2396                          locpoles[0],
2397                          local_weights_array[0][0]) ;
2398     PLib::NoDerivativeEvalPolynomial(new_parameter[1],
2399                          min_degree,
2400                          1,
2401                          min_degree,
2402                          locpoles[dimension],
2403                          local_weights_array[1][0]) ;
2404     
2405     local_poles_and_weights_array[0][0][0] = local_poles_array  [0][0][0] ;
2406     local_poles_and_weights_array[0][0][1] = local_poles_array  [0][0][1] ;
2407     local_poles_and_weights_array[0][0][2] = local_poles_array  [0][0][2] ;
2408     local_poles_and_weights_array[0][0][3] = local_weights_array[0][0]    ;
2409     
2410     local_poles_and_weights_array[0][1][0] = local_poles_array  [0][1][0] ;
2411     local_poles_and_weights_array[0][1][1] = local_poles_array  [0][1][1] ;
2412     local_poles_and_weights_array[0][1][2] = local_poles_array  [0][1][2] ;
2413     local_poles_and_weights_array[0][1][3] = local_weights_array[0][1]    ;
2414     
2415     local_poles_and_weights_array[1][0][0] = local_poles_array  [1][0][0] ;
2416     local_poles_and_weights_array[1][0][1] = local_poles_array  [1][0][1] ;
2417     local_poles_and_weights_array[1][0][2] = local_poles_array  [1][0][2] ;
2418     local_poles_and_weights_array[1][0][3] = local_weights_array[1][0]    ;
2419     
2420     local_poles_and_weights_array[1][1][0] = local_poles_array  [1][1][0] ;
2421     local_poles_and_weights_array[1][1][1] = local_poles_array  [1][1][1] ;
2422     local_poles_and_weights_array[1][1][2] = local_poles_array  [1][1][2] ;
2423     local_poles_and_weights_array[1][1][3] = local_weights_array[1][1]    ;
2424 
2425     BSplSLib::RationalDerivative(1,
2426                                  1,
2427                                  1,
2428                                  1,
2429                                  local_poles_and_weights_array[0][0][0],
2430                                  local_poles_array[0][0][0]) ;
2431   }
2432   
2433   my_point  [0] = local_poles_array              [0][0][0] ;
2434   my_vec_min[0] = inverse_min * local_poles_array[0][1][0] ;
2435   my_vec_max[0] = inverse_max * local_poles_array[1][0][0] ;
2436   
2437   my_point  [1] = local_poles_array              [0][0][1] ;
2438   my_vec_min[1] = inverse_min * local_poles_array[0][1][1] ;
2439   my_vec_max[1] = inverse_max * local_poles_array[1][0][1] ;
2440   
2441   my_point  [2] = local_poles_array              [0][0][2] ;
2442   my_vec_min[2] = inverse_min * local_poles_array[0][1][2] ;
2443   my_vec_max[2] = inverse_max * local_poles_array[1][0][2] ;
2444 }
2445 
2446 //=======================================================================
2447 //function : CacheD2
2448 //purpose  : Evaluates the polynomial cache of the Bspline Curve
2449 //           
2450 //=======================================================================
2451 
2452 void  BSplSLib::CacheD2(const Standard_Real                  UParameter,
2453                         const Standard_Real                  VParameter,
2454                         const  Standard_Integer              UDegree,
2455                         const  Standard_Integer              VDegree,
2456                         const  Standard_Real                 UCacheParameter,
2457                         const  Standard_Real                 VCacheParameter,
2458                         const  Standard_Real                 USpanLenght,
2459                         const  Standard_Real                 VSpanLenght,
2460                         const  TColgp_Array2OfPnt&           PolesArray,
2461                         const  TColStd_Array2OfReal*         WeightsArray,
2462                         gp_Pnt&                              aPoint,
2463                         gp_Vec&                              aVecU,
2464                         gp_Vec&                              aVecV,
2465                         gp_Vec&                              aVecUU,
2466                         gp_Vec&                              aVecUV,
2467                         gp_Vec&                              aVecVV)
2468 {
2469   //
2470   // the CacheParameter is where the cache polynomial was evaluated in homogeneous
2471   // form
2472   // the SpanLenght     is the normalizing factor so that the polynomial is between
2473   // 0 and 1 
2474   Standard_Integer 
2475     ii,
2476 //  jj,
2477   kk,
2478   index,
2479   dimension,
2480   min_degree,
2481   max_degree  ;
2482   
2483   Standard_Real
2484     inverse_min,
2485   inverse_max, 
2486   new_parameter[2]  ;
2487 
2488   Standard_Real * 
2489     PArray = (Standard_Real *) 
2490       &(PolesArray(PolesArray.LowerCol(), PolesArray.LowerRow())) ;
2491   Standard_Real local_poles_array[3][3][3],
2492   local_poles_and_weights_array[3][3][4],
2493   local_weights_array[3][3]  ;
2494   Standard_Real * my_vec_min,
2495   * my_vec_max,
2496   * my_vec_min_min,
2497   * my_vec_max_max,
2498   * my_vec_min_max,
2499   * my_point ;
2500   my_point = (Standard_Real *) &aPoint  ;
2501   
2502   //
2503   // initialize in case the min and max degree are less than 2
2504   //
2505   local_poles_array[0][0][0] = 0.0e0 ;
2506   local_poles_array[0][0][1] = 0.0e0 ;
2507   local_poles_array[0][0][2] = 0.0e0 ;
2508   local_poles_array[0][1][0] = 0.0e0 ;
2509   local_poles_array[0][1][1] = 0.0e0 ;
2510   local_poles_array[0][1][2] = 0.0e0 ;
2511   local_poles_array[0][2][0] = 0.0e0 ;
2512   local_poles_array[0][2][1] = 0.0e0 ;
2513   local_poles_array[0][2][2] = 0.0e0 ;
2514   
2515   local_poles_array[1][0][0] = 0.0e0 ;
2516   local_poles_array[1][0][1] = 0.0e0 ;
2517   local_poles_array[1][0][2] = 0.0e0 ;
2518   local_poles_array[1][1][0] = 0.0e0 ;
2519   local_poles_array[1][1][1] = 0.0e0 ;
2520   local_poles_array[1][1][2] = 0.0e0 ;
2521   local_poles_array[1][2][0] = 0.0e0 ;
2522   local_poles_array[1][2][1] = 0.0e0 ;
2523   local_poles_array[1][2][2] = 0.0e0 ;
2524   
2525   local_poles_array[2][0][0] = 0.0e0 ;
2526   local_poles_array[2][0][1] = 0.0e0 ;
2527   local_poles_array[2][0][2] = 0.0e0 ;
2528   local_poles_array[2][1][0] = 0.0e0 ;
2529   local_poles_array[2][1][1] = 0.0e0 ;
2530   local_poles_array[2][1][2] = 0.0e0 ;
2531   local_poles_array[2][2][0] = 0.0e0 ;
2532   local_poles_array[2][2][1] = 0.0e0 ;
2533   local_poles_array[2][2][2] = 0.0e0 ;
2534   //
2535   // initialize in case of rational evaluation
2536   // because RationalDerivative will use all
2537   // the coefficients
2538   //
2539   //
2540   if (WeightsArray != NULL) {
2541     
2542     local_poles_and_weights_array[0][0][0] = 0.0e0 ;
2543     local_poles_and_weights_array[0][0][1] = 0.0e0 ;
2544     local_poles_and_weights_array[0][0][2] = 0.0e0 ;
2545     local_poles_and_weights_array[0][1][0] = 0.0e0 ;
2546     local_poles_and_weights_array[0][1][1] = 0.0e0 ;
2547     local_poles_and_weights_array[0][1][2] = 0.0e0 ;
2548     local_poles_and_weights_array[0][2][0] = 0.0e0 ;
2549     local_poles_and_weights_array[0][2][1] = 0.0e0 ;
2550     local_poles_and_weights_array[0][2][2] = 0.0e0 ;
2551     
2552     local_poles_and_weights_array[1][0][0] = 0.0e0 ;
2553     local_poles_and_weights_array[1][0][1] = 0.0e0 ;
2554     local_poles_and_weights_array[1][0][2] = 0.0e0 ;
2555     local_poles_and_weights_array[1][1][0] = 0.0e0 ;
2556     local_poles_and_weights_array[1][1][1] = 0.0e0 ;
2557     local_poles_and_weights_array[1][1][2] = 0.0e0 ;
2558     local_poles_and_weights_array[1][2][0] = 0.0e0 ;
2559     local_poles_and_weights_array[1][2][1] = 0.0e0 ;
2560     local_poles_and_weights_array[1][2][2] = 0.0e0 ;
2561     
2562     local_poles_and_weights_array[2][0][0] = 0.0e0 ;
2563     local_poles_and_weights_array[2][0][1] = 0.0e0 ;
2564     local_poles_and_weights_array[2][0][2] = 0.0e0 ;
2565     local_poles_and_weights_array[2][1][0] = 0.0e0 ;
2566     local_poles_and_weights_array[2][1][1] = 0.0e0 ;
2567     local_poles_and_weights_array[2][1][2] = 0.0e0 ;
2568     local_poles_and_weights_array[2][2][0] = 0.0e0 ;
2569     local_poles_and_weights_array[2][2][1] = 0.0e0 ;
2570     local_poles_and_weights_array[2][2][2] = 0.0e0 ;
2571     
2572     local_poles_and_weights_array[0][0][3] =
2573       local_weights_array[0][0] = 0.0e0 ;
2574     local_poles_and_weights_array[0][1][3] =
2575       local_weights_array[0][1] = 0.0e0 ;
2576     local_poles_and_weights_array[0][2][3] =
2577       local_weights_array[0][2] = 0.0e0 ;
2578     local_poles_and_weights_array[1][0][3] =
2579       local_weights_array[1][0] = 0.0e0 ;
2580     local_poles_and_weights_array[1][1][3] =
2581       local_weights_array[1][1] = 0.0e0 ;
2582     local_poles_and_weights_array[1][2][3] =
2583       local_weights_array[1][2] = 0.0e0 ;
2584     local_poles_and_weights_array[2][0][3] =
2585       local_weights_array[2][0] = 0.0e0 ;
2586     local_poles_and_weights_array[2][1][3] =
2587       local_weights_array[2][1] = 0.0e0 ;
2588     local_poles_and_weights_array[2][2][3] =
2589       local_weights_array[2][2] = 0.0e0 ;
2590   }
2591 
2592   if (UDegree <= VDegree) {
2593     min_degree = UDegree ;
2594     max_degree = VDegree ;
2595     inverse_min = 1.0e0/USpanLenght ;
2596     inverse_max = 1.0e0/VSpanLenght ;
2597     new_parameter[0] = (VParameter - VCacheParameter) * inverse_max ; 
2598     new_parameter[1] = (UParameter - UCacheParameter) * inverse_min ;
2599     
2600     dimension = 3 * (UDegree + 1) ;
2601     my_vec_min     = (Standard_Real *) &aVecU ;
2602     my_vec_max     = (Standard_Real *) &aVecV ;
2603     my_vec_min_min = (Standard_Real *) &aVecUU ;
2604     my_vec_min_max = (Standard_Real *) &aVecUV ;
2605     my_vec_max_max = (Standard_Real *) &aVecVV ;
2606   }
2607   else {
2608     min_degree = VDegree ;
2609     max_degree = UDegree ;
2610     inverse_min = 1.0e0/VSpanLenght ;
2611     inverse_max = 1.0e0/USpanLenght ;
2612     new_parameter[0] = (UParameter - UCacheParameter) * inverse_max ;
2613     new_parameter[1] = (VParameter - VCacheParameter) * inverse_min ; 
2614     dimension = 3 * (VDegree + 1) ;
2615     my_vec_min     = (Standard_Real *) &aVecV ;
2616     my_vec_max     = (Standard_Real *) &aVecU ;
2617     my_vec_min_min = (Standard_Real *) &aVecVV ;
2618     my_vec_min_max = (Standard_Real *) &aVecUV ;
2619     my_vec_max_max = (Standard_Real *) &aVecUU ;
2620   }
2621 
2622   NCollection_LocalArray<Standard_Real> locpoles (3 * dimension);
2623   
2624   //
2625   // initialize in case min or max degree are less than 2
2626   //
2627   Standard_Integer MinIndMax = 2;
2628   if ( max_degree < 2) MinIndMax = max_degree;
2629   Standard_Integer MinIndMin = 2;
2630   if ( min_degree < 2) MinIndMin = min_degree;
2631   
2632   index =  MinIndMax * dimension ;
2633 
2634   for (ii = MinIndMax ; ii <  3 ; ii++) {
2635     
2636     for (kk = 0 ; kk < dimension ; kk++) {
2637       locpoles[index] = 0.0e0 ;
2638       index += 1 ;
2639     }
2640   }
2641   
2642   PLib::EvalPolynomial(new_parameter[0],
2643                        MinIndMax,
2644                        max_degree,
2645                        dimension,
2646                        PArray[0],
2647                        locpoles[0]) ;
2648   
2649   PLib::EvalPolynomial(new_parameter[1],
2650                        MinIndMin,
2651                        min_degree,
2652                        3,
2653                        locpoles[0],
2654                        local_poles_array[0][0][0]) ;
2655   PLib::EvalPolynomial(new_parameter[1],
2656                        1,
2657                        min_degree,
2658                        3,
2659                        locpoles[dimension],
2660                        local_poles_array[1][0][0]) ;
2661   
2662   PLib::NoDerivativeEvalPolynomial(new_parameter[1],
2663                        min_degree,
2664                        3,
2665                        (min_degree << 1) + min_degree,
2666                        locpoles[dimension + dimension],
2667                        local_poles_array[2][0][0]) ;
2668   
2669   if (WeightsArray != NULL) {
2670     dimension = min_degree + 1 ;
2671     const TColStd_Array2OfReal& refWeights = *WeightsArray;
2672     Standard_Real *
2673       WArray = (Standard_Real *) 
2674         &refWeights(refWeights.LowerCol(),refWeights.LowerRow()) ;
2675     PLib::EvalPolynomial(new_parameter[0],
2676                          MinIndMax,
2677                          max_degree,
2678                          dimension,
2679                          WArray[0],
2680                          locpoles[0]) ;
2681     
2682     PLib::EvalPolynomial(new_parameter[1],
2683                          MinIndMin,
2684                          min_degree,
2685                          1,
2686                          locpoles[0],
2687                          local_weights_array[0][0]) ;
2688     PLib::EvalPolynomial(new_parameter[1],
2689                          1,
2690                          min_degree,
2691                          1,
2692                          locpoles[dimension],
2693                          local_weights_array[1][0]) ;
2694     PLib::NoDerivativeEvalPolynomial(new_parameter[1],
2695                          min_degree,
2696                          1,
2697                          min_degree,
2698                          locpoles[dimension + dimension],
2699                          local_weights_array[2][0]) ;
2700     
2701     
2702     local_poles_and_weights_array[0][0][0] = local_poles_array[0][0][0];
2703     local_poles_and_weights_array[0][0][1] = local_poles_array[0][0][1];
2704     local_poles_and_weights_array[0][0][2] = local_poles_array[0][0][2];
2705     local_poles_and_weights_array[0][1][0] = local_poles_array[0][1][0];
2706     local_poles_and_weights_array[0][1][1] = local_poles_array[0][1][1];
2707     local_poles_and_weights_array[0][1][2] = local_poles_array[0][1][2];
2708     local_poles_and_weights_array[0][2][0] = local_poles_array[0][2][0];
2709     local_poles_and_weights_array[0][2][1] = local_poles_array[0][2][1];
2710     local_poles_and_weights_array[0][2][2] = local_poles_array[0][2][2];
2711     
2712     local_poles_and_weights_array[1][0][0] = local_poles_array[1][0][0];
2713     local_poles_and_weights_array[1][0][1] = local_poles_array[1][0][1];
2714     local_poles_and_weights_array[1][0][2] = local_poles_array[1][0][2];
2715     local_poles_and_weights_array[1][1][0] = local_poles_array[1][1][0];
2716     local_poles_and_weights_array[1][1][1] = local_poles_array[1][1][1];
2717     local_poles_and_weights_array[1][1][2] = local_poles_array[1][1][2];
2718     local_poles_and_weights_array[1][2][0] = local_poles_array[1][2][0];
2719     local_poles_and_weights_array[1][2][1] = local_poles_array[1][2][1];
2720     local_poles_and_weights_array[1][2][2] = local_poles_array[1][2][2];
2721     
2722     local_poles_and_weights_array[2][0][0] = local_poles_array[2][0][0];
2723     local_poles_and_weights_array[2][0][1] = local_poles_array[2][0][1];
2724     local_poles_and_weights_array[2][0][2] = local_poles_array[2][0][2];
2725     local_poles_and_weights_array[2][1][0] = local_poles_array[2][1][0];
2726     local_poles_and_weights_array[2][1][1] = local_poles_array[2][1][1];
2727     local_poles_and_weights_array[2][1][2] = local_poles_array[2][1][2];
2728     local_poles_and_weights_array[2][2][0] = local_poles_array[2][2][0];
2729     local_poles_and_weights_array[2][2][1] = local_poles_array[2][2][1];
2730     local_poles_and_weights_array[2][2][2] = local_poles_array[2][2][2];
2731     
2732     
2733     local_poles_and_weights_array[0][0][3] = local_weights_array[0][0];
2734     local_poles_and_weights_array[0][1][3] = local_weights_array[0][1];
2735     local_poles_and_weights_array[0][2][3] = local_weights_array[0][2];
2736     local_poles_and_weights_array[1][0][3] = local_weights_array[1][0];
2737     local_poles_and_weights_array[1][1][3] = local_weights_array[1][1];
2738     local_poles_and_weights_array[1][2][3] = local_weights_array[1][2];
2739     local_poles_and_weights_array[2][0][3] = local_weights_array[2][0];
2740     local_poles_and_weights_array[2][1][3] = local_weights_array[2][1];
2741     local_poles_and_weights_array[2][2][3] = local_weights_array[2][2];
2742     
2743     BSplSLib::RationalDerivative(2,
2744                                  2,
2745                                  2,
2746                                  2,
2747                                  local_poles_and_weights_array[0][0][0],
2748                                  local_poles_array[0][0][0]) ;
2749   }
2750   
2751 
2752   Standard_Real minmin = inverse_min * inverse_min;
2753   Standard_Real minmax = inverse_min * inverse_max;
2754   Standard_Real maxmax = inverse_max * inverse_max;
2755   
2756   my_point      [0] = local_poles_array              [0][0][0] ;
2757   my_vec_min    [0] = inverse_min * local_poles_array[0][1][0] ;
2758   my_vec_max    [0] = inverse_max * local_poles_array[1][0][0] ;
2759   my_vec_min_min[0] = minmin * local_poles_array     [0][2][0] ;
2760   my_vec_min_max[0] = minmax * local_poles_array     [1][1][0] ;
2761   my_vec_max_max[0] = maxmax * local_poles_array     [2][0][0] ;
2762 
2763   my_point      [1] = local_poles_array              [0][0][1] ;
2764   my_vec_min    [1] = inverse_min * local_poles_array[0][1][1] ;
2765   my_vec_max    [1] = inverse_max * local_poles_array[1][0][1] ;
2766   my_vec_min_min[1] = minmin * local_poles_array     [0][2][1] ;
2767   my_vec_min_max[1] = minmax * local_poles_array     [1][1][1] ;
2768   my_vec_max_max[1] = maxmax * local_poles_array     [2][0][1] ;
2769 
2770   my_point      [2] = local_poles_array              [0][0][2] ;
2771   my_vec_min    [2] = inverse_min * local_poles_array[0][1][2] ;
2772   my_vec_max    [2] = inverse_max * local_poles_array[1][0][2] ;
2773   my_vec_min_min[2] = minmin * local_poles_array     [0][2][2] ;
2774   my_vec_min_max[2] = minmax * local_poles_array     [1][1][2] ;
2775   my_vec_max_max[2] = maxmax * local_poles_array     [2][0][2] ;
2776 }
2777 
2778 //=======================================================================
2779 //function : MovePoint
2780 //purpose  : Find the new poles which allows  an old point (with a
2781 //           given  u and v as parameters) to reach a new position
2782 //=======================================================================
2783 
2784 void BSplSLib::MovePoint (const Standard_Real            U, 
2785                           const Standard_Real            V,
2786                           const gp_Vec&                  Displ,
2787                           const Standard_Integer         UIndex1,
2788                           const Standard_Integer         UIndex2,
2789                           const Standard_Integer         VIndex1,
2790                           const Standard_Integer         VIndex2,
2791                           const Standard_Integer         UDegree,
2792                           const Standard_Integer         VDegree,
2793                           const Standard_Boolean         Rational,
2794                           const TColgp_Array2OfPnt&      Poles,  
2795                           const TColStd_Array2OfReal&    Weights,
2796                           const TColStd_Array1OfReal&    UFlatKnots,
2797                           const TColStd_Array1OfReal&    VFlatKnots,
2798                           Standard_Integer&              UFirstIndex,
2799                           Standard_Integer&              ULastIndex,
2800                           Standard_Integer&              VFirstIndex,
2801                           Standard_Integer&              VLastIndex,
2802                           TColgp_Array2OfPnt&            NewPoles)
2803 {
2804   // calculate the UBSplineBasis in the parameter U
2805   Standard_Integer UFirstNonZeroBsplineIndex;
2806   math_Matrix UBSplineBasis(1, 1,
2807                             1, UDegree+1);
2808   Standard_Integer ErrorCod1 =  BSplCLib::EvalBsplineBasis(0,
2809                                                            UDegree+1,
2810                                                            UFlatKnots,
2811                                                            U,
2812                                                            UFirstNonZeroBsplineIndex,
2813                                                            UBSplineBasis);  
2814   // calculate the VBSplineBasis in the parameter V
2815   Standard_Integer VFirstNonZeroBsplineIndex;
2816   math_Matrix VBSplineBasis(1, 1,
2817                             1, VDegree+1);
2818   Standard_Integer ErrorCod2 =  BSplCLib::EvalBsplineBasis(0,
2819                                                            VDegree+1,
2820                                                            VFlatKnots,
2821                                                            V,
2822                                                            VFirstNonZeroBsplineIndex,
2823                                                            VBSplineBasis);  
2824   if (ErrorCod1 || ErrorCod2) {
2825     UFirstIndex = 0;
2826     ULastIndex = 0;
2827     VFirstIndex = 0;
2828     VLastIndex = 0;
2829     return;
2830   }
2831   
2832   // find the span which is predominant for parameter U
2833   UFirstIndex = UFirstNonZeroBsplineIndex;
2834   ULastIndex = UFirstNonZeroBsplineIndex + UDegree ;
2835   if (UFirstIndex < UIndex1) UFirstIndex = UIndex1;
2836   if (ULastIndex > UIndex2) ULastIndex = UIndex2;
2837 
2838   Standard_Real maxValue = 0.0;
2839   Standard_Integer i, ukk1=0, ukk2;
2840 
2841   for (i = UFirstIndex-UFirstNonZeroBsplineIndex+1; i <= ULastIndex-UFirstNonZeroBsplineIndex+1; i++) {
2842     if (UBSplineBasis(1,i) > maxValue) {
2843       ukk1 = i + UFirstNonZeroBsplineIndex - 1;
2844       maxValue = UBSplineBasis(1,i);
2845     }
2846   }
2847 
2848   // find a ukk2 if symetriy
2849   ukk2 = ukk1;
2850   i = ukk1 - UFirstNonZeroBsplineIndex + 2;
2851   if ((ukk1+1) <= ULastIndex) {
2852     if (Abs(UBSplineBasis(1, ukk1-UFirstNonZeroBsplineIndex+2) - maxValue) < 1.e-10) {
2853       ukk2 = ukk1+1;
2854     }
2855   }
2856 
2857   // find the span which is predominant for parameter V
2858   VFirstIndex = VFirstNonZeroBsplineIndex;
2859   VLastIndex = VFirstNonZeroBsplineIndex + VDegree ;
2860 
2861   if (VFirstIndex < VIndex1) VFirstIndex = VIndex1;
2862   if (VLastIndex > VIndex2) VLastIndex = VIndex2;
2863 
2864   maxValue = 0.0;
2865   Standard_Integer j, vkk1=0, vkk2;
2866 
2867   for (j = VFirstIndex-VFirstNonZeroBsplineIndex+1; j <= VLastIndex-VFirstNonZeroBsplineIndex+1; j++) {
2868     if (VBSplineBasis(1,j) > maxValue) {
2869       vkk1 = j + VFirstNonZeroBsplineIndex - 1;
2870       maxValue = VBSplineBasis(1,j);
2871     }
2872   }
2873 
2874   // find a vkk2 if symetriy
2875   vkk2 = vkk1;
2876   j = vkk1 - VFirstNonZeroBsplineIndex + 2;
2877   if ((vkk1+1) <= VLastIndex) {
2878     if (Abs(VBSplineBasis(1, vkk1-VFirstNonZeroBsplineIndex+2) - maxValue) < 1.e-10) {
2879       vkk2 = vkk1+1;
2880     }
2881   }
2882 
2883   // compute the vector of displacement
2884   Standard_Real D1 = 0.0;
2885   Standard_Real D2 = 0.0;
2886   Standard_Real hN, Coef, DvalU, DvalV;
2887 
2888   Standard_Integer ii, jj;
2889 
2890   for (i = 1; i <= UDegree+1; i++) {
2891     ii = i + UFirstNonZeroBsplineIndex - 1;
2892     if (ii < ukk1) {
2893       DvalU = ukk1-ii;
2894     }
2895     else if (ii > ukk2) {
2896       DvalU = ii - ukk2;
2897     }
2898     else {
2899       DvalU = 0.0;
2900     }
2901 
2902     for (j = 1; j <= VDegree+1; j++) {
2903       jj = j + VFirstNonZeroBsplineIndex - 1;
2904       if (Rational) {
2905         hN = Weights(ii, jj)*UBSplineBasis(1, i)*VBSplineBasis(1,j);
2906         D2 += hN;
2907       }
2908       else {
2909         hN = UBSplineBasis(1, i)*VBSplineBasis(1,j);
2910       }
2911       if (ii >= UFirstIndex && ii <= ULastIndex && jj >= VFirstIndex && jj <= VLastIndex) {
2912         if (jj < vkk1) {
2913           DvalV = vkk1-jj;
2914         }
2915         else if (jj > vkk2) {
2916           DvalV = jj - vkk2;
2917         }
2918         else {
2919           DvalV = 0.0;
2920         }
2921         D1 += 1./(DvalU + DvalV + 1.) * hN;
2922       }
2923     }
2924   }
2925   
2926   if (Rational) {
2927     Coef = D2/D1;
2928   }
2929   else {
2930     Coef = 1./D1;
2931   }
2932 
2933   // compute the new poles
2934 
2935   for (i=Poles.LowerRow(); i<=Poles.UpperRow(); i++) {
2936     if (i < ukk1) {
2937       DvalU = ukk1-i;
2938     }
2939     else if (i > ukk2) {
2940       DvalU = i - ukk2;
2941     }
2942     else {
2943       DvalU = 0.0;
2944     }
2945 
2946     for (j=Poles.LowerCol(); j<=Poles.UpperCol(); j++) {
2947       if (i >= UFirstIndex && i <= ULastIndex && j >= VFirstIndex && j <= VLastIndex) {
2948         if (j < vkk1) {
2949           DvalV = vkk1-j;
2950         }
2951         else if (j > vkk2) {
2952           DvalV = j - vkk2;
2953         }
2954         else {
2955           DvalV = 0.0;
2956         }
2957         NewPoles(i,j) = Poles(i,j).Translated((Coef/(DvalU + DvalV + 1.))*Displ);
2958       }
2959       else {
2960         NewPoles(i,j) = Poles(i,j);
2961       }
2962     }
2963   }
2964 }
2965 
2966 //=======================================================================
2967 // function : Resolution
2968 // purpose  : this computes an estimate for the maximum of the 
2969 // partial derivatives both in U and in V
2970 //
2971 //
2972 // The calculation resembles at the calculation of curves with 
2973 // additional index for the control point. Let Si,j be the
2974 // control points for ls surface  and  Di,j  the weights.  
2975 // The checking of upper bounds for the partial derivatives 
2976 // will be omitted and Su is the next upper bound in the polynomial case :
2977 //
2978 //
2979 //
2980 //                        |  Si,j - Si-1,j  |
2981 //          d *   Max     |  -------------  |
2982 //                i = 2,n |     ti+d - ti   |
2983 //                i=1.m
2984 //
2985 //
2986 // and in the rational case :
2987 //
2988 //
2989 //
2990 //                         Di,j * (Si,j - Sk,j) - Di-1,j * (Si-1,j - Sk,j)
2991 //   Max   Max       d  *  -----------------------------------------------
2992 // k=1,n  i dans Rj                   ti+d  - ti
2993 // j=1,m
2994 //  ----------------------------------------------------------------------
2995 //
2996 //               Min    Di,j
2997 //              i=1,n
2998 //              j=1,m
2999 //
3000 //
3001 //
3002 // with Rj = {j-d, ....,  j+d+d+1}.
3003 //
3004 //
3005 //=======================================================================
3006 
3007 void BSplSLib::Resolution(const TColgp_Array2OfPnt&      Poles,
3008                           const TColStd_Array2OfReal*    Weights,
3009                           const TColStd_Array1OfReal&    UKnots,
3010                           const TColStd_Array1OfReal&    VKnots,
3011                           const TColStd_Array1OfInteger& UMults,
3012                           const TColStd_Array1OfInteger& VMults,
3013                           const Standard_Integer         UDegree,
3014                           const Standard_Integer         VDegree,
3015                           const Standard_Boolean         URational,
3016                           const Standard_Boolean         VRational,
3017                           const Standard_Boolean         UPeriodic,
3018                           const Standard_Boolean         VPeriodic,
3019                           const Standard_Real            Tolerance3D,
3020                           Standard_Real&                 UTolerance,
3021                           Standard_Real&                 VTolerance)
3022 {
3023   Standard_Real Wij,Wmj,Wji,Wjm;
3024   Standard_Real Xij,Xmj,Xji,Xjm,Xpq,Xqp;
3025   Standard_Real Yij,Ymj,Yji,Yjm,Ypq,Yqp;
3026   Standard_Real Zij,Zmj,Zji,Zjm,Zpq,Zqp;
3027   Standard_Real factor,value,min,min_weights=0,inverse,max_derivative[2];
3028 
3029   max_derivative[0] = max_derivative[1] = 0.0e0 ; 
3030   
3031   Standard_Integer PRowLength, PColLength;
3032   Standard_Integer ii,jj,pp,qq,ii_index,jj_index,pp_index,qq_index;
3033   Standard_Integer ii_minus,upper[2],lower[2],poles_length[2];
3034   Standard_Integer num_poles[2],num_flat_knots[2];
3035   
3036   num_flat_knots[0] = 
3037     BSplCLib::KnotSequenceLength(UMults,
3038                                  UDegree,
3039                                  UPeriodic) ;
3040   num_flat_knots[1] = 
3041     BSplCLib::KnotSequenceLength(VMults,
3042                                  VDegree,
3043                                  VPeriodic) ;
3044   TColStd_Array1OfReal  flat_knots_in_u(1,num_flat_knots[0]) ;
3045   TColStd_Array1OfReal  flat_knots_in_v(1,num_flat_knots[1]) ;
3046   BSplCLib::KnotSequence(UKnots,
3047                          UMults,
3048                          UDegree,
3049                          UPeriodic,
3050                          flat_knots_in_u) ;
3051   BSplCLib::KnotSequence(VKnots,
3052                          VMults,
3053                          VDegree,
3054                          VPeriodic,
3055                          flat_knots_in_v) ;
3056   PRowLength = Poles.RowLength();
3057   PColLength = Poles.ColLength();
3058   if (URational || VRational) {
3059     Standard_Integer Wsize = PRowLength * PColLength;
3060     const TColStd_Array2OfReal& refWights = *Weights;
3061     const Standard_Real * WG = &refWights(refWights.LowerRow(), refWights.LowerCol());
3062     min_weights = WG[0];
3063     
3064     for (ii = 1 ; ii < Wsize ; ii++) {
3065       min = WG[ii];
3066       if (min_weights > min) min_weights = min;
3067     }
3068   }
3069   Standard_Integer UD1 = UDegree + 1;
3070   Standard_Integer VD1 = VDegree + 1;
3071   num_poles[0] = num_flat_knots[0] - UD1;
3072   num_poles[1] = num_flat_knots[1] - VD1;
3073   poles_length[0] = PColLength;
3074   poles_length[1] = PRowLength;
3075   if(URational) {
3076     Standard_Integer UD2 = UDegree << 1;
3077     Standard_Integer VD2 = VDegree << 1;
3078 
3079     for (ii = 2 ; ii <= num_poles[0] ; ii++) {
3080       ii_index = (ii - 1) % poles_length[0] + 1 ;
3081       ii_minus = (ii - 2) % poles_length[0] + 1 ;
3082       inverse = flat_knots_in_u(ii + UDegree) - flat_knots_in_u(ii) ;
3083       inverse = 1.0e0 / inverse ;
3084       lower[0] = ii - UD1;
3085       if (lower[0] < 1) lower[0] = 1;
3086       upper[0] = ii + UD2 + 1;
3087       if (upper[0] > num_poles[0]) upper[0] = num_poles[0];
3088 
3089       for ( jj = 1 ; jj <= num_poles[1] ; jj++) {
3090         jj_index = (jj - 1) % poles_length[1] + 1 ;
3091         lower[1] = jj - VD1;
3092         if (lower[1] < 1) lower[1] = 1;
3093         upper[1] = jj + VD2 + 1;
3094         if (upper[1] > num_poles[1]) upper[1] = num_poles[1];
3095         const gp_Pnt& Pij = Poles  .Value(ii_index,jj_index);
3096         Wij               = Weights->Value(ii_index,jj_index);
3097         const gp_Pnt& Pmj = Poles  .Value(ii_minus,jj_index);
3098         Wmj               = Weights->Value(ii_minus,jj_index);
3099         Xij = Pij.X();
3100         Yij = Pij.Y();
3101         Zij = Pij.Z();
3102         Xmj = Pmj.X();
3103         Ymj = Pmj.Y();
3104         Zmj = Pmj.Z();
3105         
3106         for (pp = lower[0] ; pp <= upper[0] ; pp++) {
3107           pp_index = (pp - 1) % poles_length[0] + 1 ;
3108 
3109           for (qq = lower[1] ; qq <= upper[1] ; qq++) {
3110             value = 0.0e0 ;
3111             qq_index = (qq - 1) % poles_length[1] + 1 ;
3112             const gp_Pnt& Ppq = Poles.Value(pp_index,qq_index);
3113             Xpq = Ppq.X();
3114             Ypq = Ppq.Y();
3115             Zpq = Ppq.Z();
3116             factor  = (Xpq - Xij) * Wij;
3117             factor -= (Xpq - Xmj) * Wmj;
3118             if (factor < 0) factor = - factor;
3119             value  += factor ;
3120             factor  = (Ypq - Yij) * Wij;
3121             factor -= (Ypq - Ymj) * Wmj;
3122             if (factor < 0) factor = - factor;
3123             value  += factor ;
3124             factor  = (Zpq - Zij) * Wij;
3125             factor -= (Zpq - Zmj) * Wmj;
3126             if (factor < 0) factor = - factor;
3127             value  += factor ;
3128             value *= inverse ;
3129             if (max_derivative[0] < value) max_derivative[0] = value ;
3130           }
3131         }
3132       }
3133     }
3134     max_derivative[0] /= min_weights ;
3135   }
3136   else {
3137 
3138     for (ii = 2 ; ii <= num_poles[0] ; ii++) {
3139       ii_index = (ii - 1) % poles_length[0] + 1 ;
3140       ii_minus = (ii - 2) % poles_length[0] + 1 ;
3141       inverse = flat_knots_in_u(ii + UDegree) - flat_knots_in_u(ii) ;
3142       inverse = 1.0e0 / inverse ;
3143 
3144       for ( jj = 1 ; jj <= num_poles[1] ; jj++) {
3145         jj_index = (jj - 1) % poles_length[1] + 1 ;
3146         value = 0.0e0 ;
3147         const gp_Pnt& Pij = Poles.Value(ii_index,jj_index);
3148         const gp_Pnt& Pmj = Poles.Value(ii_minus,jj_index);
3149         factor = Pij.X() - Pmj.X();
3150         if (factor < 0) factor = - factor;
3151         value += factor;
3152         factor = Pij.Y() - Pmj.Y();
3153         if (factor < 0) factor = - factor;
3154         value += factor;
3155         factor = Pij.Z() - Pmj.Z();
3156         if (factor < 0) factor = - factor;
3157         value += factor;
3158         value *= inverse ;
3159         if (max_derivative[0] < value) max_derivative[0] = value ;
3160       }
3161     }
3162   }
3163   max_derivative[0] *= UDegree ;
3164   if(VRational) {
3165     Standard_Integer UD2 = UDegree << 1;
3166     Standard_Integer VD2 = VDegree << 1;
3167 
3168     for (ii = 2 ; ii <= num_poles[1] ; ii++) {
3169       ii_index = (ii - 1) % poles_length[1] + 1 ;
3170       ii_minus = (ii - 2) % poles_length[1] + 1 ;
3171       inverse = flat_knots_in_v(ii + VDegree) - flat_knots_in_v(ii) ;
3172       inverse = 1.0e0 / inverse ;
3173       lower[0] = ii - VD1;
3174       if (lower[0] < 1) lower[0] = 1;
3175       upper[0] = ii + VD2 + 1;
3176       if (upper[0] > num_poles[1]) upper[0] = num_poles[1];
3177 
3178       for ( jj = 1 ; jj <= num_poles[0] ; jj++) {
3179         jj_index = (jj - 1) % poles_length[0] + 1 ;
3180         lower[1] = jj - UD1;
3181         if (lower[1] < 1) lower[1] = 1;
3182         upper[1] = jj + UD2 + 1;
3183         if (upper[1] > num_poles[0]) upper[1] = num_poles[0];
3184         const gp_Pnt& Pji = Poles  .Value(jj_index,ii_index);
3185         Wji               = Weights->Value(jj_index,ii_index);
3186         const gp_Pnt& Pjm = Poles  .Value(jj_index,ii_minus);
3187         Wjm               = Weights->Value(jj_index,ii_minus);
3188         Xji = Pji.X();
3189         Yji = Pji.Y();
3190         Zji = Pji.Z();
3191         Xjm = Pjm.X();
3192         Yjm = Pjm.Y();
3193         Zjm = Pjm.Z();
3194         
3195         for (pp = lower[1] ; pp <= upper[1] ; pp++) {
3196           pp_index = (pp - 1) % poles_length[1] + 1 ;
3197 
3198           for (qq = lower[0] ; qq <= upper[0] ; qq++) {
3199             value = 0.0e0 ;
3200             qq_index = (qq - 1) % poles_length[0] + 1 ;
3201             const gp_Pnt& Pqp = Poles.Value(qq_index,pp_index);
3202             Xqp = Pqp.X();
3203             Yqp = Pqp.Y();
3204             Zqp = Pqp.Z();
3205             factor  = (Xqp - Xji) * Wji;
3206             factor -= (Xqp - Xjm) * Wjm;
3207             if (factor < 0) factor = - factor;
3208             value += factor ;
3209             factor  = (Yqp - Yji) * Wji;
3210             factor -= (Yqp - Yjm) * Wjm;
3211             if (factor < 0) factor = - factor;
3212             value += factor ;
3213             factor  = (Zqp - Zji) * Wji;
3214             factor -= (Zqp - Zjm) * Wjm;
3215             if (factor < 0) factor = - factor;
3216             value += factor ;
3217             value *= inverse ;
3218             if (max_derivative[1] < value) max_derivative[1] = value ;
3219           }
3220         }
3221       }
3222     }
3223     max_derivative[1] /= min_weights ;
3224   }
3225   else {
3226 
3227     for (ii = 2 ; ii <= num_poles[1] ; ii++) {
3228       ii_index = (ii - 1) % poles_length[1] + 1 ;
3229       ii_minus = (ii - 2) % poles_length[1] + 1 ;
3230       inverse = flat_knots_in_v(ii + VDegree) - flat_knots_in_v(ii) ;
3231       inverse = 1.0e0 / inverse ;
3232 
3233       for ( jj = 1 ; jj <= num_poles[0] ; jj++) {
3234         jj_index = (jj - 1) % poles_length[0] + 1 ;
3235         value = 0.0e0 ;
3236         const gp_Pnt& Pji = Poles.Value(jj_index,ii_index);
3237         const gp_Pnt& Pjm = Poles.Value(jj_index,ii_minus);
3238         factor = Pji.X() - Pjm.X() ;
3239         if (factor < 0) factor = - factor;
3240         value += factor;
3241         factor = Pji.Y() - Pjm.Y() ;
3242         if (factor < 0) factor = - factor;
3243         value += factor;
3244         factor = Pji.Z() - Pjm.Z() ;
3245         if (factor < 0) factor = - factor;
3246         value += factor;
3247         value *= inverse ;
3248         if (max_derivative[1] < value) max_derivative[1] = value ;
3249       }
3250     }
3251   }
3252   max_derivative[1] *= VDegree ;
3253   max_derivative[0] *= M_SQRT2 ;
3254   max_derivative[1] *= M_SQRT2 ;
3255   if(max_derivative[0] && max_derivative[1]) { 
3256     UTolerance = Tolerance3D / max_derivative[0] ;
3257     VTolerance = Tolerance3D / max_derivative[1] ;
3258   }
3259   else { 
3260     UTolerance=VTolerance=0.0;
3261 #ifdef OCCT_DEBUG
3262     std::cout<<"ElSLib.cxx : maxderivative = 0.0 "<<std::endl;
3263 #endif
3264   }
3265 }
3266 
3267 //=======================================================================
3268 //function : Interpolate
3269 //purpose  : 
3270 //=======================================================================
3271 
3272 void BSplSLib::Interpolate(const Standard_Integer UDegree,
3273                            const Standard_Integer VDegree, 
3274                            const TColStd_Array1OfReal& UFlatKnots,
3275                            const TColStd_Array1OfReal& VFlatKnots,
3276                            const TColStd_Array1OfReal& UParameters,
3277                            const TColStd_Array1OfReal& VParameters,
3278                            TColgp_Array2OfPnt& Poles,
3279                            TColStd_Array2OfReal& Weights, 
3280                            Standard_Integer& InversionProblem)
3281 {
3282   Standard_Integer ii, jj, ll, kk, dimension;
3283   Standard_Integer ULength = UParameters.Length();
3284   Standard_Integer VLength = VParameters.Length();
3285   Standard_Real * poles_array;
3286   
3287   // extraction of iso u
3288   dimension = 4*ULength;
3289   TColStd_Array2OfReal Points(1, VLength, 
3290                               1, dimension);
3291   
3292   Handle(TColStd_HArray1OfInteger) ContactOrder = 
3293     new (TColStd_HArray1OfInteger)(1, VLength);
3294   ContactOrder->Init(0);
3295   
3296   for (ii=1; ii <= VLength; ii++) {
3297 
3298     for (jj=1, ll=1; jj<= ULength; jj++, ll+=4) {
3299       Points(ii,ll)   = Poles(jj, ii).X();
3300       Points(ii,ll+1) = Poles(jj, ii).Y();
3301       Points(ii,ll+2) = Poles(jj, ii).Z();
3302       Points(ii,ll+3) = Weights(jj,ii)  ;
3303     }
3304   }
3305 
3306   // interpolation of iso u
3307   poles_array = (Standard_Real *) &Points.ChangeValue(1,1) ;
3308   BSplCLib::Interpolate(VDegree,
3309                         VFlatKnots,
3310                         VParameters,
3311                         ContactOrder->Array1(),
3312                         dimension,
3313                         poles_array[0],
3314                         InversionProblem) ;
3315   if (InversionProblem != 0) return;
3316 
3317   // extraction of iso v
3318 
3319   dimension = VLength*4;
3320   TColStd_Array2OfReal IsoPoles(1, ULength, 
3321                                 1, dimension);
3322   
3323   ContactOrder =  new (TColStd_HArray1OfInteger)(1, ULength);
3324   ContactOrder->Init(0);
3325   poles_array = (Standard_Real *) &IsoPoles.ChangeValue(1,1) ;
3326 
3327   for (ii=1, kk=1; ii <= ULength; ii++, kk+=4) {
3328 
3329     for (jj=1, ll=1; jj<= VLength; jj++, ll+=4) {
3330       IsoPoles (ii,ll)   = Points(jj, kk);
3331       IsoPoles (ii,ll+1) = Points(jj, kk+1);
3332       IsoPoles (ii,ll+2) = Points(jj, kk+2);
3333       IsoPoles (ii,ll+3) = Points(jj, kk+3);
3334     }
3335   }
3336   // interpolation of iso v
3337   BSplCLib::Interpolate(UDegree,
3338                         UFlatKnots,
3339                         UParameters,
3340                         ContactOrder->Array1(),
3341                         dimension,
3342                         poles_array[0],
3343                         InversionProblem);
3344 
3345   // return results
3346 
3347   for (ii=1; ii <= ULength; ii++) {
3348 
3349     for (jj=1, ll=1; jj<= VLength; jj++, ll+=4) {
3350       gp_Pnt Pnt(IsoPoles(ii,ll), IsoPoles(ii,ll+1), IsoPoles(ii,ll+2));
3351       Poles.SetValue(ii, jj, Pnt);
3352       Weights.SetValue(ii,jj,IsoPoles(ii,ll+3)) ;
3353     }
3354   }
3355 }
3356 
3357 //=======================================================================
3358 //function : Interpolate
3359 //purpose  : 
3360 //=======================================================================
3361 
3362 void BSplSLib::Interpolate(const Standard_Integer UDegree,
3363                            const Standard_Integer VDegree, 
3364                            const TColStd_Array1OfReal& UFlatKnots,
3365                            const TColStd_Array1OfReal& VFlatKnots,
3366                            const TColStd_Array1OfReal& UParameters,
3367                            const TColStd_Array1OfReal& VParameters,
3368                            TColgp_Array2OfPnt& Poles,
3369                            Standard_Integer& InversionProblem)
3370 {
3371   Standard_Integer ii, jj, ll, kk, dimension;
3372   Standard_Integer ULength = UParameters.Length();
3373   Standard_Integer VLength = VParameters.Length();
3374   Standard_Real * poles_array;
3375   
3376   // extraction of iso u
3377   dimension = 3*ULength;
3378   TColStd_Array2OfReal Points(1, VLength, 
3379                               1, dimension);
3380   
3381   Handle(TColStd_HArray1OfInteger) ContactOrder = 
3382     new (TColStd_HArray1OfInteger)(1, VLength);
3383   ContactOrder->Init(0);
3384   
3385   for (ii=1; ii <= VLength; ii++) {
3386     
3387     for (jj=1, ll=1; jj<= ULength; jj++, ll+=3) {
3388       Points(ii,ll)   = Poles(jj, ii).X();
3389       Points(ii,ll+1) = Poles(jj, ii).Y();
3390       Points(ii,ll+2) = Poles(jj, ii).Z();
3391     }
3392   }
3393   
3394   // interpolation of iso u
3395   poles_array = (Standard_Real *) &Points.ChangeValue(1,1) ;
3396   BSplCLib::Interpolate(VDegree,
3397                         VFlatKnots,
3398                         VParameters,
3399                         ContactOrder->Array1(),
3400                         dimension,
3401                         poles_array[0],
3402                         InversionProblem) ;
3403   if (InversionProblem != 0) return;
3404   
3405   // extraction of iso v
3406   
3407   dimension = VLength*3;
3408   TColStd_Array2OfReal IsoPoles(1, ULength, 
3409                                 1, dimension);
3410   
3411   ContactOrder =  new (TColStd_HArray1OfInteger)(1, ULength);
3412   ContactOrder->Init(0);
3413   poles_array = (Standard_Real *) &IsoPoles.ChangeValue(1,1) ;
3414   
3415   for (ii=1, kk=1; ii <= ULength; ii++, kk+=3) {
3416 
3417     for (jj=1, ll=1; jj<= VLength; jj++, ll+=3) {
3418       IsoPoles (ii,ll)   = Points(jj, kk);
3419       IsoPoles (ii,ll+1) = Points(jj, kk+1);
3420       IsoPoles (ii,ll+2) = Points(jj, kk+2);
3421     }
3422   }
3423   // interpolation of iso v
3424   BSplCLib::Interpolate(UDegree,
3425                         UFlatKnots,
3426                         UParameters,
3427                         ContactOrder->Array1(),
3428                         dimension,
3429                         poles_array[0],
3430                         InversionProblem);
3431   
3432   // return results
3433 
3434   for (ii=1; ii <= ULength; ii++) {
3435 
3436     for (jj=1, ll=1; jj<= VLength; jj++, ll+=3) {
3437       gp_Pnt Pnt(IsoPoles(ii,ll), IsoPoles(ii,ll+1), IsoPoles(ii,ll+2));
3438       Poles.SetValue(ii, jj, Pnt);
3439     }
3440   }
3441 }
3442 
3443 //=======================================================================
3444 //function : FunctionMultiply
3445 //purpose  : 
3446 //=======================================================================
3447 
3448 void BSplSLib::FunctionMultiply
3449 (const BSplSLib_EvaluatorFunction&           Function,
3450  const Standard_Integer                      UBSplineDegree,
3451  const Standard_Integer                      VBSplineDegree,
3452  const TColStd_Array1OfReal&                 UBSplineKnots,
3453  const TColStd_Array1OfReal&                 VBSplineKnots,
3454  const TColStd_Array1OfInteger *             UMults,
3455  const TColStd_Array1OfInteger *             VMults,
3456  const TColgp_Array2OfPnt&                   Poles,
3457  const TColStd_Array2OfReal*                 Weights,
3458  const TColStd_Array1OfReal&                 UFlatKnots,
3459  const TColStd_Array1OfReal&                 VFlatKnots,
3460  const Standard_Integer                      UNewDegree,
3461  const Standard_Integer                      VNewDegree,
3462  TColgp_Array2OfPnt&                         NewNumerator,
3463  TColStd_Array2OfReal&                       NewDenominator,
3464  Standard_Integer&                           theStatus)
3465 {
3466   Standard_Integer num_uparameters,
3467 //  ii,jj,kk,
3468   ii,jj,
3469   error_code,
3470   num_vparameters ;
3471   Standard_Real    result ;
3472   
3473   num_uparameters = UFlatKnots.Length() - UNewDegree - 1 ;
3474   num_vparameters = VFlatKnots.Length() - VNewDegree - 1 ;
3475   TColStd_Array1OfReal  UParameters(1,num_uparameters) ;
3476   TColStd_Array1OfReal  VParameters(1,num_vparameters) ;
3477   
3478   if ((NewNumerator.ColLength() == num_uparameters) &&
3479       (NewNumerator.RowLength() == num_vparameters) &&
3480       (NewDenominator.ColLength() == num_uparameters) &&
3481       (NewDenominator.RowLength() == num_vparameters)) {
3482     
3483     
3484     BSplCLib::BuildSchoenbergPoints(UNewDegree,
3485                                     UFlatKnots,
3486                                     UParameters) ;
3487     
3488     BSplCLib::BuildSchoenbergPoints(VNewDegree,
3489                                     VFlatKnots,
3490                                     VParameters) ;
3491     
3492     for (ii = 1 ; ii <= num_uparameters ; ii++) {
3493 
3494       for (jj = 1 ; jj <= num_vparameters ; jj++) {
3495         HomogeneousD0(UParameters(ii),
3496                       VParameters(jj),
3497                       0,
3498                       0,
3499                       Poles,
3500                       Weights,
3501                       UBSplineKnots,
3502                       VBSplineKnots,
3503                       UMults,
3504                       VMults,
3505                       UBSplineDegree,
3506                       VBSplineDegree,
3507                       Standard_True,
3508                       Standard_True,
3509                       Standard_False,
3510                       Standard_False,
3511                       NewDenominator(ii,jj),
3512                       NewNumerator(ii,jj)) ;
3513         
3514         Function.Evaluate (0,
3515                  UParameters(ii),
3516                  VParameters(jj),
3517                  result,
3518                  error_code) ;
3519         if (error_code) {
3520           throw Standard_ConstructionError();
3521         }
3522         gp_Pnt& P = NewNumerator(ii,jj);
3523         P.SetX(P.X() * result);
3524         P.SetY(P.Y() * result);
3525         P.SetZ(P.Z() * result);
3526         NewDenominator(ii,jj) *= result ;
3527       }
3528     }
3529     Interpolate(UNewDegree,
3530                 VNewDegree, 
3531                 UFlatKnots,
3532                 VFlatKnots,
3533                 UParameters,
3534                 VParameters,
3535                 NewNumerator,
3536                 NewDenominator, 
3537                 theStatus);
3538   }
3539   else {
3540     throw Standard_ConstructionError();
3541   }
3542 }
Open CASCADE Technology repositoryRSSAtom
