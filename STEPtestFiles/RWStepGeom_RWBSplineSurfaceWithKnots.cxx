gitprojects / occt.git / blob

commit
 ? search: 
 re
summary | shortlog | log | commit | commitdiff | tree
blame | history | raw | HEAD
0033261: Data Exchange, Step Import - Empty shape after reading process
[occt.git] / src / RWStepGeom / RWStepGeom_RWBSplineSurfaceWithKnots.cxx
   1 // Copyright (c) 1999-2014 OPEN CASCADE SAS
   2 //
   3 // This file is part of Open CASCADE Technology software library.
   4 //
   5 // This library is free software; you can redistribute it and/or modify it under
   6 // the terms of the GNU Lesser General Public License version 2.1 as published
   7 // by the Free Software Foundation, with special exception defined in the file
   8 // OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
   9 // distribution for complete text of the license and disclaimer of any warranty.
  10 //
  11 // Alternatively, this file may be used under the terms of Open CASCADE
  12 // commercial license or contractual agreement.
  13 
  14 
  15 #include <Interface_Check.hxx>
  16 #include <Interface_EntityIterator.hxx>
  17 #include <Interface_ShareTool.hxx>
  18 #include <RWStepGeom_RWBSplineSurfaceWithKnots.hxx>
  19 #include <StepData_StepReaderData.hxx>
  20 #include <StepData_StepWriter.hxx>
  21 #include <StepGeom_BSplineSurfaceWithKnots.hxx>
  22 #include <StepGeom_CartesianPoint.hxx>
  23 #include <StepGeom_HArray2OfCartesianPoint.hxx>
  24 #include <StepGeom_KnotType.hxx>
  25 #include <TCollection_AsciiString.hxx>
  26 #include <TColStd_HArray1OfInteger.hxx>
  27 #include <TColStd_HArray1OfReal.hxx>
  28 
  29 // --- Enum : KnotType ---
  30 static TCollection_AsciiString ktUniformKnots(".UNIFORM_KNOTS.");
  31 static TCollection_AsciiString ktQuasiUniformKnots(".QUASI_UNIFORM_KNOTS.");
  32 static TCollection_AsciiString ktPiecewiseBezierKnots(".PIECEWISE_BEZIER_KNOTS.");
  33 static TCollection_AsciiString ktUnspecified(".UNSPECIFIED.");
  34 
  35 
  36         // --- Enum : BSplineSurfaceForm ---
  37 static TCollection_AsciiString bssfSurfOfLinearExtrusion(".SURF_OF_LINEAR_EXTRUSION.");
  38 static TCollection_AsciiString bssfPlaneSurf(".PLANE_SURF.");
  39 static TCollection_AsciiString bssfGeneralisedCone(".GENERALISED_CONE.");
  40 static TCollection_AsciiString bssfToroidalSurf(".TOROIDAL_SURF.");
  41 static TCollection_AsciiString bssfConicalSurf(".CONICAL_SURF.");
  42 static TCollection_AsciiString bssfSphericalSurf(".SPHERICAL_SURF.");
  43 static TCollection_AsciiString bssfUnspecified(".UNSPECIFIED.");
  44 static TCollection_AsciiString bssfRuledSurf(".RULED_SURF.");
  45 static TCollection_AsciiString bssfSurfOfRevolution(".SURF_OF_REVOLUTION.");
  46 static TCollection_AsciiString bssfCylindricalSurf(".CYLINDRICAL_SURF.");
  47 static TCollection_AsciiString bssfQuadricSurf(".QUADRIC_SURF.");
  48 
  49 RWStepGeom_RWBSplineSurfaceWithKnots::RWStepGeom_RWBSplineSurfaceWithKnots () {}
  50 
  51 void RWStepGeom_RWBSplineSurfaceWithKnots::ReadStep
  52         (const Handle(StepData_StepReaderData)& data,
  53          const Standard_Integer num,
  54          Handle(Interface_Check)& ach,
  55          const Handle(StepGeom_BSplineSurfaceWithKnots)& ent) const
  56 {
  57 
  58 
  59         // --- Number of Parameter Control ---
  60 
  61         if (!data->CheckNbParams(num,13,ach,"b_spline_surface_with_knots")) return;
  62 
  63         // --- inherited field : name ---
  64 
  65         Handle(TCollection_HAsciiString) aName;
  66         //szv#4:S4163:12Mar99 `Standard_Boolean stat1 =` not needed
  67         data->ReadString (num,1,"name",ach,aName);
  68 
  69         // --- inherited field : uDegree ---
  70 
  71         Standard_Integer aUDegree;
  72         //szv#4:S4163:12Mar99 `Standard_Boolean stat2 =` not needed
  73         data->ReadInteger (num,2,"u_degree",ach,aUDegree);
  74 
  75         // --- inherited field : vDegree ---
  76 
  77         Standard_Integer aVDegree;
  78         //szv#4:S4163:12Mar99 `Standard_Boolean stat3 =` not needed
  79         data->ReadInteger (num,3,"v_degree",ach,aVDegree);
  80 
  81         // --- inherited field : controlPointsList ---
  82 
  83         Handle(StepGeom_HArray2OfCartesianPoint) aControlPointsList;
  84         Handle(StepGeom_CartesianPoint) anent4;
  85         Standard_Integer nsub4;
  86         if (data->ReadSubList (num,4,"control_points_list",ach,nsub4)) {
  87           Standard_Integer nbi4 = data->NbParams(nsub4);
  88           Standard_Integer nbj4 = data->NbParams(data->ParamNumber(nsub4,1));
  89           aControlPointsList = new StepGeom_HArray2OfCartesianPoint (1, nbi4, 1, nbj4);
  90           for (Standard_Integer i4 = 1; i4 <= nbi4; i4 ++) {
  91             Standard_Integer nsi4;
  92             if (data->ReadSubList (nsub4,i4,"sub-part(control_points_list)",ach,nsi4)) {
  93               for (Standard_Integer j4 =1; j4 <= nbj4; j4 ++) {
  94                 //szv#4:S4163:12Mar99 `Standard_Boolean stat4 =` not needed
  95                 if (data->ReadEntity (nsi4, j4,"cartesian_point", ach,
  96                                       STANDARD_TYPE(StepGeom_CartesianPoint), anent4))
  97                   aControlPointsList->SetValue(i4, j4, anent4);
  98               }
  99             }
 100           }
 101         }
 102 
 103         // --- inherited field : surfaceForm ---
 104 
 105         StepGeom_BSplineSurfaceForm aSurfaceForm = StepGeom_bssfPlaneSurf;
 106         if (data->ParamType(num,5) == Interface_ParamEnum) {
 107           Standard_CString text = data->ParamCValue(num,5);
 108           if      (bssfSurfOfLinearExtrusion.IsEqual(text)) aSurfaceForm = StepGeom_bssfSurfOfLinearExtrusion;
 109           else if (bssfPlaneSurf.IsEqual(text)) aSurfaceForm = StepGeom_bssfPlaneSurf;
 110           else if (bssfGeneralisedCone.IsEqual(text)) aSurfaceForm = StepGeom_bssfGeneralisedCone;
 111           else if (bssfToroidalSurf.IsEqual(text)) aSurfaceForm = StepGeom_bssfToroidalSurf;
 112           else if (bssfConicalSurf.IsEqual(text)) aSurfaceForm = StepGeom_bssfConicalSurf;
 113           else if (bssfSphericalSurf.IsEqual(text)) aSurfaceForm = StepGeom_bssfSphericalSurf;
 114           else if (bssfUnspecified.IsEqual(text)) aSurfaceForm = StepGeom_bssfUnspecified;
 115           else if (bssfRuledSurf.IsEqual(text)) aSurfaceForm = StepGeom_bssfRuledSurf;
 116           else if (bssfSurfOfRevolution.IsEqual(text)) aSurfaceForm = StepGeom_bssfSurfOfRevolution;
 117           else if (bssfCylindricalSurf.IsEqual(text)) aSurfaceForm = StepGeom_bssfCylindricalSurf;
 118           else if (bssfQuadricSurf.IsEqual(text)) aSurfaceForm = StepGeom_bssfQuadricSurf;
 119           else ach->AddFail("Enumeration b_spline_surface_form has not an allowed value");
 120         }
 121         else ach->AddFail("Parameter #5 (surface_form) is not an enumeration");
 122 
 123         // --- inherited field : uClosed ---
 124 
 125         StepData_Logical aUClosed;
 126         //szv#4:S4163:12Mar99 `Standard_Boolean stat6 =` not needed
 127         data->ReadLogical (num,6,"u_closed",ach,aUClosed);
 128 
 129         // --- inherited field : vClosed ---
 130 
 131         StepData_Logical aVClosed;
 132         //szv#4:S4163:12Mar99 `Standard_Boolean stat7 =` not needed
 133         data->ReadLogical (num,7,"v_closed",ach,aVClosed);
 134 
 135         // --- inherited field : selfIntersect ---
 136 
 137         StepData_Logical aSelfIntersect;
 138         //szv#4:S4163:12Mar99 `Standard_Boolean stat8 =` not needed
 139         data->ReadLogical (num,8,"self_intersect",ach,aSelfIntersect);
 140 
 141         // --- own field : uMultiplicities ---
 142 
 143         Handle(TColStd_HArray1OfInteger) aUMultiplicities;
 144         Standard_Integer aUMultiplicitiesItem;
 145         Standard_Integer nsub9;
 146         if (data->ReadSubList (num,9,"u_multiplicities",ach,nsub9)) {
 147           Standard_Integer nb9 = data->NbParams(nsub9);
 148           aUMultiplicities = new TColStd_HArray1OfInteger (1, nb9);
 149           for (Standard_Integer i9 = 1; i9 <= nb9; i9 ++) {
 150             //szv#4:S4163:12Mar99 `Standard_Boolean stat9 =` not needed
 151             if (data->ReadInteger (nsub9,i9,"u_multiplicities",ach,aUMultiplicitiesItem))
 152               aUMultiplicities->SetValue(i9,aUMultiplicitiesItem);
 153           }
 154         }
 155 
 156         // --- own field : vMultiplicities ---
 157 
 158         Handle(TColStd_HArray1OfInteger) aVMultiplicities;
 159         Standard_Integer aVMultiplicitiesItem;
 160         Standard_Integer nsub10;
 161         if (data->ReadSubList (num,10,"v_multiplicities",ach,nsub10)) {
 162           Standard_Integer nb10 = data->NbParams(nsub10);
 163           aVMultiplicities = new TColStd_HArray1OfInteger (1, nb10);
 164           for (Standard_Integer i10 = 1; i10 <= nb10; i10 ++) {
 165             //szv#4:S4163:12Mar99 `Standard_Boolean stat10 =` not needed
 166             if (data->ReadInteger (nsub10,i10,"v_multiplicities",ach,aVMultiplicitiesItem))
 167               aVMultiplicities->SetValue(i10,aVMultiplicitiesItem);
 168           }
 169         }
 170 
 171         // --- own field : uKnots ---
 172 
 173         Handle(TColStd_HArray1OfReal) aUKnots;
 174         Standard_Real aUKnotsItem;
 175         Standard_Integer nsub11;
 176         if (data->ReadSubList (num,11,"u_knots",ach,nsub11)) {
 177           Standard_Integer nb11 = data->NbParams(nsub11);
 178           aUKnots = new TColStd_HArray1OfReal (1, nb11);
 179           for (Standard_Integer i11 = 1; i11 <= nb11; i11 ++) {
 180             //szv#4:S4163:12Mar99 `Standard_Boolean stat11 =` not needed
 181             if (data->ReadReal (nsub11,i11,"u_knots",ach,aUKnotsItem))
 182               aUKnots->SetValue(i11,aUKnotsItem);
 183           }
 184         }
 185 
 186         // --- own field : vKnots ---
 187 
 188         Handle(TColStd_HArray1OfReal) aVKnots;
 189         Standard_Real aVKnotsItem;
 190         Standard_Integer nsub12;
 191         if (data->ReadSubList (num,12,"v_knots",ach,nsub12)) {
 192           Standard_Integer nb12 = data->NbParams(nsub12);
 193           aVKnots = new TColStd_HArray1OfReal (1, nb12);
 194           for (Standard_Integer i12 = 1; i12 <= nb12; i12 ++) {
 195             //szv#4:S4163:12Mar99 `Standard_Boolean stat12 =` not needed
 196             if (data->ReadReal (nsub12,i12,"v_knots",ach,aVKnotsItem))
 197               aVKnots->SetValue(i12,aVKnotsItem);
 198           }
 199         }
 200 
 201         // --- own field : knotSpec ---
 202 
 203         StepGeom_KnotType aKnotSpec = StepGeom_ktUniformKnots;
 204         if (data->ParamType(num,13) == Interface_ParamEnum) {
 205           Standard_CString text = data->ParamCValue(num,13);
 206           if      (ktUniformKnots.IsEqual(text)) aKnotSpec = StepGeom_ktUniformKnots;
 207           else if (ktQuasiUniformKnots.IsEqual(text)) aKnotSpec = StepGeom_ktQuasiUniformKnots;
 208           else if (ktPiecewiseBezierKnots.IsEqual(text)) aKnotSpec = StepGeom_ktPiecewiseBezierKnots;
 209           else if (ktUnspecified.IsEqual(text)) aKnotSpec = StepGeom_ktUnspecified;
 210           else ach->AddFail("Enumeration knot_type has not an allowed value");
 211         }
 212         else ach->AddFail("Parameter #13 (knot_spec) is not an enumeration");
 213 
 214         //--- Initialisation of the read entity ---
 215 
 216 
 217         ent->Init(aName, aUDegree, aVDegree, aControlPointsList, aSurfaceForm, aUClosed, aVClosed, aSelfIntersect, aUMultiplicities, aVMultiplicities, aUKnots, aVKnots, aKnotSpec);
 218 }
 219 
 220 
 221 void RWStepGeom_RWBSplineSurfaceWithKnots::WriteStep
 222         (StepData_StepWriter& SW,
 223          const Handle(StepGeom_BSplineSurfaceWithKnots)& ent) const
 224 {
 225 
 226         // --- inherited field name ---
 227 
 228         SW.Send(ent->Name());
 229 
 230         // --- inherited field uDegree ---
 231 
 232         SW.Send(ent->UDegree());
 233 
 234         // --- inherited field vDegree ---
 235 
 236         SW.Send(ent->VDegree());
 237 
 238         // --- inherited field controlPointsList ---
 239 
 240         SW.OpenSub();
 241         for (Standard_Integer i4 = 1;  i4 <= ent->NbControlPointsListI(); i4 ++) {
 242           SW.NewLine(Standard_False);
 243           SW.OpenSub();
 244           for (Standard_Integer j4 = 1;  j4 <= ent->NbControlPointsListJ(); j4 ++) {
 245             SW.Send(ent->ControlPointsListValue(i4,j4));
 246             SW.JoinLast(Standard_False);
 247           }
 248           SW.CloseSub();
 249         }
 250         SW.CloseSub();
 251 
 252         // --- inherited field surfaceForm ---
 253 
 254         switch(ent->SurfaceForm()) {
 255           case StepGeom_bssfSurfOfLinearExtrusion : SW.SendEnum (bssfSurfOfLinearExtrusion); break;
 256           case StepGeom_bssfPlaneSurf : SW.SendEnum (bssfPlaneSurf); break;
 257           case StepGeom_bssfGeneralisedCone : SW.SendEnum (bssfGeneralisedCone); break;
 258           case StepGeom_bssfToroidalSurf : SW.SendEnum (bssfToroidalSurf); break;
 259           case StepGeom_bssfConicalSurf : SW.SendEnum (bssfConicalSurf); break;
 260           case StepGeom_bssfSphericalSurf : SW.SendEnum (bssfSphericalSurf); break;
 261           case StepGeom_bssfUnspecified : SW.SendEnum (bssfUnspecified); break;
 262           case StepGeom_bssfRuledSurf : SW.SendEnum (bssfRuledSurf); break;
 263           case StepGeom_bssfSurfOfRevolution : SW.SendEnum (bssfSurfOfRevolution); break;
 264           case StepGeom_bssfCylindricalSurf : SW.SendEnum (bssfCylindricalSurf); break;
 265           case StepGeom_bssfQuadricSurf : SW.SendEnum (bssfQuadricSurf); break;
 266         }
 267 
 268         // --- inherited field uClosed ---
 269 
 270         SW.SendLogical(ent->UClosed());
 271 
 272         // --- inherited field vClosed ---
 273 
 274         SW.SendLogical(ent->VClosed());
 275 
 276         // --- inherited field selfIntersect ---
 277 
 278         SW.SendLogical(ent->SelfIntersect());
 279 
 280         // --- own field : uMultiplicities ---
 281 
 282         SW.OpenSub();
 283         for (Standard_Integer i9 = 1;  i9 <= ent->NbUMultiplicities();  i9 ++) {
 284           SW.Send(ent->UMultiplicitiesValue(i9));
 285         }
 286         SW.CloseSub();
 287 
 288         // --- own field : vMultiplicities ---
 289 
 290         SW.OpenSub();
 291         for (Standard_Integer i10 = 1;  i10 <= ent->NbVMultiplicities();  i10 ++) {
 292           SW.Send(ent->VMultiplicitiesValue(i10));
 293         }
 294         SW.CloseSub();
 295 
 296         // --- own field : uKnots ---
 297 
 298         SW.OpenSub();
 299         for (Standard_Integer i11 = 1;  i11 <= ent->NbUKnots();  i11 ++) {
 300           SW.Send(ent->UKnotsValue(i11));
 301         }
 302         SW.CloseSub();
 303 
 304         // --- own field : vKnots ---
 305 
 306         SW.OpenSub();
 307         for (Standard_Integer i12 = 1;  i12 <= ent->NbVKnots();  i12 ++) {
 308           SW.Send(ent->VKnotsValue(i12));
 309         }
 310         SW.CloseSub();
 311 
 312         // --- own field : knotSpec ---
 313 
 314         switch(ent->KnotSpec()) {
 315           case StepGeom_ktUniformKnots : SW.SendEnum (ktUniformKnots); break;
 316           case StepGeom_ktQuasiUniformKnots : SW.SendEnum (ktQuasiUniformKnots); break;
 317           case StepGeom_ktPiecewiseBezierKnots : SW.SendEnum (ktPiecewiseBezierKnots); break;
 318           case StepGeom_ktUnspecified : SW.SendEnum (ktUnspecified); break;
 319         }
 320 }
 321 
 322 
 323 void RWStepGeom_RWBSplineSurfaceWithKnots::Share(const Handle(StepGeom_BSplineSurfaceWithKnots)& ent, Interface_EntityIterator& iter) const
 324 {
 325 
 326         Standard_Integer nbiElem1 = ent->NbControlPointsListI();
 327         Standard_Integer nbjElem1 = ent->NbControlPointsListJ();
 328         for (Standard_Integer is1=1; is1<=nbiElem1; is1 ++) {
 329           for (Standard_Integer js1=1; js1<=nbjElem1; js1 ++) {
 330             iter.GetOneItem(ent->ControlPointsListValue(is1,js1));
 331           }
 332         }
 333 
 334 }
 335 
 336 
 337 
 338 void RWStepGeom_RWBSplineSurfaceWithKnots::Check
 339   (const Handle(StepGeom_BSplineSurfaceWithKnots)& ent,
 340    const Interface_ShareTool& ,
 341    Handle(Interface_Check)& ach) const
 342 {
 343   Standard_Integer nbCPLU  = ent->NbControlPointsListI();
 344   Standard_Integer nbCPLV  = ent->NbControlPointsListJ();
 345   Standard_Integer dgBSSU  = ent->UDegree();
 346   Standard_Integer dgBSSV  = ent->VDegree();
 347   Standard_Integer nbMulU  = ent->NbUMultiplicities();
 348   Standard_Integer nbMulV  = ent->NbVMultiplicities();
 349   Standard_Integer nbKnoU  = ent->NbUKnots();
 350   Standard_Integer nbKnoV  = ent->NbVKnots();
 351   Standard_Integer sumMulU = 0;
 352   Standard_Integer sumMulV = 0;
 353   Standard_Integer i;
 354 //  std::cout << "BSplineSurfaceWithKnots: nbMulU=" << nbMulU << " nbKnoU= " << 
 355 //    nbKnoU << " nbCPLU= " << nbCPLU << " degreeU= " << dgBSSU << std::endl;
 356 //  std::cout << "                         nbMulV=" << nbMulV << " nbKnoV= " << 
 357 //    nbKnoV << " nbCPLV= " << nbCPLV << " degreeV= " << dgBSSV << std::endl;
 358   if(nbMulU != nbKnoU) {
 359     ach->AddFail("ERROR: No.of KnotMultiplicities not equal No.of Knots in U");
 360   }
 361   if(nbMulV != nbKnoV) {
 362     ach->AddFail("ERROR: No.of KnotMultiplicities not equal No.of Knots in V");
 363   }
 364 
 365   // check in U direction
 366 
 367   for(i=1; i<=nbMulU-1; i++) {
 368     sumMulU = sumMulU + ent->UMultiplicitiesValue(i);
 369   }
 370   Standard_Integer sumNonPU = nbCPLU + dgBSSU + 1;
 371   Standard_Integer mult1U = ent->UMultiplicitiesValue(1);
 372   Standard_Integer multNU = ent->UMultiplicitiesValue(nbMulU);
 373 //  std::cout << "BSplineSurfaceWithKnots: mult1U=" << mult1U << " multNU= " <<
 374 //    multNU << " sumMulU= " << sumMulU << std::endl;
 375   if((sumMulU + multNU) == sumNonPU) {
 376   }
 377   else if((sumMulU == nbCPLU) && (mult1U == multNU)) {
 378   }
 379   else {
 380     ach->AddFail("ERROR: wrong number of Knot Multiplicities in U");
 381   }
 382   for(i=2; i<=nbKnoU; i++) {
 383     Standard_Real distKn  = ent->UKnotsValue(i-1) - ent->UKnotsValue(i);
 384     if(Abs(distKn) <= RealEpsilon())
 385       ach->AddWarning("WARNING: Surface contains identical KnotsValues in U");
 386     else if(distKn > RealEpsilon())
 387       ach->AddFail("ERROR: Surface contains descending KnotsValues in U");
 388   }
 389 
 390   // check in V direction
 391 
 392   for(i=1; i<=nbMulV-1; i++) {
 393     sumMulV = sumMulV + ent->VMultiplicitiesValue(i);
 394   }
 395   Standard_Integer sumNonPV = nbCPLV + dgBSSV + 1;
 396   Standard_Integer mult1V = ent->VMultiplicitiesValue(1);
 397   Standard_Integer multNV = ent->VMultiplicitiesValue(nbMulV);
 398 //  std::cout << "                       : mult1V=" << mult1V << " multNV= " <<
 399 //    multNV << " sumMulV= " << sumMulV << std::endl;
 400   if((sumMulV + multNV) == sumNonPV) {
 401   }
 402   else if((sumMulV == nbCPLV) && (mult1V == multNV)) {
 403   }
 404   else {
 405     ach->AddFail("ERROR: wrong number of Knot Multiplicities in V");
 406   }
 407   for(i=2; i<=nbKnoV; i++) {
 408     Standard_Real distKn  = ent->VKnotsValue(i-1) - ent->VKnotsValue(i);
 409     if(Abs(distKn) <= RealEpsilon())
 410       ach->AddWarning("WARNING: Surface contains identical KnotsValues in V");
 411     else if(distKn > RealEpsilon())
 412       ach->AddFail("ERROR: Surface contains descending KnotsValues in V");
 413   }
 414 }
Open CASCADE Technology repositoryRSSAtom
