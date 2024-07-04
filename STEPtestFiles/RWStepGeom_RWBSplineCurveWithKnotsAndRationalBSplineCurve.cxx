gitprojects / occt.git / blob

commit
 ? search: 
 re
summary | shortlog | log | commit | commitdiff | tree
blame | history | raw | HEAD
0033261: Data Exchange, Step Import - Empty shape after reading process
[occt.git] / src / RWStepGeom / RWStepGeom_RWBSplineCurveWithKnotsAndRationalBSplineCurve.cxx
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
  14 // sln 04.10.2001. BUC61003. Correction of looking for items of complex entity
  15 
  16 #include <Interface_Check.hxx>
  17 #include <Interface_EntityIterator.hxx>
  18 #include <Interface_ShareTool.hxx>
  19 #include <RWStepGeom_RWBSplineCurveWithKnots.hxx>
  20 #include <RWStepGeom_RWBSplineCurveWithKnotsAndRationalBSplineCurve.hxx>
  21 #include <RWStepGeom_RWRationalBSplineCurve.hxx>
  22 #include <StepData_StepReaderData.hxx>
  23 #include <StepData_StepWriter.hxx>
  24 #include <StepGeom_BSplineCurveWithKnots.hxx>
  25 #include <StepGeom_BSplineCurveWithKnotsAndRationalBSplineCurve.hxx>
  26 #include <StepGeom_CartesianPoint.hxx>
  27 #include <StepGeom_HArray1OfCartesianPoint.hxx>
  28 #include <StepGeom_KnotType.hxx>
  29 #include <StepGeom_RationalBSplineCurve.hxx>
  30 #include <TColStd_HArray1OfInteger.hxx>
  31 #include <TColStd_HArray1OfReal.hxx>
  32 
  33 // --- Enum : BSplineCurveForm ---
  34 static TCollection_AsciiString bscfEllipticArc(".ELLIPTIC_ARC.");
  35 static TCollection_AsciiString bscfPolylineForm(".POLYLINE_FORM.");
  36 static TCollection_AsciiString bscfParabolicArc(".PARABOLIC_ARC.");
  37 static TCollection_AsciiString bscfCircularArc(".CIRCULAR_ARC.");
  38 static TCollection_AsciiString bscfUnspecified(".UNSPECIFIED.");
  39 static TCollection_AsciiString bscfHyperbolicArc(".HYPERBOLIC_ARC.");
  40 
  41         // --- Enum : KnotType ---
  42 static TCollection_AsciiString ktUniformKnots(".UNIFORM_KNOTS.");
  43 static TCollection_AsciiString ktQuasiUniformKnots(".QUASI_UNIFORM_KNOTS.");
  44 static TCollection_AsciiString ktPiecewiseBezierKnots(".PIECEWISE_BEZIER_KNOTS.");
  45 static TCollection_AsciiString ktUnspecified(".UNSPECIFIED.");
  46 
  47 RWStepGeom_RWBSplineCurveWithKnotsAndRationalBSplineCurve::RWStepGeom_RWBSplineCurveWithKnotsAndRationalBSplineCurve () {}
  48 
  49 void RWStepGeom_RWBSplineCurveWithKnotsAndRationalBSplineCurve::ReadStep
  50         (const Handle(StepData_StepReaderData)& data,
  51          const Standard_Integer num0,
  52          Handle(Interface_Check)& ach,
  53          const Handle(StepGeom_BSplineCurveWithKnotsAndRationalBSplineCurve)& ent) const
  54 {
  55 
  56 // sln 04.10.2001. BUC61003. Correction of looking for items of complex entity
  57         Standard_Integer num = 0;  // num0
  58         data->NamedForComplex("BOUNDED_CURVE", "BNDCRV",num0,num,ach);
  59 
  60 //      num = data->NextForComplex(num);
  61 // sln 04.10.2001. BUC61003. Correction of looking for items of complex entity
  62 //        num =  0; gka TRJ9
  63   data->NamedForComplex("B_SPLINE_CURVE", "BSPCR",num0,num,ach);
  64 
  65         // --- Instance of common supertype BSplineCurve ---
  66 
  67         if (!data->CheckNbParams(num,5,ach,"b_spline_curve")) return;
  68         // --- field : degree ---
  69 
  70 
  71         Standard_Integer aDegree;
  72         //szv#4:S4163:12Mar99 `Standard_Boolean stat1 =` not needed
  73         data->ReadInteger (num,1,"degree",ach,aDegree);
  74         // --- field : controlPointsList ---
  75 
  76 
  77         Handle(StepGeom_HArray1OfCartesianPoint) aControlPointsList;
  78         Handle(StepGeom_CartesianPoint) anent2;
  79         Standard_Integer nsub2;
  80         if (data->ReadSubList (num,2,"control_points_list",ach,nsub2)) {
  81           Standard_Integer nb2 = data->NbParams(nsub2);
  82           aControlPointsList = new StepGeom_HArray1OfCartesianPoint (1, nb2);
  83           for (Standard_Integer i2 = 1; i2 <= nb2; i2 ++) {
  84             //szv#4:S4163:12Mar99 `Standard_Boolean stat2 =` not needed
  85             if (data->ReadEntity (nsub2, i2,"cartesian_point", ach,
  86                                   STANDARD_TYPE(StepGeom_CartesianPoint), anent2))
  87               aControlPointsList->SetValue(i2, anent2);
  88           }
  89         }
  90 
  91         // --- field : curveForm ---
  92 
  93 
  94         StepGeom_BSplineCurveForm aCurveForm = StepGeom_bscfPolylineForm;
  95         if (data->ParamType(num,3) == Interface_ParamEnum) {
  96           Standard_CString text = data->ParamCValue(num,3);
  97           if      (bscfEllipticArc.IsEqual(text)) aCurveForm = StepGeom_bscfEllipticArc;
  98           else if (bscfPolylineForm.IsEqual(text)) aCurveForm = StepGeom_bscfPolylineForm;
  99           else if (bscfParabolicArc.IsEqual(text)) aCurveForm = StepGeom_bscfParabolicArc;
 100           else if (bscfCircularArc.IsEqual(text)) aCurveForm = StepGeom_bscfCircularArc;
 101           else if (bscfUnspecified.IsEqual(text)) aCurveForm = StepGeom_bscfUnspecified;
 102           else if (bscfHyperbolicArc.IsEqual(text)) aCurveForm = StepGeom_bscfHyperbolicArc;
 103           else ach->AddFail("Enumeration b_spline_curve_form has not an allowed value");
 104         }
 105         else ach->AddFail("Parameter #3 (curve_form) is not an enumeration");
 106         // --- field : closedCurve ---
 107 
 108 
 109         StepData_Logical aClosedCurve;
 110         //szv#4:S4163:12Mar99 `Standard_Boolean stat4 =` not needed
 111         data->ReadLogical (num,4,"closed_curve",ach,aClosedCurve);
 112         // --- field : selfIntersect ---
 113 
 114 
 115         StepData_Logical aSelfIntersect;
 116         //szv#4:S4163:12Mar99 `Standard_Boolean stat5 =` not needed
 117         data->ReadLogical (num,5,"self_intersect",ach,aSelfIntersect);
 118 
 119 //      num = data->NextForComplex(num);
 120 // sln 04.10.2001. BUC61003. Correction of looking for items of complex entity
 121 //        num =  0; //gka TRJ9
 122         data->NamedForComplex("B_SPLINE_CURVE_WITH_KNOTS", "BSCWK",num0,num,ach);
 123 
 124         // --- Instance of plex component BSplineCurveWithKnots ---
 125 
 126         if (!data->CheckNbParams(num,3,ach,"b_spline_curve_with_knots")) return;
 127 
 128         // --- field : knotMultiplicities ---
 129 
 130         Handle(TColStd_HArray1OfInteger) aKnotMultiplicities;
 131         Standard_Integer aKnotMultiplicitiesItem;
 132         Standard_Integer nsub6;
 133         if (data->ReadSubList (num,1,"knot_multiplicities",ach,nsub6)) {
 134           Standard_Integer nb6 = data->NbParams(nsub6);
 135           aKnotMultiplicities = new TColStd_HArray1OfInteger (1, nb6);
 136           for (Standard_Integer i6 = 1; i6 <= nb6; i6 ++) {
 137             //szv#4:S4163:12Mar99 `Standard_Boolean stat6 =` not needed
 138             if (data->ReadInteger (nsub6,i6,"knot_multiplicities",ach,aKnotMultiplicitiesItem))
 139               aKnotMultiplicities->SetValue(i6,aKnotMultiplicitiesItem);
 140           }
 141         }
 142 
 143         // --- field : knots ---
 144 
 145         Handle(TColStd_HArray1OfReal) aKnots;
 146         Standard_Real aKnotsItem;
 147         Standard_Integer nsub7;
 148         if (data->ReadSubList (num,2,"knots",ach,nsub7)) {
 149           Standard_Integer nb7 = data->NbParams(nsub7);
 150           aKnots = new TColStd_HArray1OfReal (1, nb7);
 151           for (Standard_Integer i7 = 1; i7 <= nb7; i7 ++) {
 152             //szv#4:S4163:12Mar99 `Standard_Boolean stat7 =` not needed
 153             if (data->ReadReal (nsub7,i7,"knots",ach,aKnotsItem))
 154               aKnots->SetValue(i7,aKnotsItem);
 155           }
 156         }
 157 
 158         // --- field : knotSpec ---
 159 
 160         StepGeom_KnotType aKnotSpec = StepGeom_ktUniformKnots;
 161         if (data->ParamType(num,3) == Interface_ParamEnum) {
 162           Standard_CString text = data->ParamCValue(num,3);
 163           if      (ktUniformKnots.IsEqual(text)) aKnotSpec = StepGeom_ktUniformKnots;
 164           else if (ktQuasiUniformKnots.IsEqual(text)) aKnotSpec = StepGeom_ktQuasiUniformKnots;
 165           else if (ktPiecewiseBezierKnots.IsEqual(text)) aKnotSpec = StepGeom_ktPiecewiseBezierKnots;
 166           else if (ktUnspecified.IsEqual(text)) aKnotSpec = StepGeom_ktUnspecified;
 167           else ach->AddFail("Enumeration knot_type has not an allowed value");
 168         }
 169         else ach->AddFail("Parameter #3 (knot_spec) is not an enumeration");
 170 
 171 //      num = data->NextForComplex(num);
 172 // sln 04.10.2001. BUC61003. Correction of looking for items of complex entity
 173 //        num =  0; gka TRJ9
 174         data->NamedForComplex("CURVE",num0,num,ach);
 175 
 176 //      num = data->NextForComplex(num);
 177 // sln 04.10.2001. BUC61003. Correction of looking for items of complex entity
 178         //num =  0;
 179         data->NamedForComplex("GEOMETRIC_REPRESENTATION_ITEM", "GMRPIT",num0,num,ach);
 180 
 181 //      num = data->NextForComplex(num);
 182 // sln 04.10.2001. BUC61003. Correction of looking for items of complex entity
 183         //num =  0;
 184         data->NamedForComplex("RATIONAL_B_SPLINE_CURVE", "RBSC",num0,num,ach);
 185 
 186         // --- Instance of plex component RationalBSplineCurve ---
 187 
 188         if (!data->CheckNbParams(num,1,ach,"rational_b_spline_curve")) return;
 189 
 190         // --- field : weightsData ---
 191 
 192         Handle(TColStd_HArray1OfReal) aWeightsData;
 193         Standard_Real aWeightsDataItem;
 194         Standard_Integer nsub9;
 195         if (data->ReadSubList (num,1,"weights_data",ach,nsub9)) {
 196           Standard_Integer nb9 = data->NbParams(nsub9);
 197           aWeightsData = new TColStd_HArray1OfReal (1, nb9);
 198           for (Standard_Integer i9 = 1; i9 <= nb9; i9 ++) {
 199             //szv#4:S4163:12Mar99 `Standard_Boolean stat9 =` not needed
 200             if (data->ReadReal (nsub9,i9,"weights_data",ach,aWeightsDataItem))
 201               aWeightsData->SetValue(i9,aWeightsDataItem);
 202           }
 203         }
 204 
 205 //      num = data->NextForComplex(num);
 206 // sln 04.10.2001. BUC61003. Correction of looking for items of complex entity
 207         //num =  0;
 208         data->NamedForComplex("REPRESENTATION_ITEM", "RPRITM",num0,num,ach);
 209 
 210         // --- Instance of plex component RepresentationItem ---
 211 
 212         if (!data->CheckNbParams(num,1,ach,"representation_item")) return;
 213 
 214         // --- field : name ---
 215 
 216         Handle(TCollection_HAsciiString) aName;
 217         //szv#4:S4163:12Mar99 `Standard_Boolean stat10 =` not needed
 218         data->ReadString (num,1,"name",ach,aName);
 219 
 220         //--- Initialisation of the red entity ---
 221 
 222         ent->Init(aName,aDegree,aControlPointsList,aCurveForm,aClosedCurve,aSelfIntersect,aKnotMultiplicities,aKnots,aKnotSpec,aWeightsData);
 223 }
 224 
 225 
 226 void RWStepGeom_RWBSplineCurveWithKnotsAndRationalBSplineCurve::WriteStep
 227         (StepData_StepWriter& SW,
 228          const Handle(StepGeom_BSplineCurveWithKnotsAndRationalBSplineCurve)& ent) const
 229 {
 230 
 231         // --- Instance of plex component BoundedCurve ---
 232 
 233         SW.StartEntity("BOUNDED_CURVE");
 234 
 235         // --- Instance of common supertype BSplineCurve ---
 236 
 237         SW.StartEntity("B_SPLINE_CURVE");
 238         // --- field : degree ---
 239 
 240         SW.Send(ent->Degree());
 241         // --- field : controlPointsList ---
 242 
 243         SW.OpenSub();
 244         for (Standard_Integer i2 = 1;  i2 <= ent->NbControlPointsList();  i2 ++) {
 245           SW.Send(ent->ControlPointsListValue(i2));
 246         }
 247         SW.CloseSub();
 248         // --- field : curveForm ---
 249 
 250         switch(ent->CurveForm()) {
 251           case StepGeom_bscfEllipticArc : SW.SendEnum (bscfEllipticArc); break;
 252           case StepGeom_bscfPolylineForm : SW.SendEnum (bscfPolylineForm); break;
 253           case StepGeom_bscfParabolicArc : SW.SendEnum (bscfParabolicArc); break;
 254           case StepGeom_bscfCircularArc : SW.SendEnum (bscfCircularArc); break;
 255           case StepGeom_bscfUnspecified : SW.SendEnum (bscfUnspecified); break;
 256           case StepGeom_bscfHyperbolicArc : SW.SendEnum (bscfHyperbolicArc); break;
 257         }
 258         // --- field : closedCurve ---
 259 
 260         SW.SendLogical(ent->ClosedCurve());
 261         // --- field : selfIntersect ---
 262 
 263         SW.SendLogical(ent->SelfIntersect());
 264 
 265         // --- Instance of plex component BSplineCurveWithKnots ---
 266 
 267         SW.StartEntity("B_SPLINE_CURVE_WITH_KNOTS");
 268         // --- field : knotMultiplicities ---
 269 
 270         SW.OpenSub();
 271         for (Standard_Integer i6 = 1;  i6 <= ent->NbKnotMultiplicities();  i6 ++) {
 272           SW.Send(ent->KnotMultiplicitiesValue(i6));
 273         }
 274         SW.CloseSub();
 275         // --- field : knots ---
 276 
 277         SW.OpenSub();
 278         for (Standard_Integer i7 = 1;  i7 <= ent->NbKnots();  i7 ++) {
 279           SW.Send(ent->KnotsValue(i7));
 280         }
 281         SW.CloseSub();
 282         // --- field : knotSpec ---
 283 
 284         switch(ent->KnotSpec()) {
 285           case StepGeom_ktUniformKnots : SW.SendEnum (ktUniformKnots); break;
 286           case StepGeom_ktQuasiUniformKnots : SW.SendEnum (ktQuasiUniformKnots); break;
 287           case StepGeom_ktPiecewiseBezierKnots : SW.SendEnum (ktPiecewiseBezierKnots); break;
 288           case StepGeom_ktUnspecified : SW.SendEnum (ktUnspecified); break;
 289         }
 290 
 291         // --- Instance of plex component Curve ---
 292 
 293         SW.StartEntity("CURVE");
 294 
 295         // --- Instance of plex component GeometricRepresentationItem ---
 296 
 297         SW.StartEntity("GEOMETRIC_REPRESENTATION_ITEM");
 298 
 299         // --- Instance of plex component RationalBSplineCurve ---
 300 
 301         SW.StartEntity("RATIONAL_B_SPLINE_CURVE");
 302         // --- field : weightsData ---
 303 
 304         SW.OpenSub();
 305         for (Standard_Integer i9 = 1;  i9 <= ent->NbWeightsData();  i9 ++) {
 306           SW.Send(ent->WeightsDataValue(i9));
 307         }
 308         SW.CloseSub();
 309 
 310         // --- Instance of plex component RepresentationItem ---
 311 
 312         SW.StartEntity("REPRESENTATION_ITEM");
 313         // --- field : name ---
 314 
 315         SW.Send(ent->Name());
 316 }
 317 
 318 
 319 void RWStepGeom_RWBSplineCurveWithKnotsAndRationalBSplineCurve::Share(const Handle(StepGeom_BSplineCurveWithKnotsAndRationalBSplineCurve)& ent, Interface_EntityIterator& iter) const
 320 {
 321 
 322         Standard_Integer nbElem1 = ent->NbControlPointsList();
 323         for (Standard_Integer is1=1; is1<=nbElem1; is1 ++) {
 324           iter.GetOneItem(ent->ControlPointsListValue(is1));
 325         }
 326 
 327 }
 328 
 329 
 330 
 331 void RWStepGeom_RWBSplineCurveWithKnotsAndRationalBSplineCurve::Check
 332   (const Handle(StepGeom_BSplineCurveWithKnotsAndRationalBSplineCurve)& ent,
 333    const Interface_ShareTool& aShto,
 334    Handle(Interface_Check)& ach) const
 335 {
 336   const Handle(StepGeom_BSplineCurveWithKnotsAndRationalBSplineCurve)& aRationalBSC = ent;
 337   Handle(StepGeom_BSplineCurveWithKnots) aBSCWK =
 338     aRationalBSC->BSplineCurveWithKnots();
 339   RWStepGeom_RWBSplineCurveWithKnots t1;
 340   t1.Check(aBSCWK,aShto,ach);
 341   Handle(StepGeom_RationalBSplineCurve) aRBSC =
 342     aRationalBSC->RationalBSplineCurve();
 343   RWStepGeom_RWRationalBSplineCurve t2;
 344   t2.Check(aRBSC,aShto,ach);
 345 }
Open CASCADE Technology repositoryRSSAtom
