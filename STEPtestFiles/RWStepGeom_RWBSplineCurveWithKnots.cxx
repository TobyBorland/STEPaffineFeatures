gitprojects / occt.git / blob

commit
 ? search: 
 re
summary | shortlog | log | commit | commitdiff | tree
blame | history | raw | HEAD
0033261: Data Exchange, Step Import - Empty shape after reading process
[occt.git] / src / RWStepGeom / RWStepGeom_RWBSplineCurveWithKnots.cxx
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
  18 #include <RWStepGeom_RWBSplineCurveWithKnots.hxx>
  19 #include <StepData_StepReaderData.hxx>
  20 #include <StepData_StepWriter.hxx>
  21 #include <StepGeom_BSplineCurveWithKnots.hxx>
  22 #include <StepGeom_CartesianPoint.hxx>
  23 #include <StepGeom_HArray1OfCartesianPoint.hxx>
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
  36         // --- Enum : BSplineCurveForm ---
  37 static TCollection_AsciiString bscfEllipticArc(".ELLIPTIC_ARC.");
  38 static TCollection_AsciiString bscfPolylineForm(".POLYLINE_FORM.");
  39 static TCollection_AsciiString bscfParabolicArc(".PARABOLIC_ARC.");
  40 static TCollection_AsciiString bscfCircularArc(".CIRCULAR_ARC.");
  41 static TCollection_AsciiString bscfUnspecified(".UNSPECIFIED.");
  42 static TCollection_AsciiString bscfHyperbolicArc(".HYPERBOLIC_ARC.");
  43 
  44 RWStepGeom_RWBSplineCurveWithKnots::RWStepGeom_RWBSplineCurveWithKnots () {}
  45 
  46 void RWStepGeom_RWBSplineCurveWithKnots::ReadStep
  47         (const Handle(StepData_StepReaderData)& data,
  48          const Standard_Integer num,
  49          Handle(Interface_Check)& ach,
  50          const Handle(StepGeom_BSplineCurveWithKnots)& ent) const
  51 {
  52 
  53 
  54         // --- Number of Parameter Control ---
  55 
  56         if (!data->CheckNbParams(num,9,ach,"b_spline_curve_with_knots")) return;
  57 
  58         // --- inherited field : name ---
  59 
  60         Handle(TCollection_HAsciiString) aName;
  61         //szv#4:S4163:12Mar99 `Standard_Boolean stat1 =` not needed
  62         data->ReadString (num,1,"name",ach,aName);
  63 
  64         // --- inherited field : degree ---
  65 
  66         Standard_Integer aDegree;
  67         //szv#4:S4163:12Mar99 `Standard_Boolean stat2 =` not needed
  68         data->ReadInteger (num,2,"degree",ach,aDegree);
  69 
  70         // --- inherited field : controlPointsList ---
  71 
  72         Handle(StepGeom_HArray1OfCartesianPoint) aControlPointsList;
  73         Handle(StepGeom_CartesianPoint) anent3;
  74         Standard_Integer nsub3;
  75         if (data->ReadSubList (num,3,"control_points_list",ach,nsub3)) {
  76           Standard_Integer nb3 = data->NbParams(nsub3);
  77     if(nb3 <1)
  78       ach->AddFail("Number of control points of the b_spline_curve_form is equal to 0");
  79     else
  80     {
  81       aControlPointsList = new StepGeom_HArray1OfCartesianPoint (1, nb3);
  82       for (Standard_Integer i3 = 1; i3 <= nb3; i3 ++) {
  83         //szv#4:S4163:12Mar99 `Standard_Boolean stat3 =` not needed
  84         if (data->ReadEntity (nsub3, i3,"cartesian_point", ach,
  85           STANDARD_TYPE(StepGeom_CartesianPoint), anent3))
  86           aControlPointsList->SetValue(i3, anent3);
  87       }
  88     }
  89         }
  90 
  91         // --- inherited field : curveForm ---
  92 
  93         StepGeom_BSplineCurveForm aCurveForm = StepGeom_bscfPolylineForm;
  94         if (data->ParamType(num,4) == Interface_ParamEnum) {
  95           Standard_CString text = data->ParamCValue(num,4);
  96           if      (bscfEllipticArc.IsEqual(text)) aCurveForm = StepGeom_bscfEllipticArc;
  97           else if (bscfPolylineForm.IsEqual(text)) aCurveForm = StepGeom_bscfPolylineForm;
  98           else if (bscfParabolicArc.IsEqual(text)) aCurveForm = StepGeom_bscfParabolicArc;
  99           else if (bscfCircularArc.IsEqual(text)) aCurveForm = StepGeom_bscfCircularArc;
 100           else if (bscfUnspecified.IsEqual(text)) aCurveForm = StepGeom_bscfUnspecified;
 101           else if (bscfHyperbolicArc.IsEqual(text)) aCurveForm = StepGeom_bscfHyperbolicArc;
 102           else ach->AddFail("Enumeration b_spline_curve_form has not an allowed value");
 103         }
 104         else ach->AddFail("Parameter #4 (curve_form) is not an enumeration");
 105 
 106         // --- inherited field : closedCurve ---
 107 
 108         StepData_Logical aClosedCurve;
 109         //szv#4:S4163:12Mar99 `Standard_Boolean stat5 =` not needed
 110         data->ReadLogical (num,5,"closed_curve",ach,aClosedCurve);
 111 
 112         // --- inherited field : selfIntersect ---
 113 
 114         StepData_Logical aSelfIntersect;
 115         //szv#4:S4163:12Mar99 `Standard_Boolean stat6 =` not needed
 116         data->ReadLogical (num,6,"self_intersect",ach,aSelfIntersect);
 117 
 118         // --- own field : knotMultiplicities ---
 119 
 120         Handle(TColStd_HArray1OfInteger) aKnotMultiplicities;
 121         Standard_Integer aKnotMultiplicitiesItem;
 122         Standard_Integer nsub7;
 123         if (data->ReadSubList (num,7,"knot_multiplicities",ach,nsub7)) {
 124           Standard_Integer nb7 = data->NbParams(nsub7);
 125           aKnotMultiplicities = new TColStd_HArray1OfInteger (1, nb7);
 126           for (Standard_Integer i7 = 1; i7 <= nb7; i7 ++) {
 127             //szv#4:S4163:12Mar99 `Standard_Boolean stat7 =` not needed
 128             if (data->ReadInteger (nsub7,i7,"knot_multiplicities",ach,aKnotMultiplicitiesItem))
 129               aKnotMultiplicities->SetValue(i7,aKnotMultiplicitiesItem);
 130           }
 131         }
 132 
 133         // --- own field : knots ---
 134 
 135         Handle(TColStd_HArray1OfReal) aKnots;
 136         Standard_Real aKnotsItem;
 137         Standard_Integer nsub8;
 138         if (data->ReadSubList (num,8,"knots",ach,nsub8)) {
 139           Standard_Integer nb8 = data->NbParams(nsub8);
 140           aKnots = new TColStd_HArray1OfReal (1, nb8);
 141           for (Standard_Integer i8 = 1; i8 <= nb8; i8 ++) {
 142             //szv#4:S4163:12Mar99 `Standard_Boolean stat8 =` not needed
 143             if (data->ReadReal (nsub8,i8,"knots",ach,aKnotsItem))
 144               aKnots->SetValue(i8,aKnotsItem);
 145           }
 146         }
 147 
 148         // --- own field : knotSpec ---
 149 
 150         StepGeom_KnotType aKnotSpec = StepGeom_ktUniformKnots;
 151         if (data->ParamType(num,9) == Interface_ParamEnum) {
 152           Standard_CString text = data->ParamCValue(num,9);
 153           if      (ktUniformKnots.IsEqual(text)) aKnotSpec = StepGeom_ktUniformKnots;
 154           else if (ktQuasiUniformKnots.IsEqual(text)) aKnotSpec = StepGeom_ktQuasiUniformKnots;
 155           else if (ktPiecewiseBezierKnots.IsEqual(text)) aKnotSpec = StepGeom_ktPiecewiseBezierKnots;
 156           else if (ktUnspecified.IsEqual(text)) aKnotSpec = StepGeom_ktUnspecified;
 157           else ach->AddFail("Enumeration knot_type has not an allowed value");
 158         }
 159         else ach->AddFail("Parameter #9 (knot_spec) is not an enumeration");
 160 
 161         //--- Initialisation of the read entity ---
 162 
 163 
 164         ent->Init(aName, aDegree, aControlPointsList, aCurveForm, aClosedCurve, aSelfIntersect, aKnotMultiplicities, aKnots, aKnotSpec);
 165 }
 166 
 167 
 168 void RWStepGeom_RWBSplineCurveWithKnots::WriteStep
 169         (StepData_StepWriter& SW,
 170          const Handle(StepGeom_BSplineCurveWithKnots)& ent) const
 171 {
 172 
 173         // --- inherited field name ---
 174 
 175         SW.Send(ent->Name());
 176 
 177         // --- inherited field degree ---
 178 
 179         SW.Send(ent->Degree());
 180 
 181         // --- inherited field controlPointsList ---
 182 
 183         SW.OpenSub();
 184         for (Standard_Integer i3 = 1;  i3 <= ent->NbControlPointsList();  i3 ++) {
 185           SW.Send(ent->ControlPointsListValue(i3));
 186         }
 187         SW.CloseSub();
 188 
 189         // --- inherited field curveForm ---
 190 
 191         switch(ent->CurveForm()) {
 192           case StepGeom_bscfEllipticArc : SW.SendEnum (bscfEllipticArc); break;
 193           case StepGeom_bscfPolylineForm : SW.SendEnum (bscfPolylineForm); break;
 194           case StepGeom_bscfParabolicArc : SW.SendEnum (bscfParabolicArc); break;
 195           case StepGeom_bscfCircularArc : SW.SendEnum (bscfCircularArc); break;
 196           case StepGeom_bscfUnspecified : SW.SendEnum (bscfUnspecified); break;
 197           case StepGeom_bscfHyperbolicArc : SW.SendEnum (bscfHyperbolicArc); break;
 198         }
 199 
 200         // --- inherited field closedCurve ---
 201 
 202         SW.SendLogical(ent->ClosedCurve());
 203 
 204         // --- inherited field selfIntersect ---
 205 
 206         SW.SendLogical(ent->SelfIntersect());
 207 
 208         // --- own field : knotMultiplicities ---
 209 
 210         SW.OpenSub();
 211         for (Standard_Integer i7 = 1;  i7 <= ent->NbKnotMultiplicities();  i7 ++) {
 212           SW.Send(ent->KnotMultiplicitiesValue(i7));
 213         }
 214         SW.CloseSub();
 215 
 216         // --- own field : knots ---
 217 
 218         SW.OpenSub();
 219         for (Standard_Integer i8 = 1;  i8 <= ent->NbKnots();  i8 ++) {
 220           SW.Send(ent->KnotsValue(i8));
 221         }
 222         SW.CloseSub();
 223 
 224         // --- own field : knotSpec ---
 225 
 226         switch(ent->KnotSpec()) {
 227           case StepGeom_ktUniformKnots : SW.SendEnum (ktUniformKnots); break;
 228           case StepGeom_ktQuasiUniformKnots : SW.SendEnum (ktQuasiUniformKnots); break;
 229           case StepGeom_ktPiecewiseBezierKnots : SW.SendEnum (ktPiecewiseBezierKnots); break;
 230           case StepGeom_ktUnspecified : SW.SendEnum (ktUnspecified); break;
 231         }
 232 }
 233 
 234 
 235 void RWStepGeom_RWBSplineCurveWithKnots::Share(const Handle(StepGeom_BSplineCurveWithKnots)& ent, Interface_EntityIterator& iter) const
 236 {
 237 
 238         Standard_Integer nbElem1 = ent->NbControlPointsList();
 239         for (Standard_Integer is1=1; is1<=nbElem1; is1 ++) {
 240           iter.GetOneItem(ent->ControlPointsListValue(is1));
 241         }
 242 
 243 }
 244 
 245 void RWStepGeom_RWBSplineCurveWithKnots::Check
 246 (const Handle(StepGeom_BSplineCurveWithKnots)& ent,
 247  const Interface_ShareTool& ,
 248  Handle(Interface_Check)& ach) const
 249 {
 250   Standard_Integer nbCPL  = ent->NbControlPointsList();
 251   Standard_Integer dgBSC  = ent->Degree();
 252   Standard_Integer nbMult = ent->NbKnotMultiplicities();
 253   Standard_Integer nbKno  = ent->NbKnots();
 254   Standard_Integer sumMult = 0;
 255 //  std::cout << "BSplineCurveWithKnots: nbMult=" << nbMult << " nbKno= " << 
 256 //    nbKno << " nbCPL= " << nbCPL << " degree= " << dgBSC << std::endl;
 257   if(nbMult != nbKno) {
 258     ach->AddFail("ERROR: No.of KnotMultiplicities not equal No.of Knots");
 259   }
 260   Standard_Integer i;//svv Jan 10 2000: porting on DEC 
 261   for (i=1; i<=nbMult-1; i++) {
 262     sumMult = sumMult + ent->KnotMultiplicitiesValue(i);
 263   }
 264   Standard_Integer sumNonP = nbCPL + dgBSC + 1;
 265   Standard_Integer mult1 = ent->KnotMultiplicitiesValue(1);
 266   Standard_Integer multN = ent->KnotMultiplicitiesValue(nbMult);
 267 //  std::cout << "BSplineCurveWithKnots: mult1=" << mult1 << " multN= " <<
 268 //    multN << " sumMult= " << sumMult << std::endl;
 269   if((sumMult + multN) == sumNonP) {
 270   }
 271   else if((sumMult == nbCPL) && (mult1 == multN)) {
 272   }
 273   else {
 274     ach->AddFail("ERROR: wrong number of Knot Multiplicities");
 275   }
 276   for(i=2; i<=nbKno; i++) {
 277     Standard_Real distKn  = ent->KnotsValue(i-1) - ent->KnotsValue(i);
 278     if(Abs(distKn) <= RealEpsilon())
 279       ach->AddWarning("WARNING: Curve contains identical KnotsValues");
 280     else if(distKn > RealEpsilon())
 281       ach->AddFail("ERROR: Curve contains descending KnotsValues");
 282   }
 283 }
 284 
Open CASCADE Technology repositoryRSSAtom
