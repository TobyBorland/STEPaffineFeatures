import FreeCAD as App
doc = App.newDocument("ConeDoc")
#doc = App.ActiveDocument
body = doc.addObject('PartDesign::Body','Body')
cone = doc.addObject('PartDesign::AdditiveCone','Cone')
bc = body.addObject(cone)
doc.recompute()
cone.Radius1 = 1
cone.Radius2 = 0.5
cone.Height = 1
dir = '/home/foobert/STEP_test_files'
Part.export([body],dir + 'Standard_Cone.step')