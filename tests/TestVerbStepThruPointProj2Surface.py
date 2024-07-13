import verb
import numpy as np
from geomdl import BSpline
from geomdl import utilities
#from geomdl.visualization import VisMPL
from geomdl import exchange
from geomdl import operations

import numpy as np
import matplotlib
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from matplotlib.patches import Ellipse, Circle
import mpl_toolkits.mplot3d.art3d as art3d

# degree = 3
# knots = [0, 0, 0, 0, 0.333, 0.666, 1, 1, 1, 1]
# pts = [ [ [0, 0, -10], 	[10, 0, 0], 	[20, 0, 0], 	[30, 0, 0] , 	[40, 0, 0], [50, 0, 0] ],
#              [ [0, -10, 0], 	[10, -10, 10], 	[20, -10, 10], 	[30, -10, 0] , [40, -10, 0], [50, -10, 0]	],
#              [ [0, -20, 0], 	[10, -20, 10], 	[20, -20, 10], 	[30, -20, 0] , [40, -20, -2], [50, -20, -12] 	],
#              [ [0, -30, 0], 	[10, -30, 0], 	[20, -30, -23], 	[30, -30, 0] , [40, -30, 0], [50, -30, 0]     ],
#              [ [0, -40, 0], 	[10, -40, 0], 	[20, -40, 0], 	[30, -40, 4] , [40, -40, -20], [50, -40, 0]     ],
#              [ [0, -50, 12], [10, -50, 0], 	[20, -50, 20], 	[30, -50, 0] , [50, -50, -10], [50, -50, -15]     ],     ]
#
# srf = verb.geom.NurbsSurface.byKnotsControlPointsWeights( degree, degree, knots, knots, pts )
#
# 	setupScene();

degree = 3
knots = [0, 0, 0, 0, 0.333, 0.666, 1, 1, 1, 1]
pts = [ [ [0, 0, -5], 	[10, 0, 0], 	[20, 0, 0], 	[30, 0, 0],     [40, 0, 0],     [50, 0, 0]      ],
		[ [0, -10, 0], 	[10, -10, 10], 	[20, -10, 10], 	[30, -10, 0],   [40, -10, 0],   [50, -10, 0]    ],
		[ [0, -20, 0], 	[10, -20, 10], 	[20, -20, 10], 	[30, -20, 0],   [40, -20, -2],  [50, -20, -12]  ],
		[ [0, -30, 0], 	[10, -30, 0], 	[20, -30, -23], [30, -30, 0],   [40, -30, 0],   [50, -30, 0]    ],
		[ [0, -40, 0], 	[10, -40, 0], 	[20, -40, 0], 	[30, -40, 4],   [40, -40, -20], [50, -40, 0]    ],
		[ [0, -50, 12], [10, -50, 0], 	[20, -50, 10], 	[30, -50, 0],   [50, -50, -3],  [50, -50, -5]   ],]

srf = verb.verb_geom_NurbsSurface.byKnotsControlPointsWeights( degree, degree, knots, knots, pts )

#lineMat = new THREE.LineBasicMaterial()

pts = []
c = 5
for i in range (0, c+1):
    for j in range(0, c + 1):
        p0 = [60 * i / (c-1) - 10, -60 * j / (c-1) + 10, 7 ]
        pts.append( p0 )
        p = srf.closestPoint( p0 )
        # l = new verb.geom.Line(p, p0)
		# addCurveToScene( l.toThreeGeometry(), lineMat )

# addPointsToScene( pts )
# addMeshToScene( srf.toThreeGeometry() )
# renderScene()

# return verb_eval_Analyze.rationalSurfaceClosestPoint(self._data,pt)
# uv = verb_eval_Analyze.rationalSurfaceClosestParam(surface,p)
# tess = verb_eval_Tess.rationalSurfaceAdaptive(surface,verb_eval_AdaptiveRefinementOptions())
# arrTrees = verb_eval_Tess.divideRationalSurfaceAdaptive(surface,options)
# x = verb_core_SurfacePoint(python_internal_ArrayImpl._get((ds[0] if 0 < len(ds) else None), 0),norm,[u, v],-1,verb_core_Vec.isZero(norm))
#