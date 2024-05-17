from geomdl import CPGen
from geomdl import BSpline
from geomdl import NURBS
import numpy as np
from geomdl import tessellate
from geomdl import operations as op
from geomdl import utilities
from geomdl import exchange
from geomdl.visualization import VisMPL
from geomdl.elements import Vertex, Triangle
import copy
import pickle

# Import and use Matplotlib's colormaps
from matplotlib import cm
import matplotlib.colors as mcolors

colors = [(1,0,0,c) for c in np.linspace(0,1,100)]
cmapred = mcolors.LinearSegmentedColormap.from_list('mycmap', colors, N=5)
colors = [(0,0,1,c) for c in np.linspace(0,1,100)]
cmapblue = mcolors.LinearSegmentedColormap.from_list('mycmap', colors, N=5)


#==============================================Fooling with displayed colourmaps
# # repackage colormap as a more transparent version
# import matplotlib.pylab as pl
# from matplotlib.colors import ListedColormap
#
# # Choose colormap
# #cmap = pl.cm.RdBu
# cmap = pl.cm.cool
#
# # Get the colormap colors
# cmap_t = cmap(np.arange(cmap.N))
#
# # Set alpha
# cmap_t[:,-1] = np.linspace(0, 1, cmap.N)
#
# # Create new colormap
# cmap_t = ListedColormap(cmap_t)

# # random surface generator
# # https://nurbs-python.readthedocs.io/en/latest/surface_generator.html
# # Generate a plane with the dimensions 50x100
# surfgrid = CPGen.Grid(25, 50)
#
# # Generate a grid of 25x30
# surfgrid.generate(30, 50)
#
# # Generate bumps on the grid
# surfgrid.bumps(num_bumps=2, bump_height=10, base_extent=8)
#
# # Create a BSpline surface instance
# surf = BSpline.Surface()
#
# # Set degrees
# surf.degree_u = 3
# surf.degree_v = 3
#
# # Get the control points from the generated grid
# surf.ctrlpts2d = surfgrid.grid
#
# # Set knot vectors
# surf.knotvector_u = utilities.generate_knot_vector(surf.degree_u, surf.ctrlpts_size_u)
# surf.knotvector_v = utilities.generate_knot_vector(surf.degree_v, surf.ctrlpts_size_v)
#
# # Set sample size
# surf.sample_size = 100

verb_EPSILON = 1E-10
machine_EPS = np.finfo(float).eps # float comparison

def radiusCentre3points(p1, p2, p3, tol=10E-14):
    # radius + centre from 3 point circle

    t = p2 - p1
    u = p3 - p1
    v = p3 - p2

    w = np.cross(t, u)  # triangle normal
    wsl = np.linalg.norm(w)
    if wsl < tol: # triangle area too small (additionally check points for colinearity)
        #return False
        return np.inf, np.array([np.inf, np.inf])
    wsl = np.dot(w, w)
    iwsl2 = 1. / (2. * wsl)
    tt = np.dot(t, t)
    uu = np.dot(u, u)

    circCenter = p1 + (u * tt * (np.dot(u, v)) - t * uu * (np.dot(t, v))) * iwsl2
    circRadius = np.sqrt(tt * uu * (np.dot(v, v)) * iwsl2 * 0.5)
    # circAxis   = w / np.sqrt(wsl)

    return circRadius, circCenter

# Lifted from verb, Piegel & Tiller
class AdaptiveRefinementOptions:
    def __init__(self):
        self.minDivsV = 1
        self.minDivsU = 1
        self.refine = True
        self.maxDepth = 10
        self.minDepth = 0
        self.normTol = 2.5e-2

class SurfacePoint: #todo: find geomdl equivalent?
    def __init__(self, point, normal, uv, id = None, degen = None):
        if (id is None):
            id = -1
        if (degen is None):
            degen = False
        self.uv = uv
        self.point = point
        self.normal = normal
        self.id = id
        self.degen = degen

class MeshData:
    def __init__(self, faces, points, normals, uvs):
        self.faces = faces
        self.points = points
        self.normals = normals
        self.uvs = uvs

    def empty():
        return MeshData([], [], [], [])

class AdaptiveRefinementNode:
    # Structure of the child nodes
    # in the adaptive refinement tree
    #
    #   v
    #   ^
    #   |
    #   +--> u
    #
    #                         neighbors[2]
    #
    #                 (u0,v1)---(u05,v1)---(u1,v1)
    #                   |           |          |
    #                   |     3     |     2    |
    #                   |           |          |
    # neighbors[3]   (u0,v05)--(u05,v05)--(u1,v05)   neighbors[1]
    #                   |           |          |
    #                   |     0     |     1    |
    #                   |           |          |
    #                 (u0,v0)---(u05,v0)---(u1,v0)
    #
    #                         neighbors[0]

    def __init__(self,srf,corners,neighbors = None):
        self.v05 = None
        self.u05 = None
        self.horizontal = None
        self.splitHoriz = None
        self.splitVert = None
        self.centerPoint = None
        self.midPoints = None
        self.children = None
        self.srf = srf
        self.neighbors = ([None, None, None, None] if ((neighbors is None)) else neighbors)
        self.corners = corners
        if (self.corners is None):
            u0 = srf.knotsU[0]
            u1 = srf.knotsU[-1]
            v0 = srf.knotsV[0]
            v1 = srf.knotsV[-1]
            self.corners = srf.evaluate_list([[u0, v0], [u1, v0], [u1, v1], [u0, v1]])

    def divide(self,options = None):
        if (options is None):
            options = AdaptiveRefinementOptions()
        if (options.normTol is None):
            options.normTol = 8.5e-2
        if (options.minDepth is None):
            options.minDepth = 0
        if (options.maxDepth is None):
            options.maxDepth = 10
        self._divide(options,0,True)

    def _divide(self, options, currentDepth, horiz):
        self.evalCorners()
        if (not self.shouldDivide(options,currentDepth)):
            return
        currentDepth += 1
        if (self.splitVert and (not self.splitHoriz)):
            horiz = False
        elif ((not self.splitVert) and self.splitHoriz):
            horiz = True
        self.horizontal = horiz

        if (self.horizontal):
            bott = [ self.corners[0], self.corners[1], self.midpoint(1), self.midpoint(3)  ]
            top = [ self.midpoint(3), self.midpoint(1), self.corners[2], self.corners[3]  ]
            self.children = [ AdaptiveRefinementNode( self.srf, bott ), AdaptiveRefinementNode( self.srf, top ) ]

            # assign neighbors to bottom node
            self.children[0].neighbors = [ self.neighbors[0], self.neighbors[1], self.children[1], self.neighbors[3] ]

            # assign neighbors to top node
            self.children[1].neighbors = [ self.children[0], self.neighbors[1], self.neighbors[2], self.neighbors[3] ]

        else:
            left = [ self.corners[0], self.midpoint(0), self.midpoint(2), self.corners[3]  ]
            right = [ self.midpoint(0), self.corners[1], self.corners[2], self.midpoint(2)  ]

            self.children = [ AdaptiveRefinementNode( self.srf, left ), AdaptiveRefinementNode( self.srf, right ) ]

            self.children[0].neighbors = [ self.neighbors[0], self.children[1], self.neighbors[2], self.neighbors[3] ]
            self.children[1].neighbors = [ self.neighbors[0], self.neighbors[1], self.neighbors[2], self.children[0] ]

        # divide all children recursively
        for child in self.children:
            #print("child divide")
            child._divide( options, currentDepth, (not horiz) )

    def isLeaf(self):
        return self.children is None

    def center(self):
        if self.centerPoint is not None:
            return self.centerPoint
        else:
            return self.evalSrf(self.u05, self.v05)

    def evalCorners(self):
        # evaluate center
        self.u05 = (self.corners[0].uv[0] + self.corners[2].uv[0]) / 2
        self.v05 = (self.corners[0].uv[1] + self.corners[2].uv[1]) / 2

        # eval corners
        for i in range(0, 4):
            # if not already evaluated
            if ( self.corners[i].point is None ):
                # evaluate it
                c = self.corners[i]
                self.evalSrf( c.uv[0], c.uv[1], c )

    def evalSrf( self, u, v, srfPt=None ):
        #derivs = Eval.rationalSurfaceDerivatives( this.srf, u, v, 1 );
        #derivs = np.array(self.srf.derivatives(u, v, 1))
        #pt = derivs[0][0]
        #norm = np.cross(  derivs[0][1], derivs[1][0] )
        norm = op.normal(self.srf, [u, v])
        pt = list(norm[0])
        degen = np.isclose( np.array(norm[1]), np.array([0.,0.,0.]), verb_EPSILON ).any()

        # if not degen:
        #     norm = norm/np.linalg.norm(norm)
        if srfPt is not None:
            srfPt.degen = degen
            srfPt.point = pt
            srfPt.normal = np.array(norm[1])#norm
            return srfPt
        else:
            return SurfacePoint( pt, np.array(norm[1]), [u,v], -1, degen )

    def getEdgeCorners(self, edgeIndex ):
        # if leaf, there are no children to obtain uvs from
        if self.isLeaf(): return [ self.corners[ edgeIndex ] ]
        if self.horizontal:
            if edgeIndex == 0: return self.children[0].getEdgeCorners(0)
            elif edgeIndex == 1: return self.children[0].getEdgeCorners(1) + self.children[1].getEdgeCorners(1)
            elif edgeIndex == 2: return self.children[1].getEdgeCorners(2)
            elif edgeIndex == 3: return self.children[1].getEdgeCorners(3) + self.children[0].getEdgeCorners(3)
            else: pass
        #edgeIndex1 = edgeIndex
        if edgeIndex == 0: return self.children[0].getEdgeCorners(0) + self.children[1].getEdgeCorners(0)
        elif edgeIndex == 1: return self.children[1].getEdgeCorners(1)
        elif edgeIndex == 2: return self.children[1].getEdgeCorners(2) + self.children[0].getEdgeCorners(2)
        elif edgeIndex == 3: return self.children[0].getEdgeCorners(3)
        else: pass
        return None

    def getAllCorners( self, edgeIndex ):
        baseArr = [ self.corners[edgeIndex] ]
        if self.neighbors[edgeIndex] is None: return baseArr

        # get opposite edges uvs
        corners = self.neighbors[edgeIndex].getEdgeCorners( ( edgeIndex + 2 ) % 4 )
        funcIndex = edgeIndex % 2
        # e = verb_EPSILON

        # range clipping functions
        def rfm0(c):
            return (c.uv[0] > (self.corners[0].uv[0] + verb_EPSILON)) and (c.uv[0] < (self.corners[2].uv[0] - verb_EPSILON))

        def rfm1(c):
            return (c.uv[1] > (self.corners[0].uv[1] + verb_EPSILON)) and (c.uv[1] < (self.corners[2].uv[1] - verb_EPSILON))

        rangeFuncMap = [rfm0, rfm1]

        # clip the range of uvs to match this one
        cornercopy = list(filter(rangeFuncMap[funcIndex], corners))
        cornercopy.reverse()
        return (baseArr + cornercopy)

    def midpoint(self, index):
        if (self.midPoints is None): self.midPoints = [None, None, None, None]
        if (self.midPoints[index] is not None): return self.midPoints[index]
        index1 = index
        if index1 == 0: self.midPoints[0] = self.evalSrf(self.u05, self.corners[0].uv[1])
        elif index1 == 1: self.midPoints[1] = self.evalSrf(self.corners[1].uv[0], self.v05)
        elif (index1 == 2): self.midPoints[2] = self.evalSrf(self.u05, self.corners[2].uv[1])
        elif (index1 == 3): self.midPoints[3] = self.evalSrf(self.corners[0].uv[0], self.v05)
        else: pass
        return self.midPoints[index]

    def hasBadNormals(self):
        #return self.corners[0].degen or self.corners[1].degen or self.corners[2].degen or self.corners[3].degen
        return any([self.corners[i].degen.any() for i in range(0, 4)])

    def fixNormals(self):
        ll = len(self.corners)
        for i in range(0, ll):
            #corn = self.corners[i]
            if (self.corners[i].degen):
                # get neighbors
                v1 = self.corners[(i + 1) % ll]
                v2 = self.corners[(i + 3) % ll]

                # correct the normal
                if v1.degen: self.corners[i].normal = v2.normal
                else: self.corners[i].normal = v1.normal

    def shouldDivide( self, options, currentDepth):
        if ( currentDepth < options.minDepth ): return True
        if ( currentDepth >= options.maxDepth ): return False
        if ( self.hasBadNormals() ):
            self.fixNormals()
            # don't divide any further when encountering a degenerate normal
            return False

        sv01 = self.corners[0].normal - self.corners[1].normal
        sv23 = self.corners[2].normal - self.corners[3].normal

        self.splitVert = (np.inner(sv01, sv01) > options.normTol) or (np.inner(sv23, sv23) > options.normTol)

        sh12 = self.corners[1].normal - self.corners[2].normal
        sh30 = self.corners[3].normal - self.corners[0].normal

        self.splitHoriz = (np.inner(sh12, sh12) > options.normTol) or (np.inner(sh30, sh30) > options.normTol) # checks out wrt verb

        if self.splitVert or self.splitHoriz: return True

        center = self.center()

        return ((np.inner(center.normal - self.corners[0].normal, center.normal - self.corners[0].normal) > options.normTol) or
                (np.inner(center.normal - self.corners[1].normal, center.normal - self.corners[1].normal) > options.normTol) or
                (np.inner(center.normal - self.corners[2].normal, center.normal - self.corners[2].normal) > options.normTol) or
                (np.inner(center.normal - self.corners[3].normal, center.normal - self.corners[3].normal) > options.normTol))

    def triangulate(self, mesh=None):
        if mesh is None: mesh = MeshData.empty()
        if self.isLeaf(): return self.triangulateLeaf( mesh )

        # recurse on the children
        for c in self.children:
            if c is None: break
            c.triangulate( mesh )
        # [c.triangulate(mesh) for c in self.children if c is not None else break]
        return mesh

    def triangulateLeaf(self, mesh):
        baseIndex = len(mesh.points)
        uvs = []
        ids = []
        splitID = 0
        for i in range(0, 4):
            edgeCorners = self.getAllCorners(i)
            if len(edgeCorners) == 2: splitID += 1
            [uvs.append(ec) for ec in edgeCorners]

        for i in range(0, len(uvs)):
            corner = uvs[i]

            if (corner.id != -1):
                ids.append(corner.id)
                continue

            # geomdl rendering doesn't like numpy arrays
            if isinstance(corner.uv, np.ndarray): mesh.uvs.append(corner.uv.tolist())
            else: mesh.uvs.append(corner.uv)

            if isinstance(corner.point, np.ndarray): mesh.points.append(corner.point.tolist())
            else: mesh.points.append(corner.point)

            #mesh.uvs.append(corner.uv)
            #mesh.points.append(corner.point)
            #mesh.normals.append(corner.normal)
            corner.id = baseIndex
            ids.append(baseIndex)
            baseIndex += 1

        if (len(uvs) == 4):
            mesh.faces.append([ids[0], ids[3], ids[1]])
            mesh.faces.append([ids[3], ids[2], ids[1]])
            return mesh

        elif (len(uvs) == 5):
            mesh.faces.append([ids[splitID], ids[(splitID + 2) % len(ids)], ids[(splitID + 1) % len(ids)]])
            mesh.faces.append([ids[(splitID + 4) % len(ids)], ids[(splitID + 3) % len(ids)], ids[splitID]])
            mesh.faces.append([ids[splitID], ids[(splitID + 3) % len(ids)], ids[(splitID + 2) % len(ids)]])
            return mesh

        # make point at center of face
        center = self.center()

        if isinstance(center.uv, np.ndarray): mesh.uvs.append(center.uv.tolist())
        else: mesh.uvs.append(center.uv)

        if isinstance(center.point, np.ndarray): mesh.points.append(center.point.tolist())
        else: mesh.points.append(center.point)

        #mesh.uvs.append(center.uv)
        #mesh.points.append(center.point)
        #mesh.normals.append(center.normal)
        centerIndex = len(mesh.points) - 1

        # build triangle fan from center
        i = 0
        j = (len(uvs) - 1)
        while (i < len(uvs)):
            mesh.faces.append([centerIndex, ids[i], ids[j]])
            i = (i + 1)
            j = (i - 1)

        return mesh

class Tess:
    def divideRationalSurfaceAdaptive(surface, options=None):
        # https://ariel.chronotext.org/dd/defigueiredo93adaptive.pdf
        def north(index, i, j, divsU, divsV, divs):
            if (i == 0): return None
            return divs[ index - divsU ]

        def south(index, i, j, divsU, divsV, divs):
            if (i == divsV - 1): return None
            return divs[ index + divsU ]

        def east(index, i, j, divsU, divsV, divs):
            if (j == divsU - 1): return None
            return divs[ index + 1 ]

        def west(index, i, j, divsU, divsV, divs):
            if (j == 0): return None
            return divs[ index - 1 ]

        if (options is None): options = AdaptiveRefinementOptions()
        options = AdaptiveRefinementOptions()
        options.minDivsU = options.minDivsU if (options.minDivsU is not None) else 1
        options.minDivsU = options.minDivsV if (options.minDivsV is not None) else 1
        options.refine = options.refine if (options.refine is not None) else True
        minU = (len(surface.ctrlpts2d) - 1) * 2
        minV = (len(surface.ctrlpts2d[0]) - 1) * 2

        divsU = options.minDivsU if options.minDivsU > minU else minU
        options.minDivsU = divsU
        divsV = options.minDivsV if options.minDivsV > minV else minV
        options.minDivsV = divsV

        # get necessary intervals
        umin = surface.knotvector_u[0]
        umax = surface.knotvector_u[-1]
        vmin = surface.knotvector_v[0]
        vmax = surface.knotvector_v[-1]

        du = (umax - umin) / divsU # NURBS division?
        dv = (vmax - vmin) / divsV
        divs = []
        pts = []

        # 1) evaluate all corners
        for i in range(0, divsV + 1):
            ptrow = []
            for j in range(0, divsU + 1):
                u = umin + (du * j)
                v = vmin + (dv * i)
                # todo: make this faster by specifying n,m
                #ds = surface.derivatives(u, v, order=1)
                norm = op.normal(surface, [u, v])
                #norm = (np.cross(ds[0][1], ds[1][0]))/np.linalg.norm(np.cross(ds[0][1], ds[1][0]))
                #ptrow.append(SurfacePoint(ds[0][0], norm, [u, v], -1, np.isclose(norm, np.array([0, 0, 0]), verb_EPSILON).any()))
                ptrow.append(SurfacePoint(np.array(norm[0]), np.array(norm[1]), [u, v], -1, np.isclose(np.array(norm[1]), np.array([0, 0, 0]), verb_EPSILON).any()))
            pts.append(ptrow)

        # 2) make all nodes
        for i in range(0, divsV):
            for j in range(0, divsU):
                corners = [ pts[divsV-i-1][j], pts[divsV-i-1][j+1],	pts[divsV-i][j+1], pts[divsV-i][j] ]
                divs.append(AdaptiveRefinementNode(surface, corners))
        if not options.refine:
            return divs

        # 3) assign all neighbors and divide
        for i in range(0, divsV):
            for j in range(0, divsU):
                ci = i * divsU + j
                n = north( ci, i, j, divsU, divsV, divs )
                e = east( ci, i, j, divsU, divsV, divs  )
                s = south( ci, i, j, divsU, divsV, divs )
                w = west( ci, i, j, divsU, divsV, divs  )
                divs[ci].neighbors = [ s, e, n, w ]
                divs[ci].divide( options )
        return divs

    def rationalSurfaceAdaptive(surface, options=None):
        if (options is None): options = AdaptiveRefinementOptions()
        arrTrees = Tess.divideRationalSurfaceAdaptive(surface, options)
        mesh = MeshData.empty()
        [aT.triangulate(mesh) for aT in arrTrees]
        return mesh

class adaptiveTessellate(tessellate.AbstractTessellate):
    """  abstract tessellation algorithm for surfaces. """
    def __init__(self, **kwargs):
        super(adaptiveTessellate, self).__init__(**kwargs)
        self._tsl_func = adaptiveTransform

    def tessellate(self, points, **kwargs):
        """ Applies tessellation.

        This function does not check if the points have already been tessellated.

        Keyword Arguments:
            * ``size_u``: number of points on the u-direction
            * ``size_v``: number of points on the v-direction

        :param points: array of points
        :type points: list, tuple
        """
        # Call parent function
        super(adaptiveTessellate, self).tessellate(points, **kwargs)

        # Apply  mesh generator function
        self._vertices, self._faces = self._tsl_func(points, **kwargs)

def adaptiveTransform(points, **kwargs):
    """ Generates adaptive mesh from an array of points.

    This function generates a triangular mesh for a NURBS or B-Spline surface on its parametric space.
    The input is the surface points and the number of points on the parametric dimensions u and v,
    indicated as row and column sizes in the function signature. This function should operate correctly if row and
    column sizes are input correctly, no matter what the points are v-ordered or u-ordered. Please see the
    documentation of ``ctrlpts`` and ``ctrlpts2d`` properties of the Surface class for more details on
    point ordering for the surfaces.

    This function accepts the following keyword arguments:

    * ``trims``: List of trim curves passed to the tessellation function
    * ``tessellate_func``: Function called for tessellation. *Default:* :func:`.tessellate.surface_tessellate`
    * ``tessellate_args``: Arguments passed to the tessellation function (as a dict)

    The tessellation function is designed to generate triangles from 4 vertices. It takes 4 :py:class:`.Vertex` objects,
    index values for setting the triangle and vertex IDs and additional parameters as its function arguments.
    It returns a tuple of :py:class:`.Vertex` and :py:class:`.Triangle` object lists generated from the input vertices.
    A default triangle generator is provided as a prototype for implementation in the source code.

    The return value of this function is a tuple containing two lists. First one is the list of vertices and the second
    one is the list of triangles.

    :param points: input points
    :type points: list, tuple

    :return: a tuple containing lists of vertices and triangles
    :rtype: tuple
    """
    # Organization of vertices in a quad element on the parametric space:
    #
    # v4      v3
    # o-------o         i
    # |       |          |
    # |       |          |
    # |       |          |_ _ _
    # o-------o                 j
    # v1      v2

    # def fix_numbering(vertex_list, triangle_list):
    #     # Initialize variables
    #     final_vertices = []
    #
    #     # Get all vertices inside the triangle list
    #     tri_vertex_ids = []
    #     for tri in triangle_list:
    #         for td in tri.data:
    #             if td not in tri_vertex_ids:
    #                 tri_vertex_ids.append(td)
    #
    #     # Find vertices used in triangles
    #     seen_vertices = []
    #     for vertex in vertex_list:
    #         if vertex.id in tri_vertex_ids and vertex.id not in seen_vertices:
    #             final_vertices.append(vertex)
    #             seen_vertices.append(vertex.id)
    #
    #     # Fix vertex numbering (automatically fixes triangle vertex numbering)
    #     vert_new_id = 0
    #     for vertex in final_vertices:
    #         vertex.id = vert_new_id
    #         vert_new_id += 1
    #
    #     return final_vertices, triangle_list

    # # Vertex spacing for triangulation
    trim_curves = kwargs.get('trims', [])

    # Tessellation algorithm
    tsl_func = kwargs.get('tessellate_func')
    if tsl_func is None:
        tsl_func = tessellate.surface_tessellate
    tsl_args = kwargs.get('tessellate_args', dict())

    if type(tsl_args) == MeshData: #AdaptiveRefinementNode:
        _tess = tsl_args
    else:
        print("adaptiveTransform(): adaptive refinement mesh translate data fail")

    # # Numbering
    # vrt_idx = 0  # vertex index numbering start
    # tri_idx = 0  # triangle index numbering start
    #
    # #     # Generate vertices directly from input points (preliminary evaluation)
    # #     vertices = [Vertex() for _ in range(varr_size_v * varr_size_u)]
    # #     u = 0.0
    # #     for i in range(0, size_u, vertex_spacing):
    # #         v = 0.0
    # #         for j in range(0, size_v, vertex_spacing):
    # #             idx = j + (i * size_v)
    # #             vertices[vrt_idx].id = vrt_idx
    # #             vertices[vrt_idx].data = points[idx]
    # #             vertices[vrt_idx].uv = [u, v]
    # #             vrt_idx += 1
    # #             v += v_jump
    # #         u += u_jump
    #
    # vertices = [Vertex() for _ in range(4 * len(ARN))]
    # triangles = []
    # for ar in ARN:
    #     for arc in ar.corners:
    #         if arc is not None:
    #             vertices[vrt_idx].id = arc.id #---------------------------------------
    #             vertices[vrt_idx].data = arc.point
    #             vertices[vrt_idx].uv = arc.uv
    #             vrt_idx += 1
    #
    #     # double node/vertex ID counting problem?
    #     if None in ar.corners:
    #         print("adaptiveTransform(): triangular quadrilateral fail")
    #
    #     vlst, tlst = tsl_func(vertices[vrt_idx-4],
    #                           vertices[vrt_idx-3],
    #                           vertices[vrt_idx-2],
    #                           vertices[vrt_idx-1],
    #                           vrt_idx,
    #                           tri_idx,
    #                           trim_curves,
    #                           tsl_args)
    #
    #     # Add tessellation results to the return lists
    #     vertices += vlst
    #     triangles += tlst
    #
    #     # Increment index values
    #     vrt_idx += len(vlst)
    #     tri_idx += len(tlst)
    #
    # # Fix vertex and triangle numbering (ID values)
    # vertices, triangles = fix_numbering(vertices, triangles)

    vertices = [Vertex() for _ in range(len(_tess.points))]
    for idx in range(len(_tess.points)):
        vertices[idx].id = idx
        vertices[idx].data = _tess.points[idx]
        vertices[idx].uv = _tess.uvs[idx]

    triangles = []
    for idx in range(len(_tess.faces)):
        tri = Triangle()
        tri.id = idx
        tri.add_vertex(vertices[_tess.faces[idx][0]],
                       vertices[_tess.faces[idx][1]],
                       vertices[_tess.faces[idx][2]])
        triangles.append(tri)

    return vertices, triangles

    # if type(tsl_args[0]) == AdaptiveRefinementNode:
    #     ARN = tsl_args
    # else:
    #     print("adaptiveTransform(): adaptive refinement mesh translate data fail")
    #
    # # Numbering
    # vrt_idx = 0  # vertex index numbering start
    # tri_idx = 0  # triangle index numbering start
    #
    # vertices = [Vertex() for _ in range(4 * len(ARN))]
    # triangles = []
    # for ar in ARN:
    #     for arc in ar.corners:
    #         if arc is not None:
    #             vertices[vrt_idx].id = arc.id # idx
    #             vertices[vrt_idx].data = arc.point #-------------------------------------------------------------------------
    #             vertices[vrt_idx].uv = arc.uv
    #             vrt_idx += 1
    #
    #     # double node/vertex ID counting problem?
    #     if None in ar.corners:
    #         print("adaptiveTransform(): triangular quadrilateral fail")
    #
    #     vlst, tlst = tsl_func(vertices[vrt_idx-4],
    #                           vertices[vrt_idx-3],
    #                           vertices[vrt_idx-2],
    #                           vertices[vrt_idx-1],
    #                           vrt_idx,
    #                           tri_idx,
    #                           trim_curves,
    #                           tsl_args)
    #
    #     # Add tessellation results to the return lists
    #     vertices += vlst
    #     triangles += tlst
    #
    #     # Increment index values
    #     vrt_idx += len(vlst)
    #     tri_idx += len(tlst)
    #
    # # Fix vertex and triangle numbering (ID values)
    # vertices, triangles = fix_numbering(vertices, triangles)
    #
    # return vertices, triangles

def rationalSurfaceClosestParam(surface, p):
    # for surfaces, minimize the following:
    #
    # f = Su(u,v) * r = 0
    # g = Sv(u,v) * r = 0
    #
    #   where r = S(u,v) - P
    #
    # Requires Newton-Raphson iteration, objective function is vector valued
    #
    #     J d = k
    #
    #     d =   [ u* - u, v* - v ] (alternatively u[i+1] -  u[i], v[i+1]- v[i])
    #     k = - [ f(u,v), g(u,v) ]
    #     J =     |Su|^2   +  Suu * r       Su*Sv  +  Suv * r
    #              Su*Sv   +  Svu * r      |Sv|^2  +  Svv * r
    #
    # halting conditions:
    #
    # point coincidence
    #
    #         |S(u,v) - p| < e1
    #
    # cosine
    #
    #    |Su(u,v)*(S(u,v) - P)|
    #    ----------------------  < e2
    #    |Su(u,v)| |S(u,v) - P|
    #
    #    |Sv(u,v)*(S(u,v) - P)|
    #    ----------------------  < e2
    #    |Sv(u,v)| |S(u,v) - P|
    #
    # 1) first check 2 & 3
    # 2) if at least one of these is not, compute new value, otherwise halt
    # 3) ensure the parameter stays within range
    #     * if not closed, don't allow outside of range a-b
    #     * if closed (e.g. circle), allow to move back to beginning
    # 4)  if |(u* - u)C'(u)| < e1, halt

    eps1 = 0.0001
    eps2 = 0.0005
    dif = None

    minu = surface.knotvector_u[0]
    maxu = surface.knotvector_u[-1]
    minv = surface.knotvector_v[0]
    maxv = surface.knotvector_v[-1]

    # check if surface closed along U and V
    closedu = np.isclose(np.array(surface.ctrlpts2d[0]), np.array(surface.ctrlpts2d[-1]), verb_EPSILON).all()  # eps_STEP_AP21
    closedv = np.isclose(np.array(surface.ctrlpts2d[0]).T, np.array(surface.ctrlpts2d[-1]).T, verb_EPSILON).all()  # eps_STEP_AP21
    cuv = None

    # todo: divide surface instead of a full on tessellation

    # approximate closest point with tessellation
    tess = Tess.rationalSurfaceAdaptive(surface, AdaptiveRefinementOptions())
    #tess = divideRationalSurfaceAdaptive(surface, AdaptiveRefinementOptions())
    #dmin = Math.POSITIVE_INFINITY
    dmin = np.inf

    for i in range(0, len(tess.points)):
        x = np.array(tess.points[i])
        d = np.inner(p - x, p - x)
        if (d < dmin):#-------------------------------<<<<<<
            dmin = d
            cuv = tess.uvs[i]

    # def f(uv):
    #     return surface.derivatives(uv[0], uv[1], 2)

    def n(uv, e, r): # np.array inputs
        #   f = Su(u,v) * r = 0
        #   g = Sv(u,v) * r = 0

        Su = e[1][0]
        Sv = e[0][1]

        Suu = e[2][0]
        Svv = e[0][2]

        Suv = e[1][1]
        Svu = e[1][1]

        f = np.dot(Su, r)
        g = np.dot(Sv, r)
        k = [-f, -g]

        J00 = np.dot( Su, Su ) + np.dot( Suu, r )
        J01 = np.dot( Su, Sv ) + np.dot( Suv, r )
        J10 = np.dot( Su, Sv ) + np.dot( Svu, r )
        J11 = np.dot( Sv, Sv ) + np.dot( Svv, r )

        # d =   [ u* - u, v* - v ]
        # k = - [ f(u,v), g(u,v) ]
        # J =   |Su|^2   +  Suu * r       Su*Sv  +  Suv * r
        #        Su*Sv   +  Svu * r      |Sv|^2  +  Svv * r

        J = [[J00, J01], [J10, J11]]
        #d = verb_core_Mat.solve(J,k)
        d = np.linalg.solve(J,k)
        return d + uv

    maxits = 5
    i = 0
    #e = None
    while (i < maxits):
        e = np.array(surface.derivatives(cuv[0], cuv[1], 2))#--------------------------------------------------------------------------------------------------------------------
        dif = e[0][0] - p

        #   point coincidence   |S(u,v) - p| < e1

        c1v = np.linalg.norm(dif)

          # cosine
          #
          #  |Su(u,v)*(S(u,v) - P)|
          #  ----------------------  < e2
          #  |Su(u,v)| |S(u,v) - P|
          #
          #  |Sv(u,v)*(S(u,v) - P)|
          #  ----------------------  < e2
          #  |Sv(u,v)| |S(u,v) - P|

        c2an = np.dot( e[1][0], dif)
        c2ad = np.linalg.norm( e[1][0] ) * c1v
        c2bn = np.dot( e[0][1], dif)
        c2bd = np.linalg.norm( e[0][1] ) * c1v
        c2av = c2an / c2ad
        c2bv = c2bn / c2bd
        c1 = c1v < eps1
        c2a = c2av < eps2
        c2b = c2bv < eps2

        # exit if all tolerance are met,
        if (c1 and c2a and c2b):
            #pass
            return cuv

        # otherwise, take a step
        ct = n(cuv, e, dif)

        #  correct for exceeding bounds
        if ct[0] < minu:
            if closedu:
                ct = [maxu - (ct[0] - minu), ct[1]]
            else:
                ct = [minu + verb_EPSILON, ct[1]]
        elif ct[0] > maxu:
            if closedu:
                ct = [minu + (ct[0] - maxu), ct[1]]
            else:
                ct = [maxu - verb_EPSILON, ct[1]]

        if ct[1] < minv:
            if closedv:
                ct = [ct[0], maxv - ( ct[1] - minv )]
            else:
                ct = [ct[0], minv + verb_EPSILON]
        elif ct[1] > maxv:
            if closedv:
                ct = [ct[0], minv + ( ct[0] - maxv )]
            else:
                ct = [ct[0], maxv - verb_EPSILON]

        c3v0 =  np.linalg.norm( (ct[0] - cuv[0]) * e[1][0] )
        c3v1 =  np.linalg.norm( (ct[1] - cuv[1]) * e[0][1] )

        # if |(u* - u) C'(u)| < e1, halt
        if (c3v0 + c3v1 < eps1):
            return cuv
            #pass

        cuv = ct
        i = (i + 1)
    return cuv
    #pass

def rationalSurfaceExtremaParam(surface, p, maxSearch=True, localExtrema=False, curvatureTest=False, eps1=0.0001, eps2=0.0005):
    # easy to find a minimum, but a maximum found through the original hillclimbing algorithm will settle on local maxima
    # for surfaces, minimize the following:
    #
    # f = Su(u,v) * r = 0
    # g = Sv(u,v) * r = 0
    #
    #   where r = S(u,v) - P
    #
    # Requires Newton-Raphson iteration, objective function is vector valued
    #
    #     J d = k
    #
    #     d =   [ u* - u, v* - v ] (alternatively u[i+1] -  u[i], v[i+1]- v[i])
    #     k = - [ f(u,v), g(u,v) ]
    #     J =     |Su|^2   +  Suu * r       Su*Sv  +  Suv * r
    #              Su*Sv   +  Svu * r      |Sv|^2  +  Svv * r
    #
    # halting conditions:
    #
    # point coincidence
    #
    #         |S(u,v) - p| < e1
    #
    # cosine
    #
    #    |Su(u,v)*(S(u,v) - P)|
    #    ----------------------  < e2
    #    |Su(u,v)| |S(u,v) - P|
    #
    #    |Sv(u,v)*(S(u,v) - P)|
    #    ----------------------  < e2
    #    |Sv(u,v)| |S(u,v) - P|
    #
    # 1) first check 2 & 3
    # 2) if at least one of these is not, compute new value, otherwise halt
    # 3) ensure the parameter stays within range
    #     * if not closed, don't allow outside of range a-b
    #     * if closed (e.g. circle), allow to move back to beginning
    # 4)  if |(u* - u)C'(u)| < e1, halt

    # eps1 = 0.0001
    # eps2 = 0.0005

    minu = surface.knotvector_u[0]
    maxu = surface.knotvector_u[-1]
    minv = surface.knotvector_v[0]
    maxv = surface.knotvector_v[-1]

    # check if surface closed along U and V
    closedu = np.isclose(np.array(surface.ctrlpts2d[0]), np.array(surface.ctrlpts2d[-1]), verb_EPSILON).all()  # eps_STEP_AP21
    closedv = np.isclose(np.array(surface.ctrlpts2d[0]).T, np.array(surface.ctrlpts2d[-1]).T, verb_EPSILON).all()  # eps_STEP_AP21

    # # verb adaptive tesselation seems pretty slow, code below not much faster
    # surface.tessellator = tessellate.TriangularTessellate() #QuadTessellate() broken--------------------trims argument unrecognised
    # surf.tessellate()
    #
    # surroundNodes = dict()
    # for stv in surface.vertices:
    #     # find faces containing vertex [x for b in a for x in b]
    #     adjTri = [stf.id for stf in surface.faces if stv.id in stf.vertex_ids]
    #     adjV = [p for atr in adjTri for p in surface.faces[atr].data]
    #     adjV = list(set(adjV))
    #     adjV.remove(stv.id)
    #     uhi = None
    #     ulo = None
    #     vhi = None
    #     vlo = None
    #
    #     for av in adjV:
    #         if np.abs(surface.vertices[av].u - stv.u) < (surface.delta_u/2): # N or S
    #             if surface.vertices[av].v - stv.v > (surface.delta_v/2):
    #                 uhi = av
    #             elif surface.vertices[av].v - stv.v < -(surface.delta_v/2):
    #                 ulo = av
    #         if np.abs(surface.vertices[av].v - stv.v) < (surface.delta_v / 2):  # E or W
    #             if surface.vertices[av].u - stv.u > (surface.delta_u / 2):
    #                 vhi = av
    #             elif surface.vertices[av].u - stv.u < -(surface.delta_u/2):
    #                 vlo = av
    #     surroundNodes[stv.id] = [ulo, uhi, vlo, vhi]


    # todo: divide surface instead of a full on tessellation, probably means quicksort?
    # approximate closest point with tessellation
    # verb adaptive tesselation
    tess = Tess.rationalSurfaceAdaptive(surface)

    # convert this tessellation triangulation to geomdl format
    surface.tessellator = adaptiveTessellate()
    surface.tessellate(tessellate_args=tess, force=True)

    # testing the minimum UV displacement value.
    #     U = [t[0] for t in tess.uvs]
    #     U = sorted(list(set(U)))
    #     U = np.diff(U).min()
    #     V = [t[1] for t in tess.uvs]
    #     V = sorted(list(set(V)))
    #     V = np.diff(V).min()
    #     UVdiffMin = min[U, V])

    nodeTol = min([surface.delta_u, surface.delta_v])/2

    # finding surrounding nodes in order to get curvature along u, v axes is complicated by adaptive mesh tesselation
    surroundNodes = dict()
    uvMax = [uvi for uvi, uv in enumerate(tess.uvs) if np.isclose(uv, np.array([maxu, maxv]), nodeTol).all()][0]
    uvMin = [uvi for uvi, uv in enumerate(tess.uvs) if np.isclose(uv, np.array([minu, minv]), nodeTol).all()][0]

    # for every u,v, get a set of surrounding u, v values
    for tuv in range(0, len(tess.uvs)):
        if not (np.isclose(tess.uvs[tuv][0], maxu, nodeTol) or np.isclose(tess.uvs[tuv][0], minu, nodeTol) or
                np.isclose(tess.uvs[tuv][1], maxv, nodeTol) or np.isclose(tess.uvs[tuv][1], minv, nodeTol)): # unlikely to work for closed surfaces
            ulo = uvMin
            uhi = uvMax
            vlo = uvMin
            vhi = uvMax
            for tuv2 in range(0, len(tess.uvs)):
                if (tuv2 != tuv):
                    # get nearest node along exact u, v values
                    if np.isclose(tess.uvs[tuv][0], tess.uvs[tuv2][0], nodeTol): # nodes along u value
                        if (tess.uvs[tuv2][1] < tess.uvs[tuv][1]) and (tess.uvs[tuv2][1] > tess.uvs[vlo][1]):
                            vlo = tuv2
                        if (tess.uvs[tuv2][1] > tess.uvs[tuv][1]) and (tess.uvs[tuv2][1] < tess.uvs[vhi][1]):
                            vhi = tuv2

                    if np.isclose(tess.uvs[tuv][1], tess.uvs[tuv2][1], nodeTol):  # nodes along v value
                        if (tess.uvs[tuv2][0] < tess.uvs[tuv][0]) and (tess.uvs[tuv2][0] > tess.uvs[ulo][0]):
                           ulo = tuv2
                        if (tess.uvs[tuv2][0] > tess.uvs[tuv][0]) and (tess.uvs[tuv2][0] < tess.uvs[uhi][0]):
                            uhi = tuv2

            surroundNodes[tuv] = [ulo, uhi, vlo, vhi]

    # create two points on either side of each point along U and V direction (obtain point equivalents)
    # use side points and point under scrutiny to determine centre of a circle passing through these 3 points

    d_min = np.inf
    localUVextremaDisp = []
    localExtremaSet = []

    for tuv in range(0, len(tess.uvs)):
        if tuv in surroundNodes:
            # [ulo, uhi, vlo, vhi]
            pulo = np.array(tess.points[surroundNodes[tuv][0]])
            puhi = np.array(tess.points[surroundNodes[tuv][1]])
            pvlo = np.array(tess.points[surroundNodes[tuv][2]])
            pvhi = np.array(tess.points[surroundNodes[tuv][3]])

            #simple test
            if not curvatureTest:
                dpuv = np.linalg.norm(p - np.array(tess.points[tuv]))
                dppulo = np.linalg.norm(p - pulo)
                dppuhi = np.linalg.norm(p - puhi)
                dppvlo = np.linalg.norm(p - pvlo)
                dppvhi = np.linalg.norm(p - pvhi)

                if maxSearch:
                    if (dpuv > dppulo) and (dpuv > dppuhi) and (dpuv > dppvlo) and (dpuv > dppvhi):
                        #localUVextremaDisp.append(tess.uvs[tuv])
                        localExtremaSet.append(tuv)
                else: #minSearch
                    if (dpuv < dppulo) and (dpuv < dppuhi) and (dpuv < dppvlo) and (dpuv < dppvhi):
                        #localUVextremaDisp.append(tess.uvs[tuv])
                        localExtremaSet.append(tuv)

            else:
                # complex test summing curvatures along U and V axes
                localCurvature = [np.inf] * len(tess.uvs)
                pur, pucc = radiusCentre3points(np.array(pulo), np.array(tess.points[tuv]), np.array(puhi))
                pvr, pvcc = radiusCentre3points(np.array(pvlo), np.array(tess.points[tuv]), np.array(pvhi))

                dpuv = np.linalg.norm(p - np.array(tess.points[tuv]))
                dpucc = np.linalg.norm(p - pucc)  # distance circle centre
                dpvcc = np.linalg.norm(p - pvcc)

                if (pur is not np.inf) and (pvr is not np.inf):  # no planar or ruled surface can have local minima within edges

                    # overall sphericity at point tuv is estimated from averaged distance of U-radius and V-radius from p
                    duvk = (dpuv - dpvcc) + (dpuv - dpvcc)
                    dpuvcc = (dpucc + dpvcc)/2 - dpuv

                    if maxSearch:
                        # spline version used a maxima of 3 points, but as outlying points are averages, this seems excessive
                        # dua = max([np.linalg.norm(p - pu) for pu in [puLower, tess.points[tuv][0], puHigher]])
                        # dva = max([np.linalg.norm(p - pv) for pv in [pvLower, tess.points[tuv][1], pvHigher]])

                        # centre of 3 point circle is nearer or further than max/min sample point from p
                        #orthoSign = dua > (dcc_p + (2*dcc)) # circle centre closer to p than circle points => maximum curve

                        # note that there is a difference between the unique minima/maxima point of any surface and a
                        # collection of dimples or bumps that identify a surface
                        localCurvature[tuv] = np.abs(duvk)
                        if (dpuvcc < dpuv) and (np.abs(duvk) < d_min): # maxima, > dpuv for minima
                            d_min = np.abs(duvk)
                            # instead rank curvature?
                            localUVextremaDisp.append(tess.uvs[tuv])
                            #localCurvature[tuv] = np.abs(duvk)

                    else: #minSearch

                        if (dpuvcc > dpuv) and (np.abs(duvk) < d_min): # minima
                            d_min = np.abs(duvk)
                            localUVextremaDisp.append(tess.uvs[tuv])
                            #localCurvature[tuv] = np.abs(duvk)

    if curvatureTest:
        for sn in surroundNodes:
            if not any(localCurvature[ssn] == np.inf for ssn in surroundNodes[sn]):
                if all(localCurvature[ssn] > localCurvature[sn] for ssn in surroundNodes[sn]):
                    localExtremaSet.append(sn)

                    # identify curvature value local minima #np.argsort(pDisps)[:window2]

                    # for each face within tesselation (arrTree) find face point position on tessa ps, point normal pn, and set of surrounding points (sps, spn)
                    # get projected point normal ppn = |pp - ps|, compare with surface point normal, np.dot(ppn, spn) for minimum value
                    # worthwhile doing exhaustive search?
                    # u||v precisely when (u⋅v)2 = (u⋅u)(v⋅v)
                    # doughnut requirement for set of local surface maxima: furthest point must be surrounded by less-far surface points

                    # import matplotlib.pyplot as plt
                    # fig = plt.figure()
                    # ax = fig.gca(projection='3d')
                    # for i in range(0, len(localCurvature)):
                    #     if localCurvature[i] == np.inf:
                    #         label = "inf"
                    #     else:
                    #         label = '%d' % (localCurvature[i])
                    #     ax.text(tess.points[i][0], tess.points[i][1], tess.points[i][2], label)
                    # # # Tweaking display region and labels
                    # # ax.set_xlim(-30, 30)
                    # # ax.set_ylim(-30, 30)
                    # # ax.set_zlim(-10, 10)
                    # # ax.set_xlabel('X axis')
                    # # ax.set_ylabel('Y axis')
                    # # ax.set_zlabel('Z axis')
                    # plt.show()



    # localSurfaceMax = []
    # localSurfaceMaxDisp = []
    # localUVmaxDisp = []
    # localSurfacePointMax = []
    # for tp in range(0, len(tess.points)):
    #     # strip out u, v values at the edge of non-closed surfaces
    #     if not(np.isclose(tess.uvs[tp][0], minu) or np.isclose(tess.uvs[tp][0], maxu) or np.isclose(tess.uvs[tp][1], minv) or np.isclose(tess.uvs[tp][1], maxv)):
    #
    #         PS = p - np.array(tess.points[tp])
    #         dispPS = np.linalg.norm(PS)
    #
    #         # surroundFaces= [tf for tf in tess.faces if tp in tf]
    #         # surroundNodes = []
    #         # [surroundNodes.append(sn) for sf in surroundFaces for sn in sf if sn != tp and sn not in surroundNodes]
    #         # dispSurroundPoints = [np.linalg.norm(p - np.array(tess.points[sn])) for sn in surroundNodes]
    #
    #         # # run a test, whether all tess.faces mentioned in connection are actually the closest points by UV
    #         # print("----------------------------")
    #         # print(surroundNodes)
    #         # uvcheck = [np.linalg.norm(np.array(tess.uvs[tp]) - np.array(tess.uvs[uv])) for uv in range(0, len(tess.points))]
    #         # nearUV = np.argsort(uvcheck)[1:10]
    #         # print(nearUV)
    #         # print([uvcheck[n] for n in nearUV])
    #
    #         dispSurroundPoints = [np.linalg.norm(p - np.array(tess.points[sn])) for sn in surroundNodes[tp]]
    #         #localMaxDisp = all((dispPS + machine_EPS) > dsp for dsp in dispSurroundPoints)
    #         localMaxDisp = all((dsp < dispPS and not np.isclose(dsp, dispPS)) for dsp in dispSurroundPoints)
    #         #print([(dsp < dispPS and not np.isclose(dsp, dispPS)) for dsp in dispSurroundPoints])
    #         #print([dispPS] + [dsp for dsp in dispSurroundPoints])
    #
    #         if localMaxDisp and tess.uvs[tp] not in localUVmaxDisp:
    #             localSurfaceMax.append(tp)
    #             localSurfaceMaxDisp.append(dispPS)
    #             localUVmaxDisp.append(tess.uvs[tp])
    #             localSurfacePointMax.append(tess.points[tp])

                # normDispPS = (PS / dispPS)
                #surroundUVS = [tess.uvs[sn] for sn in surroundNodes]
                #surroundNorms = [tess.normals[sn] for sn in surroundNodes]

                # np.dot(tess.normals[tp], normDispPS)
                #
                # #nearParallel = np.cross(normDispPS, tess.normals[tp])
                # localMaxNormal =  all(np.dot(tess.normals[tp], normDispPS) >= np.dot(normDispPS, sn) for sn in surroundNorms)
                # # if the normal at this point is not more parallel to the projection vector than surrounding point normals..
                # if not localMaxNormal:
                #     print("-------------------------------")
                #     print(tp)
                #     print(surroundNodes)
                #     print(np.dot(tess.normals[tp], normDispPS))
                #     print([np.dot(normDispPS, sn) for sn in surroundNorms])
                #     print(tess.points[tp])
                #     pass

    #dmin = np.inf
    # dmin = 0
    # for i in range(0, len(tess.points)):
    #     x = np.array(tess.points[i])
    #     d = np.inner(p - x, p - x)
    #     if (d > dmin): # d > dmin for maxima?
    #         dmin = d
    #         _cuv = tess.uvs[i]
    # print(_cuv)

    # def f(uv):
    #     return surface.derivatives(uv[0], uv[1], 2)

    def n(uv, e, r): # np.array inputs
        #   f = Su(u,v) * r = 0
        #   g = Sv(u,v) * r = 0

        Su = e[1][0]
        Sv = e[0][1]

        Suu = e[2][0]
        Svv = e[0][2]

        Suv = e[1][1]
        Svu = e[1][1]

        f = np.dot(Su, r)
        g = np.dot(Sv, r)
        k = [-f, -g]

        J00 = np.dot( Su, Su ) + np.dot( Suu, r )
        J01 = np.dot( Su, Sv ) + np.dot( Suv, r )
        J10 = np.dot( Su, Sv ) + np.dot( Svu, r )
        J11 = np.dot( Sv, Sv ) + np.dot( Svv, r )

        # d =   [ u* - u, v* - v ]
        # k = - [ f(u,v), g(u,v) ]
        # J =   |Su|^2   +  Suu * r       Su*Sv  +  Suv * r
        #        Su*Sv   +  Svu * r      |Sv|^2  +  Svv * r

        J = [[J00, J01], [J10, J11]]
        #d = verb_core_Mat.solve(J,k)
        d = np.linalg.solve(J,k)
        return d + uv

    def testConditions(cuv, maxits=5):
        i = 0
        #maxits = 5
        #e = None
        while i < maxits:
            e = np.array(surface.derivatives(cuv[0], cuv[1], 2))#--------------------------------------------------------------------------------------------------------------------
            dif = e[0][0] - p

            #   point coincidence   |S(u,v) - p| < e1
            c1v = np.linalg.norm(dif)
              # cosine
              #
              #  |Su(u,v)*(S(u,v) - P)|
              #  ----------------------  < e2
              #  |Su(u,v)| |S(u,v) - P|
              #
              #  |Sv(u,v)*(S(u,v) - P)|
              #  ----------------------  < e2
              #  |Sv(u,v)| |S(u,v) - P|

            c2an = np.dot( e[1][0], dif)
            c2ad = np.linalg.norm( e[1][0] ) * c1v
            c2bn = np.dot( e[0][1], dif)
            c2bd = np.linalg.norm( e[0][1] ) * c1v
            c2av = c2an / c2ad
            c2bv = c2bn / c2bd
            c1 = c1v < eps1
            c2a = c2av < eps2
            c2b = c2bv < eps2

            # exit if all tolerance are met,
            if (c1 and c2a and c2b):
                return cuv

            # otherwise, take a step
            ct = n(cuv, e, dif)

            #  correct for exceeding bounds
            if ct[0] < minu:
                if closedu: ct = [maxu - (ct[0] - minu), ct[1]]
                else: ct = [minu + verb_EPSILON, ct[1]]
            elif ct[0] > maxu:
                if closedu: ct = [minu + (ct[0] - maxu), ct[1]]
                else: ct = [maxu - verb_EPSILON, ct[1]]

            if ct[1] < minv:
                if closedv: ct = [ct[0], maxv - ( ct[1] - minv )]
                else: ct = [ct[0], minv + verb_EPSILON]
            elif ct[1] > maxv:
                if closedv: ct = [ct[0], minv + ( ct[0] - maxv )]
                else: ct = [ct[0], maxv - verb_EPSILON]

            c3v0 =  np.linalg.norm( (ct[0] - cuv[0]) * e[1][0] )
            c3v1 =  np.linalg.norm( (ct[1] - cuv[1]) * e[0][1] )

            # if |(u* - u) C'(u)| < e1, halt
            if (c3v0 + c3v1 < eps1):
                return cuv
                #break
            cuv = ct
            i = (i + 1)
        return cuv

    if localExtrema: # get the set of local extrema rather than global min/max of surface
        maximaUV = []
        for luv in localExtremaSet:
            maximaUV.append(testConditions(tess.uvs[luv]))

        # filter out identical values
        localExtremaUV = []
        for m in maximaUV:
            # displacement between any 2 points < eps1
            # compare u, v individually to reduce comparison pool
            if not any(np.isclose(m, u, eps1).all() for u in localExtremaUV):
                localExtremaUV.append(m)

        return localExtremaUV #, [tess.uvs[le] for le in localExtremaSet]
    else:
        dpuv = [np.linalg.norm(p - np.array(tess.points[luv])) for luv in localExtremaSet]
        return testConditions(tess.uvs[localExtremaSet[dpuv.index(max(dpuv))]])


def rationalSurfaceExtremaParam_2(S, p, maxS=True, localExtrema=False, curvatureTest=False, uv_xyz=True, eps1=0.0001, eps2=0.0005, delta=0.025):
    '''
    localExtrema: find any local extrema which differs from surrounding sample points, tolerence not sufficiently defined
    uv_xyz: return u,v normalised values or cartesian x,y,z values
    delta: interval dividing surface into points
    '''

    # easy to find a minimum, but a maximum found through the original hillclimbing algorithm will settle on local maxima
    # for surfaces, minimize the following:
    #
    # f = Su(u,v) * r = 0
    # g = Sv(u,v) * r = 0
    #
    #   where r = S(u,v) - P
    #
    # Requires Newton-Raphson iteration, objective function is vector valued
    #
    #     J d = k
    #
    #     d =   [ u* - u, v* - v ] (alternatively u[i+1] -  u[i], v[i+1]- v[i])
    #     k = - [ f(u,v), g(u,v) ]
    #     J =     |Su|^2   +  Suu * r       Su*Sv  +  Suv * r
    #              Su*Sv   +  Svu * r      |Sv|^2  +  Svv * r
    #
    # halting conditions:
    #
    # point coincidence
    #
    #         |S(u,v) - p| < e1
    #
    # cosine
    #
    #    |Su(u,v)*(S(u,v) - P)|
    #    ----------------------  < e2
    #    |Su(u,v)| |S(u,v) - P|
    #
    #    |Sv(u,v)*(S(u,v) - P)|
    #    ----------------------  < e2
    #    |Sv(u,v)| |S(u,v) - P|
    #
    # 1) first check 2 & 3
    # 2) if at least one of these is not, compute new value, otherwise halt
    # 3) ensure the parameter stays within range
    #     * if not closed, don't allow outside of range a-b
    #     * if closed (e.g. circle), allow to move back to beginning
    # 4)  if |(u* - u)C'(u)| < e1, halt

    # eps1 = 0.0001
    # eps2 = 0.0005

    # minu = S.knotvector_u[0]
    # maxu = S.knotvector_u[-1]
    # minv = S.knotvector_v[0]
    # maxv = S.knotvector_v[-1]

    # todo: take value from STEP fields
    # check if surface closed along U and V
    #closedU = np.isclose(np.array(S.ctrlpts2d[0]), np.array(S.ctrlpts2d[-1]), verb_EPSILON).all()  # eps_STEP_AP21
    #def closedU(): return np.isclose(np.array(S.ctrlpts2d[0]), np.array(S.ctrlpts2d[-1]), verb_EPSILON).all()  # eps_STEP_AP21
    #closedV = np.isclose(np.array(S.ctrlpts2d[0]).T, np.array(S.ctrlpts2d[-1]).T, verb_EPSILON).all()  # eps_STEP_AP21
    def closedV(): return np.isclose(np.array(S.ctrlpts2d[0]).T, np.array(S.ctrlpts2d[-1]).T, verb_EPSILON).all()  # eps_STEP_AP21

    S.delta = delta # set evaluation delta
    # tradeoff between sample density and feature identification

    if S.evalpts == None: # evaluate surface points---------appears redundant
        S.evaluate()

    evalptsU = S.sample_size_u
    evalptsV = S.sample_size_v

    # # find maxima/minima from surf.evalpts
    # evalptsU = int(np.round(1 / S.delta_u))
    # evalptsV = int(np.round(1 / S.delta_u))
    # if not evalptsU * evalptsV == len(surf.evalpts):
    #     print("rationalSurfaceExtremaParam_2(): unexpected evalpts division")

    d_min = np.inf
    localExtremaUV = []

    # def globalNeighbourSearch():
    #     # remove values at edges of surface
    #     edgeIndex = ([r for r in range(0, S.sample_size_v)] +
    #      [r for r in range(0, S.sample_size_v * S.sample_size_u) if r % S.sample_size_v == 0] +
    #      [r for r in range(0, S.sample_size_v * S.sample_size_u) if r % S.sample_size_v == S.sample_size_v - 1] +
    #      [r for r in range(S.sample_size_v * S.sample_size_u - S.sample_size_v, S.sample_size_v * S.sample_size_u)])
    #
    #     edgeIndex = list(set(edgeIndex))
    #
    #     edgeFreeIndex = [efi for efi in range(0, S.sample_size_v * S.sample_size_u) if efi not in edgeIndex]
    #
    #     extremaUnique = [np.linalg.norm(p - S.evalpts[efi]) for efi in edgeFreeIndex]
    #
    #     if maxS:
    #         extremaUnique = edgeFreeIndex[extremaUnique.index(max(extremaUnique))]
    #     else:
    #         extremaUnique = edgeFreeIndex[extremaUnique.index(min(extremaUnique))]
    #
    #     return (extremaUnique // S.sample_size_v, extremaUnique % S.sample_size_v)

    if not localExtrema and not curvatureTest:
        # remove values at edges of surface
        edgeIndex = ([r for r in range(0, evalptsV)] +
         [r for r in range(0, evalptsV * evalptsU) if r % evalptsV == 0] +
         [r for r in range(0, evalptsV*evalptsU) if r % evalptsV == evalptsV - 1] +
         [r for r in range(evalptsV*evalptsU - evalptsV, evalptsV*evalptsU)])
        edgeIndex = list(set(edgeIndex))
        edgeFreeIndex = [efi for efi in range(0, evalptsV*evalptsU) if efi not in edgeIndex]
        localExtremaUnique = [np.linalg.norm(p - surf.evalpts[efi]) for efi in edgeFreeIndex]

        if maxS:
            localExtremaUnique = edgeFreeIndex[localExtremaUnique.index(max(localExtremaUnique))]
        else:
            localExtremaUnique = edgeFreeIndex[localExtremaUnique.index(min(localExtremaUnique))]

        localExtremaUnique = (localExtremaUnique // evalptsV, localExtremaUnique % evalptsV)

    # def localNeighbourSearch():
    #     extremaUV = []
    #     for u in range(1, S.sample_size_u - 1):
    #         for v in range(1, S.sample_size_v - 1):
    #
    #             # print([ S.sample_size_u * (u + 1) + (v - 1), S.sample_size_u * (u + 1) + v, S.sample_size_u * (u + 1) + (v + 1)])
    #             # print([S.sample_size_u * u + (v - 1), S.sample_size_u * u + v, S.sample_size_u * u + (v + 1)])
    #             # print([ S.sample_size_u * (u - 1) + (v - 1), S.sample_size_u * (u - 1) + v, S.sample_size_u * (u - 1) + (v + 1)])
    #
    #             if v > 1:
    #                 dpuv_NW = dpuv_N
    #                 dpuv_W = dpuv
    #                 dpuv_SW = dpuv_S
    #
    #                 dpuv_N = dpuv_NE
    #                 dpuv = dpuv_E
    #                 dpuv_S = dpuv_SE
    #             else:
    #                 dpuv_NW = np.linalg.norm(p - S.evalpts[S.sample_size_u * (u + 1)])  # + (v - 1)])
    #                 dpuv_W = np.linalg.norm(p - S.evalpts[S.sample_size_u * u])  # + (v - 1)])
    #                 dpuv_SW = np.linalg.norm(p - S.evalpts[S.sample_size_u * (u - 1)])  # + (v - 1)])
    #
    #                 dpuv_N = np.linalg.norm(p - S.evalpts[(S.sample_size_u * (u + 1)) + 1])  # v])
    #                 dpuv = np.linalg.norm(p - S.evalpts[(S.sample_size_u * u) + 1])  # v])
    #                 dpuv_S = np.linalg.norm(p - S.evalpts[(S.sample_size_u * (u - 1)) + 1])  # v])
    #
    #             # only have to calculate this once
    #             dpuv_NE = np.linalg.norm(p - S.evalpts[(S.sample_size_u * (u + 1)) + v + 1])
    #             dpuv_E = np.linalg.norm(p - S.evalpts[(S.sample_size_u * u) + v + 1])
    #             dpuv_SE = np.linalg.norm(p - S.evalpts[(S.sample_size_u * (u - 1)) + v + 1])
    #
    #             # print([dpuv_NW, dpuv_N, dpuv_NE])
    #             # print([dpuv_W , dpuv ,  dpuv_E])
    #             # print([dpuv_SW, dpuv_S, dpuv_SW])
    #             # print(S.evalpts[(S.sample_size_u * u) + v])
    #             # print("========================")
    #
    #             if maxS:
    #                 if ((dpuv >= dpuv_NW) and (dpuv >= dpuv_N) and (dpuv >= dpuv_NE) and
    #                         (dpuv >= dpuv_W) and (dpuv >= dpuv_E) and
    #                         (dpuv >= dpuv_SW) and (dpuv >= dpuv_S) and (dpuv >= dpuv_SE)):
    #                     extremaUV.append((u, v,))
    #             else:  # minS
    #                 if ((dpuv <= dpuv_NW) and (dpuv <= dpuv_N) and (dpuv <= dpuv_NE) and
    #                         (dpuv <= dpuv_W) and (dpuv <= dpuv_E) and
    #                         (dpuv <= dpuv_SW) and (dpuv <= dpuv_S) and (dpuv <= dpuv_SE)):
    #                     # where a point is orthogonal to a planar surface -> sphere radius test, or accept minima
    #                     extremaUV.append((u, v,))
    #
    #     return extremaUV

    if localExtrema and not curvatureTest:
        # change to internal function
        for u in range(1, evalptsU - 1):
            for v in range(1, evalptsV - 1):

                # print([ evalptsU * (u + 1) + (v - 1), evalptsU * (u + 1) + v, evalptsU * (u + 1) + (v + 1)])
                # print([evalptsU * u + (v - 1), evalptsU * u + v, evalptsU * u + (v + 1)])
                # print([ evalptsU * (u - 1) + (v - 1), evalptsU * (u - 1) + v, evalptsU * (u - 1) + (v + 1)])

                if v > 1:
                    dpuv_NW = dpuv_N
                    dpuv_W = dpuv
                    dpuv_SW = dpuv_S

                    dpuv_N = dpuv_NE
                    dpuv = dpuv_E
                    dpuv_S = dpuv_SE
                else:
                    dpuv_NW = np.linalg.norm(p - S.evalpts[evalptsU * (u + 1)])  # + (v - 1)])
                    dpuv_W = np.linalg.norm(p - S.evalpts[evalptsU * u])  # + (v - 1)])
                    dpuv_SW = np.linalg.norm(p - S.evalpts[evalptsU * (u - 1)])  # + (v - 1)])

                    dpuv_N = np.linalg.norm(p - S.evalpts[(evalptsU * (u + 1)) + 1])  # v])
                    dpuv = np.linalg.norm(p - S.evalpts[(evalptsU * u) + 1])  # v])
                    dpuv_S = np.linalg.norm(p - S.evalpts[(evalptsU * (u - 1)) + 1])  # v])

                # only have to calculate this once
                dpuv_NE = np.linalg.norm(p - S.evalpts[(evalptsU * (u + 1)) + v + 1])
                dpuv_E = np.linalg.norm(p - S.evalpts[(evalptsU * u) + v + 1])
                dpuv_SE = np.linalg.norm(p - S.evalpts[(evalptsU * (u - 1)) + v + 1])

                # print([dpuv_NW, dpuv_N, dpuv_NE])
                # print([dpuv_W , dpuv ,  dpuv_E])
                # print([dpuv_SW, dpuv_S, dpuv_SW])
                # print("========================")

                if maxS:
                    if ((dpuv >= dpuv_NW) and (dpuv >= dpuv_N) and (dpuv >= dpuv_NE) and
                            (dpuv >= dpuv_W) and (dpuv >= dpuv_E) and
                            (dpuv >= dpuv_SW) and (dpuv >= dpuv_S) and (dpuv >= dpuv_SE)):
                        # where a point is orthogonal to a
                        localExtremaUV.append((u, v,))
                else:  # minS
                    if ((dpuv <= dpuv_NW) and (dpuv <= dpuv_N) and (dpuv <= dpuv_NE) and
                            (dpuv <= dpuv_W) and (dpuv <= dpuv_E) and
                            (dpuv <= dpuv_SW) and (dpuv <= dpuv_S) and (dpuv <= dpuv_SE)):

                        # print([dpuv_NW, dpuv_N, dpuv_NE])
                        # print([dpuv_W, dpuv, dpuv_E])
                        # print([dpuv_SW, dpuv_S, dpuv_SW])
                        # print(S.evalpts[(evalptsU * u) + v])
                        # print("========================")
                        localExtremaUV.append((u, v,))

    # def curvatureSearch():
    #     discreteK = [np.inf] * len(S.evalpts)
    #
    #     if closedU():
    #         rangeU = (0, S.sample_size_u)
    #     else:
    #         rangeU = (1, S.sample_size_u - 1)
    #     if closedV():
    #         rangeV = (0, S.sample_size_v)
    #     else:
    #         rangeV = (1, S.sample_size_v - 1)
    #     # get all local curvatures
    #     for u in range(rangeU[0], rangeU[1]):
    #         for v in range(rangeV[0], rangeV[1]):
    #
    #             if v > 1:
    #                 puv_NW = puv_N
    #                 puv_W = puv
    #                 puv_SW = puv_S
    #
    #                 puv_N = puv_NE
    #                 puv = puv_E
    #                 puv_S = puv_SE
    #             else:
    #                 puv_NW = np.array(S.evalpts[S.sample_size_u * (u + 1)])  # + (v - 1)])
    #                 puv_W = np.array(S.evalpts[S.sample_size_u * u])  # + (v - 1)])
    #                 puv_SW = np.array(S.evalpts[S.sample_size_u * (u - 1)])  # + (v - 1)])
    #
    #                 puv_N = np.array(S.evalpts[(S.sample_size_u * (u + 1)) + 1])  # v])
    #                 puv = np.array(S.evalpts[(S.sample_size_u * u) + 1])  # v])
    #                 puv_S = np.array(S.evalpts[(S.sample_size_u * (u - 1)) + 1])  # v])
    #
    #             # only have to calculate this once
    #             puv_NE = np.array(S.evalpts[(S.sample_size_u * (u + 1)) + v + 1])
    #             puv_E = np.array(S.evalpts[(S.sample_size_u * u) + v + 1])
    #             puv_SE = np.array(S.evalpts[(S.sample_size_u * (u - 1)) + v + 1])
    #
    #             pur, pucc = radiusCentre3points(puv_W, puv, puv_E)
    #             pvr, pvcc = radiusCentre3points(puv_S, puv, puv_N)
    #
    #             # S.evalpts[(S.sample_size_u * u) + v]
    #
    #             # find media point of U curvature centre and V curvature, kuv
    #             if pur is np.inf and pvr is not np.inf:
    #                 kuv = pvcc - puv
    #             elif pur is not np.inf and pvr is np.inf:
    #                 kuv = pucc - puv
    #             elif pur is np.inf and pvr is np.inf:
    #                 discreteK[(S.sample_size_u * u) + v] = np.inf
    #                 break
    #             else:
    #                 kuv = ((pucc - puv) + (pvcc - puv)) / 2
    #
    #             dkuv = np.linalg.norm(kuv - puv)
    #
    #             if np.linalg.norm(p - puv) > np.linalg.norm(kuv - p):
    #                 # p is on the same side of the surface as the median curvature
    #                 discreteK[(S.sample_size_u * u) + v] = dkuv
    #
    #             if np.linalg.norm(p - puv) < np.linalg.norm(kuv - p):
    #                 # smallest radius, with curvature centre on far side of surface from p
    #                 discreteK[(S.sample_size_u * u) + v] = -dkuv
    #
    #     return discreteK
    #
    # def globalCurvatureExtrema(discreteK):
    #     d_min = np.inf
    #     for u in range(0, S.sample_size_u):
    #         for v in range(0, S.sample_size_v):
    #
    #             dkuv = discreteK[(S.sample_size_u * u) + v]
    #
    #             if maxS and (dkuv < 0.0): # maxima
    #                 if -dkuv < d_min:
    #                     d_min = dkuv
    #                     leu = (u, v,)
    #
    #             elif not maxS and (dkuv > 0.0): # minima
    #                 if dkuv < d_min:
    #                     d_min = dkuv
    #                     leu = (u, v,)
    #     return leu
    #
    #                 # if not localExtrema:
    #                 #     if maxS:
    #                 #         # spline version used a maxima of 3 points, but as outlying points are averages, this seems excessive
    #                 #         # dua = max([np.linalg.norm(p - pu) for pu in [puLower, tess.points[tuv][0], puHigher]])
    #                 #         # dva = max([np.linalg.norm(p - pv) for pv in [pvLower, tess.points[tuv][1], pvHigher]])
    #                 #
    #                 #         # centre of 3 point circle is nearer or further than max/min sample point from p
    #                 #         #orthoSign = dua > (dcc_p + (2*dcc)) # circle centre closer to p than circle points => maximum curve
    #                 #
    #                 #         # note that there is a difference between the unique minima/maxima point of any surface and a
    #                 #         # collection of dimples or bumps that identify a surface
    #                 #
    #                 #         if (dpuvcc < dpuv) and (np.abs(duvk) < d_min): # maxima, > dpuv for minima
    #                 #             d_min = np.abs(duvk)
    #                 #             # instead rank curvature?
    #                 #             #localExtremaUV.append((u * S.delta_u, v * S.delta_v,))
    #                 #             localExtremaUnique = (u, v,)
    #                 #
    #                 #     else: #minS
    #                 #         if (dpuvcc > dpuv) and (np.abs(duvk) < d_min): # minima
    #                 #             d_min = np.abs(duvk)
    #                 #             #localExtremaUV.append((u * S.delta_u, v * S.delta_v,))
    #                 #             localExtremaUnique = (u, v,)

    # def localCurvatureExtrema(discreteK):
    #     leu = []
    #
    #     if closedU():
    #         rangeU = (0, S.sample_size_u)
    #     else:
    #         rangeU = (1, S.sample_size_u - 1)
    #     if closedV():
    #         rangeV = (0, S.sample_size_v)
    #     else:
    #         rangeV = (1, S.sample_size_v - 1)
    #     # get all local curvatures
    #     for u in range(rangeU[0], rangeU[1]):
    #         for v in range(rangeV[0], rangeV[1]):
    #
    #             surroundNodes = [
    #                 # S.sample_size_u * (u + 1) + (v - 1),
    #                 S.sample_size_u * (u + 1) + v,
    #                 # S.sample_size_u * (u + 1) + (v + 1),
    #                 S.sample_size_u * u + (v - 1),
    #                 # S.sample_size_u * u + v,
    #                 S.sample_size_u * u + (v + 1),
    #                 # S.sample_size_u * (u - 1) + (v - 1),
    #                 S.sample_size_u * (u - 1) + v,
    #                 # S.sample_size_u * (u - 1) + (v + 1)
    #             ]
    #
    #             if not any(discreteK[ssn] == np.inf for ssn in surroundNodes):
    #                 if maxS:
    #                     if all((np.abs(discreteK[ssn]) >= discreteK[(S.sample_size_u * u) + v]) for ssn in
    #                            surroundNodes) and (discreteK[(S.sample_size_u * u) + v] > 0):
    #                         # this version excludes negative curves from surrounding curvature values
    #                         # if all(np.abs(discreteK[ssn]) > discreteK[S.sample_size_u * u + v] or
    #                         # (discreteK[ssn]<0) for ssn in surroundNodes) and
    #                         # (discreteK[S.sample_size_u * u + v] > 0):
    #                         leu.append((u, v,))
    #                 else:
    #                     if all((np.abs(discreteK[ssn]) >= discreteK[(S.sample_size_u * u) + v]) for ssn in
    #                            surroundNodes) and (discreteK[(S.sample_size_u * u) + v] < 0):
    #                         leu.append((u, v,))
    #     return leu

    if curvatureTest: # complex test summing curvatures along U and V axes
        localCurvature = [np.inf] * len(S.evalpts)

        # get all local curvatures
        for u in range(1, evalptsU - 1):
            for v in range(1, evalptsV - 1):

                if v > 1:
                    puv_NW = puv_N
                    puv_W = puv
                    puv_SW = puv_S

                    puv_N = puv_NE
                    puv = puv_E
                    puv_S = puv_SE
                else:
                    puv_NW = np.array(S.evalpts[evalptsU * (u + 1)])  # + (v - 1)])
                    puv_W = np.array(S.evalpts[evalptsU * u])  # + (v - 1)])
                    puv_SW = np.array(S.evalpts[evalptsU * (u - 1)])  # + (v - 1)])

                    puv_N = np.array(S.evalpts[(evalptsU * (u + 1)) + 1])  # v])
                    puv = np.array(S.evalpts[(evalptsU * u) + 1])  # v])
                    puv_S = np.array(S.evalpts[(evalptsU * (u - 1)) + 1])  # v])

                # only have to calculate this once
                puv_NE = np.array(S.evalpts[(evalptsU * (u + 1)) + v + 1])
                puv_E = np.array(S.evalpts[(evalptsU * u) + v + 1])
                puv_SE = np.array(S.evalpts[(evalptsU * (u - 1)) + v + 1])

                if not all(type(t[0])==float and type(t[1])==float and type(t[2])==float for t in [puv_N, puv_S, puv_E, puv_W, puv]):
                    _=1

                #print([puv_W, puv, puv_E])

                pur, pucc = radiusCentre3points(puv_W, puv, puv_E)
                pvr, pvcc = radiusCentre3points(puv_S, puv, puv_N)

                #S.evalpts[(evalptsU * u) + v]

                # find media point of U curvature centre and V curvature, kuv
                if pur is np.inf and pvr is not np.inf:
                    kuv = pvcc - puv
                elif pur is not np.inf and pvr is np.inf:
                    kuv = pucc - puv
                elif pur is np.inf and pvr is np.inf:
                    localCurvature[(evalptsU * u) + v] = np.inf
                    break
                else:
                    kuv = ((pucc - puv) + (pvcc - puv))/2

                dkuv = np.linalg.norm(kuv - puv)

                if np.linalg.norm(p - puv) > np.linalg.norm(kuv - p):
                    # p is on the same side of the surface as the median curvature
                    localCurvature[(evalptsU * u) + v] = dkuv
                    if maxS and dkuv < d_min:  # maxima
                        d_min = dkuv
                        # instead rank curvature?
                        localExtremaUnique = (u, v,)

                if np.linalg.norm(p - puv) < np.linalg.norm(kuv - p):
                    # smallest radius, with curvature centre on far side of surface from p
                    localCurvature[(evalptsU * u) + v] = -dkuv
                    if not maxS and dkuv < d_min:
                        d_min = dkuv
                        # instead rank curvature?
                        localExtremaUnique = (u, v,)

                    # if not localExtrema:
                    #     if maxS:
                    #         # spline version used a maxima of 3 points, but as outlying points are averages, this seems excessive
                    #         # dua = max([np.linalg.norm(p - pu) for pu in [puLower, tess.points[tuv][0], puHigher]])
                    #         # dva = max([np.linalg.norm(p - pv) for pv in [pvLower, tess.points[tuv][1], pvHigher]])
                    #
                    #         # centre of 3 point circle is nearer or further than max/min sample point from p
                    #         #orthoSign = dua > (dcc_p + (2*dcc)) # circle centre closer to p than circle points => maximum curve
                    #
                    #         # note that there is a difference between the unique minima/maxima point of any surface and a
                    #         # collection of dimples or bumps that identify a surface
                    #
                    #         if (dpuvcc < dpuv) and (np.abs(duvk) < d_min): # maxima, > dpuv for minima
                    #             d_min = np.abs(duvk)
                    #             # instead rank curvature?
                    #             #localExtremaUV.append((u * S.delta_u, v * S.delta_v,))
                    #             localExtremaUnique = (u, v,)
                    #
                    #     else: #minS
                    #         if (dpuvcc > dpuv) and (np.abs(duvk) < d_min): # minima
                    #             d_min = np.abs(duvk)
                    #             #localExtremaUV.append((u * S.delta_u, v * S.delta_v,))
                    #             localExtremaUnique = (u, v,)

        # test individual
        if localExtrema:

            for u in range(1, evalptsU - 1):
                for v in range(1, evalptsV - 1):
                    surroundNodes = [ #evalptsU * (u + 1) + (v - 1),
                                      evalptsU * (u + 1) + v,
                                      #evalptsU * (u + 1) + (v + 1),
                                      evalptsU * u + (v - 1),
                                      #evalptsU * u + v,
                                      evalptsU * u + (v + 1),
                                      #evalptsU * (u - 1) + (v - 1),
                                      evalptsU * (u - 1) + v,
                                      #evalptsU * (u - 1) + (v + 1)
                                    ]

                    if not any(localCurvature[ssn] == np.inf for ssn in surroundNodes):
                        if maxS:

                            if all((np.abs(localCurvature[ssn]) >= localCurvature[evalptsU * u + v]) for ssn in surroundNodes) and (localCurvature[evalptsU * u + v] > 0):
                                # this version excludes negative curves from surrounding curvature values
                            #if all(np.abs(localCurvature[ssn]) > localCurvature[evalptsU * u + v] or (localCurvature[ssn]<0) for ssn in surroundNodes) and (localCurvature[evalptsU * u + v] > 0):
                                localExtremaUV.append((u, v,))
                        else:
                            if all((np.abs(localCurvature[ssn]) >= localCurvature[evalptsU * u + v]) for ssn in surroundNodes) and (localCurvature[evalptsU * u + v] < 0):
                                localExtremaUV.append((u, v,))

    def subSearchUV(mpU, mpV, tol):
        tolFlag = 0
        divU = S.delta_u
        divV = S.delta_v
        dpuvlocal = np.linalg.norm(p - S.evaluate_single((mpU, mpV,)))

        while not tolFlag:
            # create simple subdivisions around a minimum point as simple hillclimb search

            divU /= 2
            divV /= 2

            subDivs = [(mpU-divU, mpV+divV),
                       (mpU,      mpV+divV),
                       (mpU+divU, mpV+divV),
                       (mpU-divU, mpV),
                       (mpU,      mpV),
                       (mpU+divU, mpV),
                       (mpU-divU, mpV-divV),
                       (mpU,      mpV-divV),
                       (mpU+divU, mpV-divV)]

            pSubDivs = S.evaluate_list(subDivs)

            dpuv_NW = np.linalg.norm(p - pSubDivs[0])
            dpuv_N = np.linalg.norm(p - pSubDivs[1])
            dpuv_NE = np.linalg.norm(p - pSubDivs[2])

            dpuv_W = np.linalg.norm(p - pSubDivs[3])
            #dpuv = np.linalg.norm(p - pSubDivs[4])
            dpuv_E = np.linalg.norm(p - pSubDivs[5])

            dpuv_SW = np.linalg.norm(p - pSubDivs[6])
            dpuv_S = np.linalg.norm(p - pSubDivs[7])
            dpuv_SE = np.linalg.norm(p - pSubDivs[8])

            dispUV = [dpuv_NW,
                      dpuv_N,
                      dpuv_NE,
                      dpuv_W,
                      dpuvlocal,
                      dpuv_E,
                      dpuv_SW,
                      dpuv_S,
                      dpuv_SE]

            if maxS:
                if np.abs(max(dispUV) - dpuvlocal) >= tol:
                    maxUVindex = dispUV.index(max(dispUV))
                    mpU = subDivs[maxUVindex][0]
                    mpV = subDivs[maxUVindex][1]
                    dpuvlocal = max(dispUV)
                else:
                    tolFlag +=1
            else:
                if np.abs(min(dispUV) - dpuvlocal) >= tol:
                    minUVindex = dispUV.index(min(dispUV))
                    mpU = subDivs[minUVindex][0]
                    mpV = subDivs[minUVindex][1]
                    dpuvlocal = min(dispUV)
                else:
                    tolFlag +=1
        return (mpU, mpV,)

    def testConditions(cuv, maxits=5):

        # def f(uv):
        #     return surface.derivatives(uv[0], uv[1], 2)

        def n(uv, e, r):  # np.array inputs
            #   f = Su(u,v) * r = 0
            #   g = Sv(u,v) * r = 0

            Su = e[1][0]
            Sv = e[0][1]

            Suu = e[2][0]
            Svv = e[0][2]

            Suv = e[1][1]
            Svu = e[1][1]

            f = np.dot(Su, r)
            g = np.dot(Sv, r)
            k = [-f, -g]

            J00 = np.dot(Su, Su) + np.dot(Suu, r)
            J01 = np.dot(Su, Sv) + np.dot(Suv, r)
            J10 = np.dot(Su, Sv) + np.dot(Svu, r)
            J11 = np.dot(Sv, Sv) + np.dot(Svv, r)

            # d =   [ u* - u, v* - v ]
            # k = - [ f(u,v), g(u,v) ]
            # J =   |Su|^2   +  Suu * r       Su*Sv  +  Suv * r
            #        Su*Sv   +  Svu * r      |Sv|^2  +  Svv * r

            J = [[J00, J01], [J10, J11]]
            # d = verb_core_Mat.solve(J,k)
            d = np.linalg.solve(J, k)
            return d + uv

        i = 0
        #e = None
        while i < maxits:
            e = np.array(S.derivatives(cuv[0], cuv[1], 2))
            dif = e[0][0] - p

            #   point coincidence   |S(u,v) - p| < e1
            c1v = np.linalg.norm(dif)
              # cosine
              #
              #  |Su(u,v)*(S(u,v) - P)|
              #  ----------------------  < e2
              #  |Su(u,v)| |S(u,v) - P|
              #
              #  |Sv(u,v)*(S(u,v) - P)|
              #  ----------------------  < e2
              #  |Sv(u,v)| |S(u,v) - P|

            c2an = np.dot( e[1][0], dif)
            c2ad = np.linalg.norm( e[1][0] ) * c1v
            c2bn = np.dot( e[0][1], dif)
            c2bd = np.linalg.norm( e[0][1] ) * c1v
            c2av = c2an / c2ad
            c2bv = c2bn / c2bd
            c1 = c1v < eps1
            c2a = c2av < eps2
            c2b = c2bv < eps2

            # exit if all tolerance are met,
            if (c1 and c2a and c2b):
                return cuv

            # otherwise, take a step
            ct = n(cuv, e, dif)

            #  correct for exceeding bounds
            if ct[0] < S.knotvector_u[0]:
                if closedU(): ct = [S.knotvector_u[-1]- (ct[0] - S.knotvector_u[0]), ct[1]]
                else: ct = [S.knotvector_u[0] + verb_EPSILON, ct[1]]
            elif ct[0] > S.knotvector_u[-1]:
                if closedU(): ct = [S.knotvector_u[0] + (ct[0] - S.knotvector_u[-1]), ct[1]]
                else: ct = [S.knotvector_u[-1] - verb_EPSILON, ct[1]]

            if ct[1] < S.knotvector_v[0]:
                if closedV(): ct = [ct[0], S.knotvector_v[-1] - ( ct[1] - S.knotvector_v[0] )]
                else: ct = [ct[0], S.knotvector_v[0] + verb_EPSILON]
            elif ct[1] > S.knotvector_v[-1]:
                if closedV(): ct = [ct[0], S.knotvector_v[0] + ( ct[0] - S.knotvector_v[-1] )]
                else: ct = [ct[0], S.knotvector_v[-1] - verb_EPSILON]

            c3v0 =  np.linalg.norm( (ct[0] - cuv[0]) * e[1][0] )
            c3v1 =  np.linalg.norm( (ct[1] - cuv[1]) * e[0][1] )

            # if |(u* - u) C'(u)| < e1, halt
            if (c3v0 + c3v1 < eps1):
                return cuv
                #break
            cuv = ct
            i = (i + 1)

        if i == maxits: # Newton-Raphson fails
            return (np.inf, np.inf)
        return cuv


    # if not localExtrema and not curvatureTest:
    #     localExtremaUnique = globalNeighbourSearch()
    #     if localExtrema:
    #         localExtremaUV = localNeighbourSearch()
    #
    # if curvatureTest: # discrete grid summing curvatures along U and V axes
    #     localCurvature = curvatureSearch()
    #     if localExtrema:
    #         localExtremaUV = localCurvatureExtrema(localCurvature)
    #     else:
    #         localExtremaUV = globalCurvatureExtrema(localCurvature)

    if localExtrema: # get the set of local extrema rather than global min/max of surface
        extremaUV = [] #
        for luv in localExtremaUV:
            print(luv[0] * S.delta_u, luv[1] * S.delta_v,)
            # tc = testConditions((luv[0] * S.delta_u, luv[1] * S.delta_v,))
            # print(S.evaluate_single(tc))
            mUV = testConditions((luv[0] * S.delta_u, luv[1] * S.delta_v,))
            if np.inf not in mUV:
                extremaUV.append(mUV)
            else: # Newton-Raphson fail, try hillclimb
                mUV = subSearchUV(luv[0] * S.delta_u, luv[1] * S.delta_v, verb_EPSILON)
                extremaUV.append(mUV)
            print(S.evalpts[(S.sample_size_u * luv[0]) + luv[1]])
        # filter out identical values
        localExtremaUnique = []
        for m in extremaUV:
            # displacement between any 2 points < eps1
            # compare u, v individually to reduce comparison pool
            if not any(np.isclose(m, u, eps1).all() for u in localExtremaUnique):
                localExtremaUnique.append(m)
        if uv_xyz:
            return localExtremaUnique #, [tess.uvs[le] for le in localExtremaSet]
        else:
            return S.evaluate_list(localExtremaUnique)
    else:
        print(localExtremaUnique[0] * S.delta_u, localExtremaUnique[1] * S.delta_v, )
        if uv_xyz:
            return testConditions((localExtremaUnique[0] * S.delta_u, localExtremaUnique[1] * S.delta_v,))
        else:
            localExtremaUnique = testConditions((localExtremaUnique[0] * S.delta_u, localExtremaUnique[1] * S.delta_v,))
            return S.evaluate_single(localExtremaUnique)



def rationalSurfaceExtremaParam_3(S, p, maxS=True, localExtrema=False, curvatureTest=False, uv_xyz=True, eps1=0.0001, eps2=0.0005, delta=0.025):
    '''
    Use either Newton-Raphson, or a hillclimbing neighbour search to find minimal or maximal distance on a surface, S, from a supplied point, p
    - not tested with surfaces closed along U & V
    maxS: search for most distant point from p if true, else nearest if false
    localExtrema: find any local extrema which differs from surrounding sample points, tolerance not sufficiently defined
    curvatureTest: rather than displacement of sampled points to provide seed locations for N-R/hillclimb search,
    use a combination of the centres of curvature at any point across U-axis and V-axis to identify C2 inflection point.
    uv_xyz: return u,v normalised values or cartesian x,y,z values
    eps1: point coincidence limit, halt N-R on reaching this minima
    eps2: cosine angle limit, halt N-R on reaching this minima
    delta: sample interval dividing surface into points
    '''

    # easy to find a minimum, but a maximum found through the original hillclimbing algorithm will settle on local maxima
    # for surfaces, minimize the following:
    #
    # f = Su(u,v) * r = 0
    # g = Sv(u,v) * r = 0
    #
    #   where r = S(u,v) - P
    #
    # Requires Newton-Raphson iteration, objective function is vector valued
    #
    #     J d = k
    #
    #     d =   [ u* - u, v* - v ] (alternatively u[i+1] -  u[i], v[i+1]- v[i])
    #     k = - [ f(u,v), g(u,v) ]
    #     J =     |Su|^2   +  Suu * r       Su*Sv  +  Suv * r
    #              Su*Sv   +  Svu * r      |Sv|^2  +  Svv * r
    #
    # halting conditions:
    # point coincidence
    #
    #         |S(u,v) - p| < e1
    #
    # cosine
    #
    #    |Su(u,v)*(S(u,v) - P)|
    #    ----------------------  < e2
    #    |Su(u,v)| |S(u,v) - P|
    #
    #    |Sv(u,v)*(S(u,v) - P)|
    #    ----------------------  < e2
    #    |Sv(u,v)| |S(u,v) - P|
    #
    # 1) first check 2 & 3
    # 2) if at least one of these is not, compute new value, otherwise halt
    # 3) ensure the parameter stays within range
    #     * if not closed, don't allow outside of range a-b
    #     * if closed (e.g. circle), allow to move back to beginning
    # 4)  if |(u* - u)C'(u)| < e1, halt

    # todo: take value from STEP fields
    # check if surface closed along U and V
    def closedU(): return np.isclose(np.array(S.ctrlpts2d[0]), np.array(S.ctrlpts2d[-1]), verb_EPSILON).all()  # eps_STEP_AP21
    def closedV(): return np.isclose(np.array(S.ctrlpts2d[0]).T, np.array(S.ctrlpts2d[-1]).T, verb_EPSILON).all()  # eps_STEP_AP21

    S.delta = delta # set evaluation delta
    # tradeoff between sample density and feature identification

    if S.evalpts == None: # evaluate surface points---------appears redundant
        S.evaluate()

    def globalNeighbourSearch():
        # remove values at edges of surface
        if closedU() and not closedV():
            edgeIndex = ([r for r in range(0, S.sample_size_v)] +
                         [r for r in range(S.sample_size_v * S.sample_size_u - S.sample_size_v, S.sample_size_v * S.sample_size_u)])
        elif not closedU() and closedV():
            edgeIndex = ([r for r in range(0, S.sample_size_v * S.sample_size_u) if r % S.sample_size_v == 0] +
                         [r for r in range(0, S.sample_size_v * S.sample_size_u) if r % S.sample_size_v == S.sample_size_v - 1])
        elif not closedU() and not closedV():
            edgeIndex = ([r for r in range(0, S.sample_size_v)] +
             [r for r in range(0, S.sample_size_v * S.sample_size_u) if r % S.sample_size_v == 0] +
             [r for r in range(0, S.sample_size_v * S.sample_size_u) if r % S.sample_size_v == S.sample_size_v - 1] +
             [r for r in range(S.sample_size_v * S.sample_size_u - S.sample_size_v, S.sample_size_v * S.sample_size_u)])

        edgeIndex = list(set(edgeIndex))

        edgeFreeIndex = [efi for efi in range(0, S.sample_size_v * S.sample_size_u) if efi not in edgeIndex]

        extremaUnique = [np.linalg.norm(p - S.evalpts[efi]) for efi in edgeFreeIndex]

        if maxS:
            extremaUnique = edgeFreeIndex[extremaUnique.index(max(extremaUnique))]
        else:
            extremaUnique = edgeFreeIndex[extremaUnique.index(min(extremaUnique))]

        return (extremaUnique // S.sample_size_v, extremaUnique % S.sample_size_v)

    def localNeighbourSearch():
        extremaUV = []

        if closedU():
            rangeU = (0, S.sample_size_u)
        else:
            rangeU = (1, S.sample_size_u - 1)
        if closedV():
            rangeV = (0, S.sample_size_v)
        else:
            rangeV = (1, S.sample_size_v - 1)
        # get all local curvatures
        for u in range(rangeU[0], rangeU[1]):
            for v in range(rangeV[0], rangeV[1]):

                # print([ S.sample_size_u * (u + 1) + (v - 1), S.sample_size_u * (u + 1) + v, S.sample_size_u * (u + 1) + (v + 1)])
                # print([S.sample_size_u * u + (v - 1), S.sample_size_u * u + v, S.sample_size_u * u + (v + 1)])
                # print([ S.sample_size_u * (u - 1) + (v - 1), S.sample_size_u * (u - 1) + v, S.sample_size_u * (u - 1) + (v + 1)])

                if v > 1:
                    dpuv_NW = dpuv_N
                    dpuv_W = dpuv
                    dpuv_SW = dpuv_S

                    dpuv_N = dpuv_NE
                    dpuv = dpuv_E
                    dpuv_S = dpuv_SE
                else:
                    dpuv_NW = np.linalg.norm(p - S.evalpts[S.sample_size_u * (u + 1)])  # + (v - 1)])
                    dpuv_W = np.linalg.norm(p - S.evalpts[S.sample_size_u * u])  # + (v - 1)])
                    dpuv_SW = np.linalg.norm(p - S.evalpts[S.sample_size_u * (u - 1)])  # + (v - 1)])

                    dpuv_N = np.linalg.norm(p - S.evalpts[(S.sample_size_u * (u + 1)) + 1])  # v])
                    dpuv = np.linalg.norm(p - S.evalpts[(S.sample_size_u * u) + 1])  # v])
                    dpuv_S = np.linalg.norm(p - S.evalpts[(S.sample_size_u * (u - 1)) + 1])  # v])

                # only have to calculate this once
                dpuv_NE = np.linalg.norm(p - S.evalpts[(S.sample_size_u * (u + 1)) + v + 1])
                dpuv_E = np.linalg.norm(p - S.evalpts[(S.sample_size_u * u) + v + 1])
                dpuv_SE = np.linalg.norm(p - S.evalpts[(S.sample_size_u * (u - 1)) + v + 1])

                # print([dpuv_NW, dpuv_N, dpuv_NE])
                # print([dpuv_W , dpuv ,  dpuv_E])
                # print([dpuv_SW, dpuv_S, dpuv_SW])
                # print(S.evalpts[(S.sample_size_u * u) + v])
                # print("========================")

                if maxS:
                    if ((dpuv >= dpuv_NW) and (dpuv >= dpuv_N) and (dpuv >= dpuv_NE) and
                            (dpuv >= dpuv_W) and (dpuv >= dpuv_E) and
                            (dpuv >= dpuv_SW) and (dpuv >= dpuv_S) and (dpuv >= dpuv_SE)):
                        extremaUV.append((u, v,))
                else:  # minS
                    if ((dpuv <= dpuv_NW) and (dpuv <= dpuv_N) and (dpuv <= dpuv_NE) and
                            (dpuv <= dpuv_W) and (dpuv <= dpuv_E) and
                            (dpuv <= dpuv_SW) and (dpuv <= dpuv_S) and (dpuv <= dpuv_SE)):
                        # where a point is orthogonal to a planar surface -> sphere radius test, or accept minima
                        extremaUV.append((u, v,))

        return extremaUV

    def curvatureSearch():
        discreteK = [np.inf] * len(S.evalpts)

        if closedU():
            rangeU = (0, S.sample_size_u)
        else:
            rangeU = (1, S.sample_size_u - 1)
        if closedV():
            rangeV = (0, S.sample_size_v)
        else:
            rangeV = (1, S.sample_size_v - 1)
        # get all local curvatures
        for u in range(rangeU[0], rangeU[1]):
            for v in range(rangeV[0], rangeV[1]):

                if v > 1:
                    puv_NW = puv_N
                    puv_W = puv
                    puv_SW = puv_S

                    puv_N = puv_NE
                    puv = puv_E
                    puv_S = puv_SE
                else:
                    puv_NW = np.array(S.evalpts[S.sample_size_u * (u + 1)])  # + (v - 1)])
                    puv_W = np.array(S.evalpts[S.sample_size_u * u])  # + (v - 1)])
                    puv_SW = np.array(S.evalpts[S.sample_size_u * (u - 1)])  # + (v - 1)])

                    puv_N = np.array(S.evalpts[(S.sample_size_u * (u + 1)) + 1])  # v])
                    puv = np.array(S.evalpts[(S.sample_size_u * u) + 1])  # v])
                    puv_S = np.array(S.evalpts[(S.sample_size_u * (u - 1)) + 1])  # v])

                # only have to calculate this once
                puv_NE = np.array(S.evalpts[(S.sample_size_u * (u + 1)) + v + 1])
                puv_E = np.array(S.evalpts[(S.sample_size_u * u) + v + 1])
                puv_SE = np.array(S.evalpts[(S.sample_size_u * (u - 1)) + v + 1])

                pur, pucc = radiusCentre3points(puv_W, puv, puv_E)
                pvr, pvcc = radiusCentre3points(puv_S, puv, puv_N)

                # S.evalpts[(S.sample_size_u * u) + v]

                # find media point of U curvature centre and V curvature, kuv
                if pur is np.inf and pvr is not np.inf:
                    kuv = pvcc - puv
                elif pur is not np.inf and pvr is np.inf:
                    kuv = pucc - puv
                elif pur is np.inf and pvr is np.inf:
                    discreteK[(S.sample_size_u * u) + v] = np.inf
                    break
                else:
                    kuv = ((pucc - puv) + (pvcc - puv)) / 2

                dkuv = np.linalg.norm(kuv - puv)

                if np.linalg.norm(p - puv) > np.linalg.norm(kuv - p):
                    # p is on the same side of the surface as the median curvature
                    discreteK[(S.sample_size_u * u) + v] = dkuv

                if np.linalg.norm(p - puv) < np.linalg.norm(kuv - p):
                    # smallest radius, with curvature centre on far side of surface from p
                    discreteK[(S.sample_size_u * u) + v] = -dkuv

        return discreteK

    def globalCurvatureExtrema(discreteK):
        d_min = np.inf
        for u in range(0, S.sample_size_u):
            for v in range(0, S.sample_size_v):

                dkuv = discreteK[(S.sample_size_u * u) + v]

                if maxS and (dkuv < 0.0): # maxima
                    if -dkuv < d_min:
                        d_min = dkuv
                        leu = (u, v,)

                elif not maxS and (dkuv > 0.0): # minima
                    if dkuv < d_min:
                        d_min = dkuv
                        leu = (u, v,)
        return leu

    def localCurvatureExtrema(discreteK):
        leu = []

        if closedU():
            rangeU = (0, S.sample_size_u)
        else:
            rangeU = (1, S.sample_size_u - 1)
        if closedV():
            rangeV = (0, S.sample_size_v)
        else:
            rangeV = (1, S.sample_size_v - 1)
        # get all local curvatures
        for u in range(rangeU[0], rangeU[1]):
            for v in range(rangeV[0], rangeV[1]):

                surroundNodes = [
                    # S.sample_size_u * (u + 1) + (v - 1),
                    S.sample_size_u * (u + 1) + v,
                    # S.sample_size_u * (u + 1) + (v + 1),
                    S.sample_size_u * u + (v - 1),
                    # S.sample_size_u * u + v,
                    S.sample_size_u * u + (v + 1),
                    # S.sample_size_u * (u - 1) + (v - 1),
                    S.sample_size_u * (u - 1) + v,
                    # S.sample_size_u * (u - 1) + (v + 1)
                ]

                if not any(discreteK[ssn] == np.inf for ssn in surroundNodes):
                    if maxS:
                        if all((np.abs(discreteK[ssn]) >= discreteK[(S.sample_size_u * u) + v]) for ssn in
                               surroundNodes) and (discreteK[(S.sample_size_u * u) + v] > 0):
                            # this version excludes negative curves from surrounding curvature values
                            # if all(np.abs(discreteK[ssn]) > discreteK[S.sample_size_u * u + v] or
                            # (discreteK[ssn]<0) for ssn in surroundNodes) and
                            # (discreteK[S.sample_size_u * u + v] > 0):
                            leu.append((u, v,))
                    else:
                        if all((np.abs(discreteK[ssn]) >= discreteK[(S.sample_size_u * u) + v]) for ssn in
                               surroundNodes) and (discreteK[(S.sample_size_u * u) + v] < 0):
                            leu.append((u, v,))
        return leu

    def subSearchUV(mpU, mpV, tol):
        # create simple subdivisions around a minimum point as simple hillclimb search
        tolFlag = 0
        divU = S.delta_u
        divV = S.delta_v
        dpuvlocal = np.linalg.norm(p - S.evaluate_single((mpU, mpV,)))

        while not tolFlag:
            divU /= 2
            divV /= 2

            subDivs = [(mpU-divU, mpV+divV),
                       (mpU,      mpV+divV),
                       (mpU+divU, mpV+divV),
                       (mpU-divU, mpV),
                       (mpU,      mpV),
                       (mpU+divU, mpV),
                       (mpU-divU, mpV-divV),
                       (mpU,      mpV-divV),
                       (mpU+divU, mpV-divV)]

            pSubDivs = S.evaluate_list(subDivs)

            dpuv_NW = np.linalg.norm(p - pSubDivs[0])
            dpuv_N = np.linalg.norm(p - pSubDivs[1])
            dpuv_NE = np.linalg.norm(p - pSubDivs[2])

            dpuv_W = np.linalg.norm(p - pSubDivs[3])
            #dpuv = np.linalg.norm(p - pSubDivs[4])
            dpuv_E = np.linalg.norm(p - pSubDivs[5])

            dpuv_SW = np.linalg.norm(p - pSubDivs[6])
            dpuv_S = np.linalg.norm(p - pSubDivs[7])
            dpuv_SE = np.linalg.norm(p - pSubDivs[8])

            dispUV = [dpuv_NW,
                      dpuv_N,
                      dpuv_NE,
                      dpuv_W,
                      dpuvlocal,
                      dpuv_E,
                      dpuv_SW,
                      dpuv_S,
                      dpuv_SE]

            if maxS:
                if np.abs(max(dispUV) - dpuvlocal) >= tol:
                    maxUVindex = dispUV.index(max(dispUV))
                    mpU = subDivs[maxUVindex][0]
                    mpV = subDivs[maxUVindex][1]
                    dpuvlocal = max(dispUV)
                else:
                    tolFlag +=1
            else:
                if np.abs(min(dispUV) - dpuvlocal) >= tol:
                    minUVindex = dispUV.index(min(dispUV))
                    mpU = subDivs[minUVindex][0]
                    mpV = subDivs[minUVindex][1]
                    dpuvlocal = min(dispUV)
                else:
                    tolFlag +=1
        return (mpU, mpV,)

    def testConditions(cuv, maxits=5):

        # def f(uv):
        #     return surface.derivatives(uv[0], uv[1], 2)

        def n(uv, e, r):  # np.array inputs
            #   f = Su(u,v) * r = 0
            #   g = Sv(u,v) * r = 0

            Su = e[1][0]
            Sv = e[0][1]

            Suu = e[2][0]
            Svv = e[0][2]

            Suv = e[1][1]
            Svu = e[1][1]

            f = np.dot(Su, r)
            g = np.dot(Sv, r)
            k = [-f, -g]

            J00 = np.dot(Su, Su) + np.dot(Suu, r)
            J01 = np.dot(Su, Sv) + np.dot(Suv, r)
            J10 = np.dot(Su, Sv) + np.dot(Svu, r)
            J11 = np.dot(Sv, Sv) + np.dot(Svv, r)

            # d =   [ u* - u, v* - v ]
            # k = - [ f(u,v), g(u,v) ]
            # J =   |Su|^2   +  Suu * r       Su*Sv  +  Suv * r
            #        Su*Sv   +  Svu * r      |Sv|^2  +  Svv * r

            J = [[J00, J01], [J10, J11]]
            # d = verb_core_Mat.solve(J,k)
            d = np.linalg.solve(J, k)
            return d + uv

        i = 0
        #e = None
        while i < maxits:
            e = np.array(S.derivatives(cuv[0], cuv[1], 2))
            dif = e[0][0] - p

            #   point coincidence   |S(u,v) - p| < e1
            c1v = np.linalg.norm(dif)
              # cosine
              #
              #  |Su(u,v)*(S(u,v) - P)|
              #  ----------------------  < e2
              #  |Su(u,v)| |S(u,v) - P|
              #
              #  |Sv(u,v)*(S(u,v) - P)|
              #  ----------------------  < e2
              #  |Sv(u,v)| |S(u,v) - P|

            c2an = np.dot( e[1][0], dif)
            c2ad = np.linalg.norm( e[1][0] ) * c1v
            c2bn = np.dot( e[0][1], dif)
            c2bd = np.linalg.norm( e[0][1] ) * c1v
            c2av = c2an / c2ad
            c2bv = c2bn / c2bd
            c1 = c1v < eps1
            c2a = c2av < eps2
            c2b = c2bv < eps2

            # exit if all tolerance are met,
            if (c1 and c2a and c2b):
                return cuv

            # otherwise, take a step
            ct = n(cuv, e, dif)

            #  correct for exceeding bounds
            if ct[0] < S.knotvector_u[0]:
                if closedU(): ct = [S.knotvector_u[-1]- (ct[0] - S.knotvector_u[0]), ct[1]]
                else: ct = [S.knotvector_u[0] + verb_EPSILON, ct[1]]
            elif ct[0] > S.knotvector_u[-1]:
                if closedU(): ct = [S.knotvector_u[0] + (ct[0] - S.knotvector_u[-1]), ct[1]]
                else: ct = [S.knotvector_u[-1] - verb_EPSILON, ct[1]]

            if ct[1] < S.knotvector_v[0]:
                if closedV(): ct = [ct[0], S.knotvector_v[-1] - ( ct[1] - S.knotvector_v[0] )]
                else: ct = [ct[0], S.knotvector_v[0] + verb_EPSILON]
            elif ct[1] > S.knotvector_v[-1]:
                if closedV(): ct = [ct[0], S.knotvector_v[0] + ( ct[0] - S.knotvector_v[-1] )]
                else: ct = [ct[0], S.knotvector_v[-1] - verb_EPSILON]

            c3v0 =  np.linalg.norm( (ct[0] - cuv[0]) * e[1][0] )
            c3v1 =  np.linalg.norm( (ct[1] - cuv[1]) * e[0][1] )

            # if |(u* - u) C'(u)| < e1, halt
            if (c3v0 + c3v1 < eps1):
                return cuv
                #break
            cuv = ct
            i = (i + 1)

        if i == maxits: # Newton-Raphson fails
            return (np.inf, np.inf)
        return cuv

    if not localExtrema and not curvatureTest:
        globalExtrema = globalNeighbourSearch()
    if localExtrema:
        localExtremaUV = localNeighbourSearch()

    if curvatureTest: # discrete grid summing curvatures along U and V axes
        localCurvature = curvatureSearch()
        if localExtrema:
            localExtremaUV = localCurvatureExtrema(localCurvature)
        else:
            localExtremaUV = globalCurvatureExtrema(localCurvature)

    if localExtrema: # get the set of local extrema rather than global min/max of surface
        extremaUV = []
        for luv in localExtremaUV:
            #print(luv[0] * S.delta_u, luv[1] * S.delta_v,)
            # tc = testConditions((luv[0] * S.delta_u, luv[1] * S.delta_v,))
            # print(S.evaluate_single(tc))
            mUV = testConditions((luv[0] * S.delta_u, luv[1] * S.delta_v,))
            if np.inf not in mUV:
                extremaUV.append(mUV)
            else: # Newton-Raphson fail, try hillclimb
                mUV = subSearchUV(luv[0] * S.delta_u, luv[1] * S.delta_v, verb_EPSILON)
                extremaUV.append(mUV)
            #print(S.evalpts[(S.sample_size_u * luv[0]) + luv[1]])

        # filter out identical values
        extremaUnique = []
        for m in extremaUV:
            # displacement between any 2 points < eps1
            # compare u, v individually to reduce comparison pool
            if not any(np.isclose(m, u, eps1).all() for u in extremaUnique):
                extremaUnique.append(m)

        #print(S.evaluate_list(extremaUnique))

        if uv_xyz:
            return extremaUnique #, [tess.uvs[le] for le in localExtremaSet]
        else:
            return S.evaluate_list(extremaUnique)
    else:
        #print(globalExtrema[0] * S.delta_u, globalExtrema[1] * S.delta_v, )
        if uv_xyz:
            return testConditions((globalExtrema[0] * S.delta_u, globalExtrema[1] * S.delta_v,))
        else:
            globalExtrema = testConditions((globalExtrema[0] * S.delta_u, globalExtrema[1] * S.delta_v,))
            return S.evaluate_single(globalExtrema)



#=======================================================================================================================

# # Control points
# ctrlpts = [
#     [[-25.0, -25.0, -10.0], [-25.0, -15.0, -5.0], [-25.0, -5.0, 0.0], [-25.0, 5.0, 0.0], [-25.0, 15.0, -5.0], [-25.0, 25.0, -10.0]],
#     [[-15.0, -25.0, -8.0], [-15.0, -15.0, -4.0], [-15.0, -5.0, -4.0], [-15.0, 5.0, -4.0], [-15.0, 15.0, -4.0], [-15.0, 25.0, -8.0]],
#     [[-5.0, -25.0, -5.0], [-5.0, -15.0, -3.0], [-5.0, -5.0, -8.0], [-5.0, 5.0, -8.0], [-5.0, 15.0, -3.0], [-5.0, 25.0, -5.0]],
#     [[5.0, -25.0, -3.0], [5.0, -15.0, -2.0], [5.0, -5.0, -8.0], [5.0, 5.0, -8.0], [5.0, 15.0, -2.0], [5.0, 25.0, -3.0]],
#     [[15.0, -25.0, -8.0], [15.0, -15.0, -4.0], [15.0, -5.0, -4.0], [15.0, 5.0, -4.0], [15.0, 15.0, -4.0], [15.0, 25.0, -8.0]],
#     [[25.0, -25.0, -10.0], [25.0, -15.0, -5.0], [25.0, -5.0, 2.0], [25.0, 5.0, 2.0], [25.0, 15.0, -5.0], [25.0, 25.0, -10.0]]
# ]
# # ctrlpts = [
# #     [-25.0, -25.0, -10.0], [-25.0, -15.0, -5.0], [-25.0, -5.0, 0.0], [-25.0, 5.0, 0.0], [-25.0, 15.0, -5.0], [-25.0, 25.0, -10.0],
# #     [-15.0, -25.0, -8.0], [-15.0, -15.0, -4.0], [-15.0, -5.0, -4.0], [-15.0, 5.0, -4.0], [-15.0, 15.0, -4.0], [-15.0, 25.0, -8.0],
# #     [-5.0, -25.0, -5.0], [-5.0, -15.0, -3.0], [-5.0, -5.0, -8.0], [-5.0, 5.0, -8.0], [-5.0, 15.0, -3.0], [-5.0, 25.0, -5.0],
# #     [5.0, -25.0, -3.0], [5.0, -15.0, -2.0], [5.0, -5.0, -8.0], [5.0, 5.0, -8.0], [5.0, 15.0, -2.0], [5.0, 25.0, -3.0],
# #     [15.0, -25.0, -8.0], [15.0, -15.0, -4.0], [15.0, -5.0, -4.0], [15.0, 5.0, -4.0], [15.0, 15.0, -4.0], [15.0, 25.0, -8.0],
# #     [25.0, -25.0, -10.0], [25.0, -15.0, -5.0], [25.0, -5.0, 2.0], [25.0, 5.0, 2.0], [25.0, 15.0, -5.0], [25.0, 25.0, -10.0]
# # ]
#
# # Create a BSpline surface
# #surface = BSpline.Surface(normalize_kv=False)
# #surf = BSpline.Surface(tessellate='make_quad_mesh')#tessellate
# #surf = BSpline.Surface()
# surf = NURBS.Surface()
#
# # Set degrees
# surf.degree_u = 3
# surf.degree_v = 3
#
# # STEP only uses unitary weights
#
# wctrlpts = copy.deepcopy(ctrlpts)
# # non-nested list
# if all([isinstance(wcc, float) for wc in wctrlpts for wcc in wc]):
#     ctrlpts = [wc+[1.0] for wc in wctrlpts if isinstance(wc, float)]
# else:
#     ctrlpts = [[wcc+[1.0] for wcc in wc if all([isinstance(wccc, float) for wccc in wcc])] for wc in wctrlpts]
#
# # Set control points
# surf.ctrlpts2d = ctrlpts # BSpline
# #surf.set_ctrlpts(ctrlpts, 6, 6) #NURBS homogeneous weighted control points
#
# # Set knot vectors
# surf.knotvector_u = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0]
# surf.knotvector_v = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0]
#
# # # Set evaluation delta
# # surf.delta = 0.025 # this seems to be the minima delta under adaptive tesselation
#
# testPoint = np.array([0., 0., 40.])

#=======================================================================================================================


# random surface generator
# https://nurbs-python.readthedocs.io/en/latest/surface_generator.html
# Generate a plane with the dimensions 50x100
surfgrid = CPGen.Grid(20, 20)

# Generate a grid of 25x30
surfgrid.generate(40, 40)

# Generate bumps on the grid
surfgrid.bumps(num_bumps=4, bump_height=-5, base_extent=5)

# p_file = open('surfgrid.pkl', 'wb')
# pickle.dump(surfgrid, p_file)
# p_file.close()

#f = open('surfgrid2.pkl','rb')
f = open('surfgrid.pkl','rb')
new_model = pickle.load(f)

# Create a BSpline surface instance
surf = BSpline.Surface()

# Set degrees
surf.degree_u = 3
surf.degree_v = 3

# Get the control points from the generated grid
surf.ctrlpts2d = surfgrid.grid

# Set knot vectors
surf.knotvector_u = utilities.generate_knot_vector(surf.degree_u, surf.ctrlpts_size_u)
surf.knotvector_v = utilities.generate_knot_vector(surf.degree_v, surf.ctrlpts_size_v)

# Set sample size
surf.sample_size = 200

#testPoint = np.array([0., 0., 40.])
testPoint = np.array([10., 10., 30.])

#=======================================================================================================================


xyzExtrema = rationalSurfaceExtremaParam_3(surf,
                                           testPoint,
                                           maxS=0,
                                           localExtrema=0,
                                           uv_xyz=0,
                                           curvatureTest=0)

# xyzExtremaG = rationalSurfaceExtremaParam_2(surf,
#                                            testPoint,
#                                            maxS=0,
#                                            localExtrema=False,
#                                            uv_xyz=0,
#                                            curvatureTest=0)

#xyzExtremaG in xyzExtrema

# [21.368049898892505, 1.996158427441941e-08, -0.4104889766482239]
# 0.925 0.47500000000000003

#==============================================

# minSegLE = rationalSurfaceExtremaParam_2(surf,
#                                    np.array([0., 0., 40.]),
#                                    maxS=False,
#                                    localExtrema=True,
#                                    uv_xyz=False,
#                                    curvatureTest=False)

# # 0.125 0.47500000000000003
# # 0.125 0.5
# # 0.55 0.2
# # 0.55 0.775
# # 0.925 0.47500000000000003
# # 0.925 0.5
# minKgen = rationalSurfaceExtremaParam_2(surf,
#                                    np.array([0., 0., 40.]),
#                                    maxS=False,
#                                    localExtrema=False,
#                                    uv_xyz=False,
#                                    curvatureTest=True)
# # 0.45 0.7250000000000001
# minKgen = rationalSurfaceExtremaParam_2(surf,
#                                    np.array([0., 0., 40.]),
#                                    maxS=False,
#                                    localExtrema=True,
#                                    uv_xyz=False,
#                                    curvatureTest=True)



# uvExtrema = rationalSurfaceExtremaParam(surf, p, maxSearch=True, localExtrema=True)
# MinPoints = [surf.derivatives(uvl[0], uvl[1], order=0)[0][0] for uvl in uvExtrema]
# # test global minima/maxima
# uvExtrema = rationalSurfaceExtremaParam(surf, p, maxSearch=True, localExtrema=False)
# MinPoint = surf.derivatives(uvExtrema[0], uvExtrema[1], order=0)[0][0]

# get global maxima from a local maxima set
dispsExtrema = [np.linalg.norm(np.array(m)- testPoint) for m in xyzExtrema]

MaxPointL = [xyzExtrema[dispsExtrema.index(min(dispsExtrema))], testPoint.tolist()]
NonMaxPointsL = sum([[m, testPoint.tolist()] for m in xyzExtrema], [])

# SearchPoints = [surf.derivatives(uvl[0], uvl[1], order=0)[0][0] for uvl in uvCoarse]
# SearchPointsL = sum([[m, p.tolist()] for m in SearchPoints], [])[:-1]

# Plot the control points grid and the evaluated surface
# Create a visualization configuration instance with no legend, no axes and set the resolution to 120 dpi

vis_config = VisMPL.VisConfig(alpha=0.4, display_axes=False, display_ctrlpts=False)
surf.vis = VisMPL.VisSurface(vis_config)
surf.vis.mconf['others']='points'
surf.vis.mconf['alpha']=0.1

pNonMaxPoints = dict(points=NonMaxPointsL, name="NonMaxPoints", color="green", size=1)
#pSearchPoints = dict(points=SearchPointsL, name="pSearchPoints", color="black", size=1)
pMaxPoints = dict(points=MaxPointL, name="MaxPoint", color="red", size=1)

#surf.render(extras=[pMinPoints, pSearchPoints], colormap=cm.cool)
surf.render(extras=[pNonMaxPoints, pMaxPoints], colormap=cm.cool)

pass

