from geomdl import CPGen
#from geomdl import BSpline
from geomdl import NURBS
import numpy as np
from geomdl import tessellate
from geomdl import operations as op
from geomdl import utilities
from geomdl.visualization import VisMPL
from geomdl.elements import Vertex, Triangle
import copy

# Import and use Matplotlib's colormaps
from matplotlib import cm
import matplotlib.colors as mcolors

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

def radiusCentre3points(p1, p2, p3):
    # radius + centre from 3 point circle

    t = p2 - p1
    u = p3 - p1
    v = p3 - p2

    w = np.cross(t, u)  # triangle normal
    wsl = np.linalg.norm(w)
    if (wsl < 10e-14): return False  # triangle area too small (additionally check points for colinearity)

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



# Control points
ctrlpts = [
    [[-25.0, -25.0, -10.0], [-25.0, -15.0, -5.0], [-25.0, -5.0, 0.0], [-25.0, 5.0, 0.0], [-25.0, 15.0, -5.0], [-25.0, 25.0, -10.0]],
    [[-15.0, -25.0, -8.0], [-15.0, -15.0, -4.0], [-15.0, -5.0, -4.0], [-15.0, 5.0, -4.0], [-15.0, 15.0, -4.0], [-15.0, 25.0, -8.0]],
    [[-5.0, -25.0, -5.0], [-5.0, -15.0, -3.0], [-5.0, -5.0, -8.0], [-5.0, 5.0, -8.0], [-5.0, 15.0, -3.0], [-5.0, 25.0, -5.0]],
    [[5.0, -25.0, -3.0], [5.0, -15.0, -2.0], [5.0, -5.0, -8.0], [5.0, 5.0, -8.0], [5.0, 15.0, -2.0], [5.0, 25.0, -3.0]],
    [[15.0, -25.0, -8.0], [15.0, -15.0, -4.0], [15.0, -5.0, -4.0], [15.0, 5.0, -4.0], [15.0, 15.0, -4.0], [15.0, 25.0, -8.0]],
    [[25.0, -25.0, -10.0], [25.0, -15.0, -5.0], [25.0, -5.0, 2.0], [25.0, 5.0, 2.0], [25.0, 15.0, -5.0], [25.0, 25.0, -10.0]]
]
# ctrlpts = [
#     [-25.0, -25.0, -10.0], [-25.0, -15.0, -5.0], [-25.0, -5.0, 0.0], [-25.0, 5.0, 0.0], [-25.0, 15.0, -5.0], [-25.0, 25.0, -10.0],
#     [-15.0, -25.0, -8.0], [-15.0, -15.0, -4.0], [-15.0, -5.0, -4.0], [-15.0, 5.0, -4.0], [-15.0, 15.0, -4.0], [-15.0, 25.0, -8.0],
#     [-5.0, -25.0, -5.0], [-5.0, -15.0, -3.0], [-5.0, -5.0, -8.0], [-5.0, 5.0, -8.0], [-5.0, 15.0, -3.0], [-5.0, 25.0, -5.0],
#     [5.0, -25.0, -3.0], [5.0, -15.0, -2.0], [5.0, -5.0, -8.0], [5.0, 5.0, -8.0], [5.0, 15.0, -2.0], [5.0, 25.0, -3.0],
#     [15.0, -25.0, -8.0], [15.0, -15.0, -4.0], [15.0, -5.0, -4.0], [15.0, 5.0, -4.0], [15.0, 15.0, -4.0], [15.0, 25.0, -8.0],
#     [25.0, -25.0, -10.0], [25.0, -15.0, -5.0], [25.0, -5.0, 2.0], [25.0, 5.0, 2.0], [25.0, 15.0, -5.0], [25.0, 25.0, -10.0]
# ]

# Create a BSpline surface
#surface = BSpline.Surface(normalize_kv=False)
#surf = BSpline.Surface(tessellate='make_quad_mesh')#tessellate
#surf = BSpline.Surface()
surf = NURBS.Surface()

# Set degrees
surf.degree_u = 3
surf.degree_v = 3

# STEP only uses unitary weights

wctrlpts = copy.deepcopy(ctrlpts)
# non-nested list
if all([isinstance(wcc, float) for wc in wctrlpts for wcc in wc]):
    ctrlpts = [wc+[1.0] for wc in wctrlpts if isinstance(wc, float)]
else:
    ctrlpts = [[wcc+[1.0] for wcc in wc if all([isinstance(wccc, float) for wccc in wcc])] for wc in wctrlpts]


# Set control points
surf.ctrlpts2d = ctrlpts # BSpline
#surf.set_ctrlpts(ctrlpts, 6, 6) #NURBS homogeneous weighted control points

# Set knot vectors
surf.knotvector_u = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0]
surf.knotvector_v = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0]

# Set evaluation delta
surf.delta = 0.025 # this seems to be the minima delta under adaptive tesselation

# Evaluate surface points
surf.evaluate()

# Tessellate surface
#surf.tessellate()


# # checking organisation of tessalation neighbours
#
# #from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.pyplot as plt
# fig = plt.figure()
# ax = fig.gca(projection='3d')
#
# tess_test = Tess.rationalSurfaceAdaptive(surf, AdaptiveRefinementOptions())
#
# for i in range(0, len(tess_test.points)):
#     label = '%d' % (i)
#     ax.text(tess_test.points[i][0], tess_test.points[i][1], tess_test.points[i][2], label)
#
# # # Tweaking display region and labels
# ax.set_xlim(-30, 30)
# ax.set_ylim(-30, 30)
# ax.set_zlim(-10, 10)
# ax.set_xlabel('X axis')
# ax.set_ylabel('Y axis')
# ax.set_zlabel('Z axis')
# plt.show()

# vsurface = verb.verb_geom_NurbsSurface.byKnotsControlPointsWeights( surface.degree_u,
#                                                                    surface.degree_v,
#                                                                    surface.knotvector_u,
#                                                                    surface.knotvector_v,
#                                                                    ctrlpts )
p = np.array([0, 0, 40])
#p = np.array([45, 30, -20])

#pv = vsurface.closestPoint( p )


# arrTrees = Tess.divideRationalSurfaceAdaptive(surf, AdaptiveRefinementOptions())
# tess = MeshData.empty()
# [aT.triangulate(tess) for aT in arrTrees]
# from matplotlib import cm
# from matplotlib.ticker import LinearLocator
# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# X = np.array([x[0] for x in tess.points])
# Y = np.array([y[1] for y in tess.points])
# Z = np.array([z[2] for z in tess.points])
# surf = ax.plot_trisurf(X, Y, Z, linewidth=10, cmap=plt.cm.CMRmap, antialiased=False)
# plt.show()
#
# surf.tessellator = tessellate.TriangularTessellate()

#uvg = rationalSurfaceClosestParam(surface, p) # [1.49683032, 1.5       ]

# test local minima/maxima
uvExtrema = rationalSurfaceExtremaParam(surf, p, maxSearch=True, localExtrema=True)
MinPoints = [surf.derivatives(uvl[0], uvl[1], order=0)[0][0] for uvl in uvExtrema]

# test global minima/maxima
uvExtrema = rationalSurfaceExtremaParam(surf, p, maxSearch=True, localExtrema=False)
MinPoint = surf.derivatives(uvExtrema[0], uvExtrema[1], order=0)[0][0]

mmax = [np.linalg.norm(np.array(m)- p) for m in MinPoints]
MaxPointL = [MinPoints[mmax.index(max(mmax))], p.tolist()]

MinPointsL = sum([[m, p.tolist()] for m in MinPoints], [])[:-1]

# SearchPoints = [surf.derivatives(uvl[0], uvl[1], order=0)[0][0] for uvl in uvCoarse]
# SearchPointsL = sum([[m, p.tolist()] for m in SearchPoints], [])[:-1]



colors = [(1,0,0,c) for c in np.linspace(0,1,100)]
cmapred = mcolors.LinearSegmentedColormap.from_list('mycmap', colors, N=5)
colors = [(0,0,1,c) for c in np.linspace(0,1,100)]

cmapblue = mcolors.LinearSegmentedColormap.from_list('mycmap', colors, N=5)

#==============================================

# Plot the control points grid and the evaluated surface

# Create a visualization configuration instance with no legend, no axes and set the resolution to 120 dpi
vis_config = VisMPL.VisConfig(alpha=0.4)

surf.vis = VisMPL.VisSurface(vis_config)
surf.vis.mconf['others']='points'
surf.vis.mconf['alpha']=0.1

pMinPoints = dict(points=MinPointsL, name="pMinPoints", color="green", size=1)
#pSearchPoints = dict(points=SearchPointsL, name="pSearchPoints", color="black", size=1)
pMaxPoints = dict(points=MaxPointL, name="pSearchPoints", color="red", size=1)

#surf.render(extras=[pMinPoints, pSearchPoints], colormap=cm.cool)
surf.render(extras=[pMinPoints, pMaxPoints], colormap=cm.cool)

pass

