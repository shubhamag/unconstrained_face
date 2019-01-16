import math
import scipy as sp
import scipy.linalg
import scipy.sparse
import scipy.io

class Mesh(object):
  def __init__(self):
    self.halfedges = []
    self.vertices = []
    self.edges = []
    self.faces = []
    self.boundaries = []

    self.positions = []
    self.normals = []
    self.indices = []

  def clearTmpData(self):
    self.positions = []
    self.normals = []
    self.indices = []

  def loadObjFile(self, filename):
    self.clearTmpData()
    fi = open(filename)
    if not self.readObjData(fi):
      return False
    if not self.buildMesh():
      return False
    self.indexElements()
    return True

  def loadMatFile(self, filename):
    self.clearTmpData()
    if not self.readMatData(filename):
      return False
    if not self.buildMesh():
      return False
    self.indexElements()

  def indexElements(self):
    for i, vertex in enumerate(self.vertices):
      vertex.index = i
    for i, edge in enumerate(self.edges):
      edge.index = i
      edge.computeEdgeLength()
    for i, face in enumerate(self.faces):
      face.index = i

  def computation(self, filename):
    data = []
    row = []
    col = []

    # Border Laplacian
    dataB = []
    rowB = []
    colB = []

    num_vertices = len(self.vertices)
    kappa = sp.zeros(num_vertices, float)
    border = sp.zeros(num_vertices, sp.uint8)

    for he in self.halfedges:
      border[he.vertex.index] = he.onBoundary
      coeff = (he.cotan() + he.flip.cotan()) / 2.0
      col.append(he.vertex.index)
      row.append(he.flip.vertex.index)
      data.append(-coeff)  # L_i,j = -coeff
      col.append(he.vertex.index)
      row.append(he.vertex.index)  # L_i,i = coeff
      data.append(coeff)

      if he.onBoundary:
        val = 0.5 / he.edge.e_length  # put half of the size on each
        colB.append(he.vertex.index)
        rowB.append(he.vertex.index)
        dataB.append(val)  # L_i,i = 1/e_ij
        colB.append(he.vertex.index)
        rowB.append(he.flip.vertex.index)
        dataB.append(-val)  # L_i,j = 1/e_ij

        colB.append(he.flip.vertex.index)  # he on border on goes in 1 direction, so we need to put this on both
        rowB.append(he.vertex.index)
        dataB.append(-val)  # L_i,j = 1/e_ij
        colB.append(he.flip.vertex.index)
        rowB.append(he.flip.vertex.index)
        dataB.append(val)  # L_i,j = 1/e_ij

      # H = 1/4 * (cot(alpha) + cot(beta)) * e_ij dot (n_j - n_i)
      kappa[he.vertex.index] += coeff * sp.dot(he.edgeVector(), (he.flip.vertex.normal - he.vertex.normal)) / 2.0

    mass = sp.zeros(num_vertices, float)
    for face in self.faces:
      # Place area of face equally in all 3 vertices
      a = face.area()
      mass[face.he.vertex.index] += a / 3.0
      mass[face.he.next.vertex.index] += a / 3.0
      mass[face.he.next.next.vertex.index] += a / 3.0

    laplacian = sp.sparse.coo_matrix((data, (row, col)))
    laplacian_border = sp.sparse.coo_matrix((dataB, (rowB, colB)), shape=(num_vertices, num_vertices))
    sp.io.savemat(filename, {
          'L': laplacian,
          'Lb': laplacian_border,
          'kappa': kappa,
          'M': mass,
          'isB': border,
        })

  # ================================================
  # Mesh creation functions
  # ================================================
  def readObjData(self, fi):
    for line in fi:
      tokens = line.split()

      if len(tokens) == 0:
        continue  # blank line
      elif tokens[0] == 'v':
        self.addPosition(tokens[1:])  # vertex
      elif tokens[0] == 'vn':
        self.addNormal(tokens[1:])  # vertex normal
      elif tokens[0] == 'f':
        self.addObjFace(tokens[1:])  # face
        continue  # face
      else:
        print 'Error: does not appear to be a valid Wavefront OBJ file!'
        print '(Offending line: %s)' % line
        return False
    return True

  def readMatData(self, filename):
    try:
      matdata = sp.io.loadmat(filename)
      X = matdata['X']  # vertex positions
      N = matdata['N']  # vertex normals
      T = matdata['T']  # Triangles
      num_vertices = len(X[1,:])
      for i in xrange(num_vertices):
        self.addPosition(X[:,i])
        self.addNormal(N[:,i])
      num_triangles = len(T[:,1])
      for i in xrange(num_triangles):
        self.addFace(T[i,:])
    except KeyError:
      print 'Error: provide X, N, and T'
      return False
    return True


  def buildMesh(self):
    """
    Build the mesh from the data previously read.
    """
    # Data structures used for reference while building mesh
    edgeCount = {}  # Dict from (int, int) -> int
    existingHalfEdges = {}  # Dict from (int, int) -> HalfEdge
    hasFlipEdge = {}  # Dict from HalfEdge -> bool

    # Clear current mesh
    self.halfedges = []
    self.edges = []
    self.faces = []
    self.boundaries = []

    # Allocate a vertex for each position
    self.vertices = [Vertex(position, normal) for (position, normal) in zip(self.positions, self.normals)]

    # Insert each face into the mesh
    degenerateFaces = False
    for faceIndices in self.indices:
      N = len(faceIndices)

      if N < 3:
        print 'Error: face %d is degenerate (fewer than three vertices)!' % (len(self.faces) + 1)
        degenerateFaces = True
        continue

      # Create a new face
      newFace = Face()
      self.faces.append(newFace)

      # Create a new half edge for each edge of the face
      hes = [HalfEdge() for i in xrange(N)]

      for i, he in enumerate(hes):
        v_i = faceIndices[i].position
        v_j = faceIndices[(i + 1) % N].position

        # Point HalfEdge to next edge on face
        self.halfedges.append(he)
        he.next = hes[(i + 1) % N]

        # Associate vertex and half edge to each other
        he.vertex = self.vertices[v_i]
        self.vertices[v_i].he = he

        # Keep track of which edges have flip edges
        hasFlipEdge[he] = False

        # Associate face and half edge to each other
        he.face = newFace
        newFace.he = he

        # Create edge if it doesn't exist
        # If we've created an edge between a and b it is the flip edge
        if v_i > v_j:
          v_i, v_j = v_j, v_i
        if (v_i, v_j) in existingHalfEdges:
          he.flip = existingHalfEdges[(v_i, v_j)]
          he.flip.flip = he
          he.edge = he.flip.edge
          hasFlipEdge[he] = True
          hasFlipEdge[he.flip] = True
        else:
          edge = Edge()
          self.edges.append(edge)
          he.edge = edge
          edge.he = he
          edgeCount[(v_i, v_j)] = 0

        # Record the face that we created the half edge
        existingHalfEdges[(v_i, v_j)] = he

        # Check for nonmanifold edges.
        edgeCount[(v_i, v_j)] += 1
        if edgeCount[(v_i, v_j)] > 2:
          print "Error: edge (%d, %d) is nonmanifold (more than two faces sharing a single edge)." % v_i, v_j
          return False

    # Give up now if there are degenerate faces
    if degenerateFaces:
      return False

    # Insert extra faces for each boundary cycle
    boundaryCycle = []
    for he in self.halfedges:
      # If we find a half edge without a flip edge, it is a boundary
      if not hasFlipEdge[he]:
        newBoundary = Face()
        self.boundaries.append(newBoundary)

        # Walk around this boundary
        startHe = he
        while True:
          # Create a new half edge on the boundary
          newHe = HalfEdge()
          self.halfedges.append(newHe)

          # Only half edges on the boundary face should be marked
          newHe.onBoundary = True

          # Link the current half edge to its new flip edge
          he.flip = newHe

          # Find the next half edge along the boundary
          nextHe = he.next
          while hasFlipEdge[nextHe]:
            nextHe = nextHe.flip.next

          # Set attributes for the boundary edge
          newHe.flip = he
          newHe.vertex = nextHe.vertex
          newHe.edge = he.edge
          newHe.face = newBoundary

          # Point the new face to this half edge
          newBoundary.he = newHe

          boundaryCycle.append(newHe)

          # walk along the cycle until we come to the beginning
          he = nextHe
          if he == startHe:
            break

        # link together the new cycle
        N = len(boundaryCycle)
        for i, he in enumerate(boundaryCycle):
          he.next = boundaryCycle[(i + N - 1) % N]
          hasFlipEdge[he] = True
          hasFlipEdge[he.flip] = True

    return True


  def addPosition(self, tokens):
    self.positions.append(sp.array([
        float(tokens[0]),
        float(tokens[1]),
        float(tokens[2])
      ]))

  def addNormal(self, tokens):
    self.normals.append(sp.array([
        float(tokens[0]),
        float(tokens[1]),
        float(tokens[2])
      ]))

  def addFace(self, tokens):
    faceIndices = []
    for token in tokens:
      faceIndices.append(Index(token, -1, token))
    self.indices.append(faceIndices)

  def addObjFace(self, tokens):
    """
    Face in format with n vertices
    p/[t]/[n] ...

    p: index position
    t: texture position
    n: normal position
    """

    faceIndices = []
    for token in tokens:
      indices = token.split('/')
      faceIndices.append(Index(*indices))
    self.indices.append(faceIndices)



class Vertex(object):
  index = 0
  he = None

  def __init__(self, position, normal):
    self.position = position
    self.normal = normal


class HalfEdge(object):
  next = None  # next halfedge around the current face
  flip = None  # opposite halfedge of this edge
  vertex = None  # vertex at the tail of this halfedge h_ij the tail is v_i
  edge = None
  face = None
  onBoundary = False

  def cotan(self):
    if self.onBoundary:
      return 0.0
    p0 = self.next.next.vertex.position
    p1 = self.vertex.position
    p2 = self.next.vertex.position

    u = p1 - p0
    v = p2 - p0

    val = sp.dot(u, v) / sp.linalg.norm(sp.cross(u, v))

    # Correct numerical issues
    if math.fabs(val) > 1e6:
      return 1e6
    return val

  def edgeVector(self):
    pi = self.vertex.position
    pj = self.flip.vertex.position

    return pj - pi

class Edge(object):
  he = None
  index = 0
  e_length = 0.0

  def computeEdgeLength(self):
    p_i = self.he.vertex.position
    p_j = self.he.flip.vertex.position
    self.e_length = sp.linalg.norm(p_i - p_j)


class Face(object):
  he = None
  index = 0

  def area(self):
    """
    Area of face, assuming it is a triangle.
    """
    a = self.he.vertex.position
    b = self.he.next.vertex.position
    c = self.he.next.next.vertex.position
    ab = b - a
    ac = c - a
    return sp.linalg.norm(sp.cross(ab,ac)) / 2.0


class Index(object):
  def __init__(self, p, t, n):
    # Decrement since OBJ files are 1-indexed
    self.position = int(p) - 1
    self.normal = int(n) - 1
