import mesh_lib

mesh = mesh_lib.Mesh();
mesh.loadMatFile('Python/meshin');
mesh.computation('Python/meshout');
