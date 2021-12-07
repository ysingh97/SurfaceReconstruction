# Surface Reconstruction

Performing isosurface extraction using the marching cubes algorithm, working with Professor Greg Turk. Given the position of numerous particles, I use marching cubes to extract a surface. I experiment with a variety of kernel functions to generate the scalar field required for marching cubes, including gaussian, quartic, cauchy, and metaballs functions. Ultimately I decide to use metaballs.

I was given data reflecting the positions of particles being moved by a control policy developed through reinforcement learning. Each row of the data corresponded to a particular frame of the animation. An animation with just the particles can be viewed below:

https://user-images.githubusercontent.com/18231852/145081455-32aad24b-e8fc-486a-891c-9ad5794337d2.mp4

For each frame and its corresponding point cloud, I compute a scalar value for every individual particle. For a particular particle P, I do this by calculating the value of the metaball kernel function for each pair of particles that includes P and summing the results. I do this for every particle to generate the scalar field that is used in the marching cubes algorithm to generate an isosurface. I also take the gradient of the metaball function in order to calculate the normals for the resulting polygon mesh. For each row, or frame, of particle data, I generate a .obj polygon mesh with vertexes, normals, and faces. I then read these .obj files using OpenGL to generate a visualization for each frame. Below is the resulting animation.

https://user-images.githubusercontent.com/18231852/145081925-b71f9d56-6db9-4f4e-8689-36fda1ecd70f.mp4

Compile the mesh generator with g++ -std=c++11 cubes.cpp implicitFunctions.cpp utils.cpp -o cubes -lGL -lGLU -lglut -lfreeimageplus.
Compile the mesh reader with g++ -std=c++11 parser.cpp -o parser -lGL -lGLU -lglut -lfreeimageplus

Run ./cubes to generate .obj polygon meshes for each frame.
Run ./parser to read the meshes and generate an animation.

Uses OpenGL 3.0
