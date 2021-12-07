# SurfaceReconstruction

Performing isosurface extraction using the marching cubes algorithm. Given the position of numerous particles, I use marching cubes to extract a surface. I experiment with a variety of kernel functions to generate the scalar field required for marching cubes, including gaussian, quartic, cauchy, and metaballs functions. Ultimately I decide to use metaballs.

I was given data reflecting the positions of particles being moved by a control policy developed through reinforcement learning. Each row of the data corresponded to a particular frame of the animation. An animation with just the particles can be viewed below:

https://user-images.githubusercontent.com/18231852/145081455-32aad24b-e8fc-486a-891c-9ad5794337d2.mp4

For each frame, and each corresopnding point cloud, I calculate a scalar value for each particle by summing the 


Compile with g++ -std=c++11 cubes.cpp implicitFunctions.cpp utils.cpp -o cubes -lGL -lGLU -lglut -lfreeimageplus
and g++ -std=c++11 parser.cpp -o parser -lGL -lGLU -lglut -lfreeimageplus

Run with ./cubes to 


Animation after surface reconsruction

https://user-images.githubusercontent.com/18231852/145081925-b71f9d56-6db9-4f4e-8689-36fda1ecd70f.mp4



Uses OpenGL 3.0
