/*
 * cubes.cpp: Contains vertex interpolation code and marching cubes code
 */
//include <windows.h>  // for MS Windows
#include <GL/glut.h>  // GLUT, include glu.h and gl.h
#include <cmath>
#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "implicitFunctions.h"
#include "utils.h"
#include "cubes.h"
#include <FreeImagePlus.h>
#include <map>

unsigned int window_width = 1000, window_height = 700;

char title[] = "Marching Cubes";

typedef double (*ImplicitFunctions)(XYZ point, XYZ poi);
ImplicitFunctions functions[] = {sphere, gaussian, cauchy, inverse, inverseSquared, metaballs, modifiedMetaballs, softObjects, quartic};

int function = 6;

bool showTriangles = true;
bool showPoints = false;
bool rotate = false;
bool simpleCloud = true;

double cubeSize = .1;
double regIsovalue = .8;
double smallCloudIsoValue = .05;

int frame = 400;

int numVertexes = 0;

double t = 0;
int refreshMills = 15;

double coord = .8;
XYZ a(coord, coord, coord, 0, 0, 0);

int numPolygons = 0;

//data about scene
int numGoals = 0;
std::vector<XYZ> goalsPositions;
double bar_half_extent_x = 0;
double bar_half_extent_y = 0;
double bar_half_extent_z = 0;

int minX = (int)(-2/cubeSize);
int maxX = (int)(2/cubeSize);
int minY = (int)(-1/cubeSize);
int maxY = (int)(1/cubeSize);
int minZ = (int)(-3/cubeSize);
int maxZ = (int)(3/cubeSize);

double functionWidth = 400;

// XYZ b(-coord, -coord, -coord);

std::vector<XYZ> pointCloud;
std::vector<int> barStatus;
std::vector<double> barPosition;
std::vector<double> barRotation;

std::string frameDirectory = "./videos/correct/fullVideo17/";
std::string objDirectory = "./obj/fullVideo17/";

std::vector<std::string> lines;
std::vector<std::string> sceneData;

std::map<std::string, int> vertexMap;
std::map<std::vector<int>, int> ijkMap;
std::map<std::vector<int>, std::vector<XYZ>> subdivisionGrid;

std::vector<XYZ> vertexList;
std::vector<std::vector<int>> faceList;

void writeToImage(int frame) {
	std::string fileName = "Frame_" + std::to_string(frame);

	BYTE* pixels = new BYTE[3 * window_width * window_height];
	glReadPixels(0, 0, window_width, window_height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	// Convert to FreeImage format & save to file
	FIBITMAP* image = FreeImage_ConvertFromRawBits(pixels, window_width, window_height, 3 * window_width, 24, 0x0000FF, 0xFF0000, 0x00FF00, false);
	FreeImage_Save(FIF_BMP, image, (frameDirectory + fileName).c_str(), 0);

	// Free resources
	FreeImage_Unload(image);
	delete [] pixels;
}

void writeMeshToFile(int frame) {
	std::ofstream fileStream;
	std::string filename = objDirectory + "Frame_" + std::to_string(frame) + ".obj";
	fileStream.open(filename);
	if (fileStream.is_open()) {

		fileStream << "# Bar_Position: " << barPosition[0] << " " << barPosition[1] << " " << barPosition[2] << "\n";
		fileStream << "# Bar_Rotation: " << barRotation[0] << " " << barRotation[2] << "\n";
		fileStream << "# Bar_Extents: " << bar_half_extent_x << " " << bar_half_extent_y << " " << bar_half_extent_z << "\n";

		for (int i = 0; i < goalsPositions.size(); i++) {
			fileStream << "# Goal " << goalsPositions[i].toString() << "\n";
		}
 		 
		for (int i = 0; i < vertexList.size(); i++) {
			XYZ p = vertexList[i];
			fileStream << "v " << p.x << " " << p.y << " " << p.z << "\n";
		}
		// fileStream << "\n";
		for (int i = 0; i < vertexList.size(); i++) {
			XYZ p = vertexList[i];
			fileStream << "vn " << p.normal[0] << " " << p.normal[1] << " " << p.normal[2] << "\n";
		}
		// fileStream << "\n";
		for (int i = 0; i < faceList.size(); i++) {
			std::vector<int> face = faceList[i];
			fileStream << "f " << face[0] << "//" << face[0] << " " << face[1] << "//" << face[1] << " " << face[2] << "//" << face[2] << "\n";
		}
		// fileStream << "\n";
		fileStream.close();
	}
}

void debug(std::string message) {
	std::ofstream debugging;
	debugging.open("debugging.txt", std::ofstream::app);
	if (debugging.is_open()) {
		debugging << message + "\n";
		debugging.close();
	}
}

void readData() {
	std::ifstream data("sample_data.csv");
	std::string line;
	if (data.good()) {
		while (std::getline(data, line)) {
			lines.push_back(line);
		}
		data.close();
	}
	
}

void parseFullLines() {
	// std::cout << "reading lines\n";
	for (int i = 0; i < sceneData.size(); i++) {
		std::string line = sceneData[i];
		std::istringstream s(line);
		std::string field;
		int count = 0;
		std::vector<double> goalPosition;
		// std::cout << "line: " << line << "\n";
		switch (i) {
			case 0:
				count = 0;
				while(getline(s, field, ',')) {
					if (count == 1) {
						numGoals = std::stoi(field);
					}
					count++;
				}
				break;
			case 1:
				count = 0;
				while(getline(s, field, ',')) {
					// std::cout << "data: " << field << "\n";
					if (count > 0) {
						goalPosition.push_back(std::stod(field));
						if (count % 2 == 0) {
							XYZ point;
							point.x = goalPosition[0];
							point.y = 0;
							point.z = goalPosition[1];
							goalsPositions.push_back(point);
							goalPosition.clear();
						}
					}
					count++;
				}
				break;
			case 2:
				count = 0;
				while(getline(s, field, ',')) {
					if (count == 1) {
						bar_half_extent_x = std::stod(field);
					}
					count++;
				}
				break;
			case 3:
				count = 0;
				while(getline(s, field, ',')) {
					if (count == 1) {
						bar_half_extent_y = std::stod(field);
					}
					count++;
				}
				break;
			case 4:
				count = 0;
				while(getline(s, field, ',')) {
					if (count == 1) {
						bar_half_extent_z = std::stod(field);
					}
					count++;
				}
				break;		
		}
	}
	for(int i = 0; i < goalsPositions.size(); i++) {
		std::cout << "goal " << i << ": " << goalsPositions[i].toString() << "\n";
	}
}

std::vector<int> generateSubdivision(XYZ point) {
	double minXPosition = minX * cubeSize;
	double minYPosition = minY * cubeSize;
	double minZPosition = minZ * cubeSize;
	double scaledFunctionWidth = functionWidth/modifier;
	std::vector<int> subdivision = {(int)floor((point.x - minXPosition)/scaledFunctionWidth), (int)floor((point.y - minYPosition)/scaledFunctionWidth), (int)floor((point.z - minZPosition)/scaledFunctionWidth)};
	return subdivision;
}

void generateFrameData(int frame) {
	pointCloud.clear();
	barPosition.clear();
	barRotation.clear();
	barStatus.clear();

	std::string line = lines[frame];
	std::istringstream s(line);
	std::string field;
	int count = 0;
	XYZ point;
	// std::cout << "Full frame data\n";
	while(getline(s, field , ',')) {			
		double data = std::stod(field);
		if (count < 3) {
			// std::cout << "bar position\n";
			barPosition.push_back(data);
		} else if (count >= 3 && count <= 5) {
			// std::cout << "bar rotation\n";
			barRotation.push_back(data);
		} else if (count == 6) {
			// std::cout << "barStatus\n";
			data = std::stoi(field);
			barStatus.push_back(data);
		} else {
			// std::cout << "points: " << field << "\n";
			switch (count % 3) {
				case 1:
					point.x = data;
					break;
				case 2:
					point.y = data;
					break;
				case 0:
					point.z = data;
					std::vector<int> subdivision = generateSubdivision(point);
					if (subdivisionGrid.count(subdivision) <= 0) {
						std::vector<XYZ> gridPoints = {point};
						subdivisionGrid[subdivision] = gridPoints;
					} else {
						std::vector<XYZ> gridPoints = subdivisionGrid[subdivision];
						gridPoints.push_back(point);
						subdivisionGrid[subdivision] = gridPoints;
					}
					pointCloud.push_back(point);
					break;
			}
		} 
		count++;
	}
}

void generatePoints(int frame) {
	pointCloud.clear();
	barPosition.clear();
	barRotation.clear();
	std::string line = lines[frame];
	std::istringstream s(line);
	std::string field;
	int count = 0;
	XYZ point;
	while(getline(s, field , ',')) {			
		double data = std::stod(field);
		if (count % 3 == 0) {
			point.x = data;
		} else if (count % 3 == 1) {
			point.y = data;
		} else if (count % 3 == 2) {
			point.z = data;
			pointCloud.push_back(point);
		}
		count++;
	}
}

void readFullData() {
	std::ifstream data("sample_data_full.csv");
	std::string line;
	int lineNumber = 0;
	if (data.good()) {
		while (std::getline(data, line)) {
			if (lineNumber < 5) {
				sceneData.push_back(line);
			} 
			if (lineNumber > 5) {
				lines.push_back(line);
			}	
			lineNumber++;
		}
		data.close();
	}

	parseFullLines();
}

void generateSmallCloud() {
	pointCloud.resize(0);
	pointCloud.push_back(a);
}

void drawPoints() {
	for (int j = 0; j < pointCloud.size(); j++) {
		XYZ point = pointCloud[j];
		glPointSize(10);
		glBegin(GL_POINTS);
			glColor3f(1.0, 0.0, 0.0);
			glVertex3f(point.x, point.y, point.z);
		glEnd();
	}
}

std::vector<int> findEdgeIJK(XYZ p1, XYZ p2) {
		std::vector<int> ijk;
		bool p1First = true;
		if (p1.ijk[0] > p2.ijk[0]) {
			p1First = false;
		} else if(p1.ijk[0] == p2.ijk[0]) {
				if (p1.ijk[1] > p2.ijk[1]) {
					p1First = false;
				} else if (p1.ijk[1] == p2.ijk[1]) {
						if (p1.ijk[2] > p2.ijk[2]) {
							p1First = false;
						}
				}
		}

		ijk.reserve(ijk.size() + p1.ijk.size() + p2.ijk.size());
		if (p1First) {
				ijk.insert(ijk.end(), p1.ijk.begin(), p1.ijk.end());
				ijk.insert(ijk.end(), p2.ijk.begin(), p2.ijk.end());
		} else {
				ijk.insert(ijk.end(), p2.ijk.begin(), p2.ijk.end());
				ijk.insert(ijk.end(), p1.ijk.begin(), p1.ijk.end());	
		}
		return ijk;
}

/*
	Linearly interpolate the position where an isosurface cuts
	an edge between two vertices, each with their own scalar value
*/
XYZ VertexInterp(double isolevel, XYZ p1, XYZ p2, double valp1, double valp2)
{
	// printf("val 1 val 2: %f, %f\n", valp1, valp2);
	double mu;
	XYZ p;
	std::vector<int> newIJK = findEdgeIJK(p1, p2);
	p.ijk = newIJK;
	
	if (fabs(isolevel-valp1) < 0.00001)
		return(p1);
	if (fabs(isolevel-valp2) < 0.00001)
		return(p2);
	if (fabs(valp1-valp2) < 0.00001) 
		return(p1);
		
	mu = (isolevel - valp1) / (valp2 - valp1);
	p.x = p1.x + mu * (p2.x - p1.x);
	p.y = p1.y + mu * (p2.y - p1.y);
	p.z = p1.z + mu * (p2.z - p1.z);
	
	double distanceInterp = distance(p1, p)/distance(p1, p2);
	p.normal[0] = p1.normal[0] + distanceInterp * (p2.normal[0] - p1.normal[0]);
	p.normal[1] = p1.normal[1] + distanceInterp * (p2.normal[1] - p1.normal[1]);
	p.normal[2] = p1.normal[2] + distanceInterp * (p2.normal[2] - p1.normal[2]);
	p.normal = normalize(p.normal);
	
	return(p);
}

/*
	Given a grid cell and an isolevel, calculate the triangular
	facets required to represent the isosurface through the cell.
	Return the number of triangular facets, the array "triangles"
	will be loaded up with the vertices at most 5 triangular facets.
	0 will be returned if the grid cell is either totally above
	of totally below the isolevel.
*/

int Polygonise(GRIDCELL grid,double isolevel,TRIANGLE *triangles)
{
	// printf("polygonising\n");
	int i,ntriang;
	int cubeindex;
	XYZ vertlist[12];

	/*
		Determine the index into the edge table which
		tells us which vertices are inside of the surface
	*/
	cubeindex = 0;
	if (grid.val[0] < isolevel) cubeindex |= 1;
	if (grid.val[1] < isolevel) cubeindex |= 2;
	if (grid.val[2] < isolevel) cubeindex |= 4;
	if (grid.val[3] < isolevel) cubeindex |= 8;
	if (grid.val[4] < isolevel) cubeindex |= 16;
	if (grid.val[5] < isolevel) cubeindex |= 32;
	if (grid.val[6] < isolevel) cubeindex |= 64;
	if (grid.val[7] < isolevel) cubeindex |= 128;

	/* Cube is entirely in/out of the surface */
	if (edgeTable[cubeindex] == 0)
		// printf("return 0\n");
		return(0);

	// printf("Cubeindex: %d\n", cubeindex);

	/* Find the vertices where the surface intersects the cube */
	if (edgeTable[cubeindex] & 1)
		vertlist[0] =
			VertexInterp(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);
	if (edgeTable[cubeindex] & 2)
		vertlist[1] =
			VertexInterp(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);
	if (edgeTable[cubeindex] & 4)
		vertlist[2] =
			VertexInterp(isolevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);
	if (edgeTable[cubeindex] & 8)
		vertlist[3] =
			VertexInterp(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
	if (edgeTable[cubeindex] & 16)
		vertlist[4] =
			VertexInterp(isolevel,grid.p[4],grid.p[5],grid.val[4],grid.val[5]);
	if (edgeTable[cubeindex] & 32)
		vertlist[5] =
			VertexInterp(isolevel,grid.p[5],grid.p[6],grid.val[5],grid.val[6]);
	if (edgeTable[cubeindex] & 64)
		vertlist[6] =
			VertexInterp(isolevel,grid.p[6],grid.p[7],grid.val[6],grid.val[7]);
	if (edgeTable[cubeindex] & 128)
		vertlist[7] =
			VertexInterp(isolevel,grid.p[7],grid.p[4],grid.val[7],grid.val[4]);
	if (edgeTable[cubeindex] & 256)
		vertlist[8] =
			VertexInterp(isolevel,grid.p[0],grid.p[4],grid.val[0],grid.val[4]);
	if (edgeTable[cubeindex] & 512)
		vertlist[9] =
			VertexInterp(isolevel,grid.p[1],grid.p[5],grid.val[1],grid.val[5]);
	if (edgeTable[cubeindex] & 1024)
		vertlist[10] =
			VertexInterp(isolevel,grid.p[2],grid.p[6],grid.val[2],grid.val[6]);
	if (edgeTable[cubeindex] & 2048)
		vertlist[11] =
			VertexInterp(isolevel,grid.p[3],grid.p[7],grid.val[3],grid.val[7]);

	/* Create the triangle */
	ntriang = 0;
	// printf("before polygonise loop\n");
	for (i=0;triTable[cubeindex][i]!=-1;i+=3) {
		triangles[ntriang].p[0] = vertlist[triTable[cubeindex][i  ]];
		triangles[ntriang].p[1] = vertlist[triTable[cubeindex][i+1]];
		triangles[ntriang].p[2] = vertlist[triTable[cubeindex][i+2]];

		XYZ vertexZero = triangles[ntriang].p[0];
		XYZ vertexOne = triangles[ntriang].p[1];
		XYZ vertexTwo = triangles[ntriang].p[2];

		double p0Normal[] = {vertexZero.normal[0], vertexZero.normal[1], vertexZero.normal[2]};
		double p1Normal[] = {vertexOne.normal[0], vertexOne.normal[1], vertexOne.normal[2]};
		double p2Normal[] = {vertexTwo.normal[0], vertexTwo.normal[1], vertexTwo.normal[2]};

		if (ijkMap.count(vertexZero.ijk) <= 0) {
				ijkMap.insert({vertexZero.ijk, numVertexes});
				vertexList.push_back(vertexZero);
				// normalList.push_back(p0Normal);
				numVertexes++;
				
		}
		if (ijkMap.count(vertexOne.ijk) <= 0) {
				ijkMap.insert({vertexOne.ijk, numVertexes});
				vertexList.push_back(vertexOne);
				// normalList.push_back(p1Normal);
				numVertexes++;
				
				// printf("adding vertextwo to map\n");
		}
		if (ijkMap.count(vertexTwo.ijk) <= 0) {
				ijkMap.insert({vertexTwo.ijk, numVertexes});
				vertexList.push_back(vertexTwo);
				// normalList.push_back(p2Normal);
				numVertexes++;
				// printf("adding vertex3ToMap\n");
		}
		int faceZeroIndex = ijkMap.at(vertexZero.ijk);
		int faceOneIndex = ijkMap.at(vertexOne.ijk);
		int faceTwoIndex = ijkMap.at(vertexTwo.ijk);
		std::vector<int> face;
		face.push_back(faceZeroIndex);
		face.push_back(faceOneIndex);
		face.push_back(faceTwoIndex);
		faceList.push_back(face);

		glBegin(GL_TRIANGLES);
			glNormal3dv(p0Normal);
			glVertex3f(triangles[ntriang].p[0].x, triangles[ntriang].p[0].y, triangles[ntriang].p[0].z);
			glNormal3dv(p1Normal);
			glVertex3f(triangles[ntriang].p[1].x, triangles[ntriang].p[1].y, triangles[ntriang].p[1].z);
			glNormal3dv(p2Normal);
			glVertex3f(triangles[ntriang].p[2].x, triangles[ntriang].p[2].y, triangles[ntriang].p[2].z);
		glEnd();
		ntriang++;
	}
	
	return(ntriang);
}

void generateTriangles() {
	// printf("Begin generating traingles\n");
	std::map<std::vector<int>, double> valueMap;
	std::map<std::vector<int>, std::vector<double>> gradientMap;

	GLfloat mat_ambient[] = { .2, .2, .2, .2};
	// GLfloat mat_diffuse[] = { 1.0, 1.0, 1.0, 1.0};
	GLfloat mat_specular[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat mat_shininess[] = { 0.0 };
	glPushMatrix();

		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
		// glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

	for(int x = minX; x < maxX; x++) {
		for(int y = minY; y < maxY; y++) {
			for(int z = minZ; z < maxZ; z++) {
				GRIDCELL cell;
				TRIANGLE triangles[5];
				XYZ a {x*cubeSize, y*cubeSize, z*cubeSize, x, y, z};
				XYZ b {(x + 1)*cubeSize, y*cubeSize, z*cubeSize, x + 1, y, z};
				XYZ c {(x + 1)*cubeSize, y*cubeSize, (z + 1)*cubeSize, x + 1, y, z + 1};
				XYZ d {x*cubeSize, y*cubeSize, (z + 1)*cubeSize, x, y, z + 1};
				XYZ e {x*cubeSize, (y + 1)*cubeSize, z*cubeSize, x, y + 1, z};
				XYZ f {(x + 1)*cubeSize, (y + 1)*cubeSize, z*cubeSize, x + 1, y + 1, z};
				XYZ g {(x + 1)*cubeSize, (y + 1)*cubeSize, (z + 1)*cubeSize, x + 1, y + 1, z + 1};	
				XYZ h {x*cubeSize, (y + 1)*cubeSize, (z + 1)*cubeSize, x, y + 1, z + 1};
				// printf("after creating gridcell points\n");
				std::array<XYZ, 8> points = {a, b, c, d, e, f, g, h};
				std::array<double, 8> values;
				std::array<std::vector<double>, 8> gradients;
				// printf("a: %f, %f, %f\n", a.x, a.y, a.z);
				/* 	0 - basic sphere
					1 - gaussian
					2 - cauchy
					3 - inverse
					4 - inverse squared
					5 - metaballs
					6 - soft objects
					7 - quartic
				*/
				//speed up here, hash
				for (int i = 0; i < 8; i++) {
					double value = 0;
					XYZ gridPoint = points[i];
					std::vector<int> ijk = gridPoint.ijk;
					if (valueMap.count(ijk) <= 0) {
						std::vector<double> totalGradient = {0, 0, 0};

						std::vector<int> subdivision = generateSubdivision(gridPoint);
						for (int x = -1; x <= 1; x++) {
							for (int y = -1; y <= 1; y++) {
								for (int z = -1; z <= 1; z++) {
									std::vector<int> gridNeighbor = {subdivision[0] + x, subdivision[1] + y, subdivision[2] + z};
									// std::cout << "grid neighbor: ";
									// std::cout << gridNeighbor[0] << " " << gridNeighbor[1] << " " << gridNeighbor[2] << "\n";
									if (subdivisionGrid.count(gridNeighbor)) {
										std::vector<XYZ> neighborPoints = subdivisionGrid[gridNeighbor];
										// std::cout << "Points: " << points.size() << "\n";
										for (int i = 0; i < neighborPoints.size(); i++) {
											// std::cout << "Point: " << neighborPoints[i].toString() << "\n";
											double increase = functions[function](gridPoint, neighborPoints[i]);
											// std::cout << increase << "\n";
											value += increase;
											if (increase > 0) {
												std::vector<double> gradient = modifiedMetaballGradient(gridPoint, neighborPoints[i]);
												// printf("gradient: %f, %f, %f\n", gradient[0], gradient[1], gradient[2]);
												totalGradient[0] -= gradient[0];
												totalGradient[1] -= gradient[1];
												totalGradient[2] -= gradient[2];
											}	
										}
									}
								}
							}
						}
						// std::cout << "value: " << value << "\n";
						// for (int j = 0; j < pointCloud.size(); j++) {
						// 	double increase = functions[function](gridPoint, pointCloud[j]);
						// 	value += increase;
						// 	if (increase > 0) {
						// 		std::vector<double> gradient = modifiedMetaballGradient(gridPoint, pointCloud[j]);
						// 		// printf("gradient: %f, %f, %f\n", gradient[0], gradient[1], gradient[2]);
						// 		totalGradient[0] -= gradient[0];
						// 		totalGradient[1] -= gradient[1];
						// 		totalGradient[2] -= gradient[2];
						// 	}	
						// }
						values[i] = value;
						gradients[i] = normalize(totalGradient);
						if (value > 0) {
							// printf("IJK: %d, %d, %d, Gradient: {%f, %f, %f}, Value: %f\n", ijk[0], ijk[1], ijk[2], gradients[i][0], gradients[i][1], gradients[i][2], value);
						}
						
						// printf("%f, %f, %f\n", gradients[i][0], gradients[i][1], gradients[i][2]);
						valueMap.insert({ijk, value});
						gradientMap.insert({ijk, gradients[i]});
						points[i].normal = gradients[i];
						// printf("not in map: \n");
						// }

					} else {
						// printf("already in map\n");
						values[i] = valueMap[ijk];
						gradients[i] = gradientMap[ijk];
						points[i].normal = gradientMap[ijk];
					}
					
					// bool normalCalculated = gradientMap.count(ijk) q== 1;
					// printf("after checkign if normal has already been taken\n");
					
				}
				// printf("before assigning gridcell\n");
				cell.p = points;
				cell.val = values;
				cell.normals = gradients;
				double isovalue = simpleCloud ? smallCloudIsoValue : regIsovalue;
				Polygonise(cell, isovalue, triangles);
			}
		}
	}
	glPopMatrix();
}
 
/* Initialize OpenGL Graphics */
void initGL() {
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Set background color to black and opaque
	glClearDepth(1.0f);                   // Set background depth to farthest
	glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
	glDepthFunc(GL_LEQUAL);    // Set the type of depth-test
	glShadeModel(GL_SMOOTH);   // Enable flat shading
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Nice perspective corrections

	GLfloat light_position[] = { 1.0, 1.0, 10.0, 0.0 };

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
}

void drawBar(int frame) {

	double theta = std::atan2(barRotation[2], barRotation[0]);

	if (frame == 20) {
		std::cout << "Angle: " << theta << "\n";
	}

	double minX = barPosition[0] - bar_half_extent_x;
	double maxX = barPosition[0] + bar_half_extent_x;
	double minY = barPosition[1] - bar_half_extent_y;
	double maxY = barPosition[1] + bar_half_extent_y;
	double minZ = barPosition[2] - bar_half_extent_z;
	double maxZ = barPosition[2] + bar_half_extent_z;

	GLfloat mat_ambient[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat mat_diffuse[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat mat_specular[] = { 1.0, 0.0, 0.0, 1.0 };
	GLfloat mat_shininess[] = { 50.0 };  
	
	glPushMatrix();

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
	glTranslatef(barPosition[0], barPosition[1], barPosition[2]);
	glRotatef(theta * (180/M_PI), 0, 1, 0);
	glTranslatef(-barPosition[0], -barPosition[1], -barPosition[2]);

	glBegin(GL_QUADS);
		//face 1
		glVertex3f(minX, maxY, maxZ);
		glVertex3f(minX, maxY, minZ);
		glVertex3f(maxX, maxY, minZ);
		glVertex3f(maxX, maxY, maxZ);

		//face 2
		glVertex3f(maxX, maxY, maxZ);
		glVertex3f(maxX, maxY, minZ);
		glVertex3f(maxX, minY, minZ);
		glVertex3f(maxX, minY, maxZ);

		//face 3
		glVertex3f(minX, maxY, maxZ);
		glVertex3f(maxX, maxY, maxZ);
		glVertex3f(maxX, minY, maxZ);
		glVertex3f(minX, minY, maxZ);

		//face 4
		glVertex3f(minX, maxY, minZ);
		glVertex3f(minX, minY, minZ);
		glVertex3f(minX, minY, maxZ);
		glVertex3f(minX, maxY, maxZ);

		//face 5
		glVertex3f(maxX, maxY, minZ);
		glVertex3f(maxX, minY, minZ);
		glVertex3f(minX, minY, minZ);
		glVertex3f(minX, maxY, minZ);

		//face 6
		glVertex3f(minX, minY, minZ);
		glVertex3f(maxX, minY, minZ);
		glVertex3f(maxX, minY, maxZ);
		glVertex3f(minX, minY, maxZ);

	glEnd();
	glPopMatrix();
}

void drawBackground() {
	GLfloat mat_ambient[] = {1.0, 0.0, 0.0, 1.0};
	GLfloat mat_diffuse[] = {1.0, 0.0, 0.0, 1.0};
	GLfloat mat_specular[] = { 1.0, 0.0, 0.0, 1.0 };
	GLfloat mat_shininess[] = { 50.0 };  
	
	glPushMatrix();
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_ambient);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
		glBegin(GL_QUADS);
			glVertex3f( 5.0f,  -3.0f, 5.0f);
			glVertex3f( -5.0f,  -3.0f, 5.0f);
			glVertex3f( -5.0f, -3.0f, -5.0f);
			glVertex3f( 5.0f, -3.0f, -5.0f);
		glEnd();
	glPopMatrix();
}

/* Handler for window-repaint event. Called back when the window first appears and
	whenever the window needs to be re-painted. */
void display() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear color and depth buffers
	glMatrixMode(GL_MODELVIEW);     // To operate on model-view matrix
 
	glLoadIdentity();                  // Reset the model-view matrix
	glTranslatef(-1.5f, 0.0f, -15.0f);  // Move left and into the screen
	glRotatef(90, 1, 0, 0);
	// glRotatef(t, 0, 1, 0);
	numVertexes = 0;
	vertexMap.clear();
	ijkMap.clear();
	vertexList.clear();
	faceList.clear();
	subdivisionGrid.clear();

	if (simpleCloud) {
		generateSmallCloud();
		generateTriangles();
		writeMeshToFile(frame);
	} else {
		printf("%d\n", frame);
		// printf("%f\n", modifier);
		if (frame < lines.size()) {
			// generatePoints(frame);
			generateFrameData(frame);
			// std::cout << "numpoints: " + pointCloud.size() << "\n";
			// std::cout << "generated frame data \n";
			drawBar(frame);
			generateTriangles();
			// drawBackground();
			for (int i = 0; i < goalsPositions.size(); i++) {
				glPushMatrix();
				glTranslatef(goalsPositions[i].x, goalsPositions[i].y, goalsPositions[i].z);
				glutSolidSphere(.5, 6, 6);
				glPopMatrix();
			}
			
			writeMeshToFile(frame);
			writeToImage(frame);
		}
	}

	if (showPoints)
		drawPoints();
	 
	glutSwapBuffers();  // Swap the front and back frame buffers (double buffering)
	if (rotate)
		t += .5;
	
	frame++;
}

void MyKeyboardFunc(unsigned char key, int x, int y) {
	switch (key){
		case 't':
			showTriangles = !showTriangles;
			break;
		case 'r':
			rotate = !rotate;
			break;
		case 'p':
			showPoints = !showPoints;
			break;
	}
}

void SpecialKeyboardFunc(int key, int x, int y) {
	switch (key) {
		case GLUT_KEY_LEFT:
			if (!rotate)
				t -= 1.5;
			break;
		case GLUT_KEY_RIGHT:
			if (!rotate)
				t += 1.5;
			break;
	}
}
 
/* Handler for window re-size event. Called back when the window first appears and
	whenever the window is re-sized with its new width and height */
void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
	// Compute aspect ratio of the new window
	if (height == 0) height = 1;                // To prevent divide by 0
	GLfloat aspect = (GLfloat)width / (GLfloat)height;
 
	// Set the viewport to cover the new window
	glViewport(0, 0, width, height);
 
	// Set the aspect ratio of the clipping volume to match the viewport
	glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
	glLoadIdentity();             // Reset
	// Enable perspective projection with fovy, aspect, zNear and zFar
	gluPerspective(45.0f, aspect, 0.1f, 100.0f);
}

void timer(int value) {
	glutPostRedisplay();      // Post re-paint request to activate display()
	glutTimerFunc(refreshMills, timer, 0); // next timer call milliseconds later
}
 
/* Main function: GLUT runs as a console application starting at main() */
int main(int argc, char** argv) {
	glutInit(&argc, argv);            // Initialize GLUT
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // Enable double buffered mode
	glutInitWindowSize(window_width, window_height);   // Set the window's initial width & height
	glutInitWindowPosition(50, 50); // Position the window's initial top-left corner
	glutCreateWindow(title);          // Create window with the given title
	glutDisplayFunc(display);       // Register callback handler for window re-paint event
	glutReshapeFunc(reshape);       // Register callback handler for window re-size event
	glutKeyboardFunc(MyKeyboardFunc);
	glutSpecialFunc(SpecialKeyboardFunc);
	initGL();                       // Our own OpenGL initialization
	if (!simpleCloud) {
		// readData();
		readFullData();
	}
	glutTimerFunc(0, timer, 0); 
	glutMainLoop();                 // Enter the infinite event-processing loop
	return 0;
}