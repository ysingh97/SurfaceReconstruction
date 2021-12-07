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
#include <iterator>
#include <istream>
#include <string>
#include <FreeImagePlus.h>
 
/* Global variables */

unsigned int window_width = 1000, window_height = 700;

char title[] = "Marching Cubes";

int refreshMills = 15;
int frame = 0;

std::string frameDirectory = "./videos/parsed/fullVideo20/";
std::string objDirectory = "./obj/fullVideo20/";
std::vector<std::vector<double>> vertexes;
std::vector<std::vector<double>> normals;
std::vector<std::vector<int>> faces;
std::vector<std::vector<double>> goalsPositions;
std::vector<std::vector<std::vector<double>>> barFaces;
std::vector<double> barPosition;
std::vector<double> barRotation;
std::vector<double> barExtents;


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

bool parseFile(int frame) {
    std::string filename = objDirectory + "Frame_" + std::to_string(frame) + ".obj";
    std::ifstream data(filename);
	vertexes.clear();
    normals.clear();
    faces.clear();
    barPosition.clear();
    barRotation.clear();
    barExtents.clear();
    std::string line;
    if(data.good()) {
        int i = 0;
        while (std::getline(data, line)) {
            std::istringstream splitStream(line);
            // std::cout << "After splitting" << "\n";
            std::vector<std::string> words((std::istream_iterator<std::string>(splitStream)), std::istream_iterator<std::string>());
            // std::cout << line << "\n";
            if (words[0].compare("#") == 0 && words[1].compare("Goal") == 0) {
                std::vector<double> goalPosition;
                for(int i = 2; i < 5; i++) {
                    goalPosition.push_back(std::stod(words[i]));
                }
                goalsPositions.push_back(goalPosition);
            }
            if (words[0].compare("#") == 0 && words[1].compare("Bar_Position:") == 0) {
                barPosition.push_back(std::stod(words[2]));
                barPosition.push_back(std::stod(words[3]));
                barPosition.push_back(std::stod(words[4]));
            }
            if (words[0].compare("#") == 0 && words[1].compare("Bar_Rotation:") == 0) {
                barRotation.push_back(std::stod(words[2]));
                barRotation.push_back(std::stod(words[3]));
            }
            if (words[0].compare("#") == 0 && words[1].compare("Bar_Extents:") == 0) {
                barExtents.push_back(std::stod(words[2]));
                barExtents.push_back(std::stod(words[3]));
                barExtents.push_back(std::stod(words[4]));
            }
            // std::cout << "after comments: " << "\n";
            if (words[0].compare("v") == 0) {
                // std::cout << "v\n";
                std::vector<double> vert;
                for(int i = 1; i <= 3; i++) {
                    // std::cout << words[i] << "\n";
                    vert.push_back(std::stod(words[i]));
                }
                vertexes.push_back(vert);
            }
            if (words[0].compare("vn") == 0) {
                // std::cout << "vn\n";
                std::vector<double> norm;
                for(int j = 1; j <= 3; j++) {
                    // std::cout << words[j] << "\n";
                    norm.push_back(std::stod(words[j]));
                }
                normals.push_back(norm);
            }
            if (words[0].compare("f") == 0) {
                // std::cout << "f\n";
                std::vector<int> face;
                for(int k = 1; k <= 3; k++) {
                    std::size_t slashIndex = words[k].find("/");
                    face.push_back(std::stoi(words[k].substr(0, slashIndex)));
                }
                faces.push_back(face);
            }
            i++;
	    }
	    data.close();
        return true;
    }
    return false;	
}
 
std::string normalToString(double normal[]) {
	std::string normalString;
	normalString += std::to_string(normal[0]) + " " + std::to_string(normal[1]) + " " + std::to_string(normal[2]);
	return normalString;
}

void drawGoals() {
    for (int i = 0; i < goalsPositions.size(); i++) {
            glPushMatrix();
            glTranslatef(goalsPositions[i][0], goalsPositions[i][1], goalsPositions[i][2]);
            glutSolidSphere(.5, 6, 6);
            glPopMatrix();
    }
}

void drawPolygons() {
    GLfloat mat_ambient[] = { .2, .2, .2, .2};
    // GLfloat mat_diffuse[] = { 1.0, 1.0, 1.0, 1.0};
    GLfloat mat_specular[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat mat_shininess[] = { 0.0 };
    glPushMatrix();
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
	// glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

    for (int i = 0; i < faces.size(); i++) {
        int p0Index = faces[i][0];
        int p1Index = faces[i][1];
        int p2Index = faces[i][2];
        double p0Normal[] = {normals[p0Index][0], normals[p0Index][1], normals[p0Index][2]};
        double p1Normal[] = {normals[p1Index][0], normals[p1Index][1], normals[p1Index][2]};
        double p2Normal[] = {normals[p2Index][0], normals[p2Index][1], normals[p2Index][2]};
        
        glBegin(GL_TRIANGLES);
            glColor3f(1.0f, 0.0f, 0.0f);
            glNormal3dv(p0Normal);
            glVertex3f(vertexes[p0Index][0], vertexes[p0Index][1], vertexes[p0Index][2]);
            glNormal3dv(p1Normal);
            glVertex3f(vertexes[p1Index][0], vertexes[p1Index][1], vertexes[p1Index][2]);
            glNormal3dv(p2Normal);
            glVertex3f(vertexes[p2Index][0], vertexes[p2Index][1], vertexes[p2Index][2]);
        glEnd();
    }
    glPopMatrix();
}

void drawBar() {
    double theta = std::atan2(barRotation[1], barRotation[0]);
    if (frame == 20)
        std::cout << "Angle: " << theta << "\n";

	double minX = barPosition[0] - barExtents[0];
	double maxX = barPosition[0] + barExtents[0];
	double minY = barPosition[1] - barExtents[1];
	double maxY = barPosition[1] + barExtents[1];
	double minZ = barPosition[2] - barExtents[2];
	double maxZ = barPosition[2] + barExtents[2];

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

/* Initialize OpenGL Graphics */
void initGL() {
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Set background color to black and opaque
	glClearDepth(1.0f);                   // Set background depth to farthest
	glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
	glDepthFunc(GL_LEQUAL);    // Set the type of depth-test
	glShadeModel(GL_SMOOTH);   // Enable flat shading
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Nice perspective corrections

	// shading
	GLfloat light_position[] = { 1.0, 1.0, 10.0, 0.0 };

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
}
 
/* Handler for window-repaint event. Called back when the window first appears and
	whenever the window needs to be re-painted. */
void display() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear color and depth buffers
	glMatrixMode(GL_MODELVIEW);     // To operate on model-view matrix
 
	glLoadIdentity();                  // Reset the model-view matrix
	glTranslatef(-1.5f, 0.0f, -15.0f);  // Move left and into the screen
	glRotatef(90, 1, 0, 0);

    bool fileParsed = parseFile(frame);
    if (fileParsed) {
        // printf("file parsed");
        drawPolygons();
        drawBar();
        drawGoals();
        writeToImage(frame);
    }
	 
	glutSwapBuffers();  // Swap the front and back frame buffers (double buffering)
	
	frame++;
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
	initGL();                       // Our own OpenGL initialization
	glutTimerFunc(0, timer, 0); 
	glutMainLoop();                 // Enter the infinite event-processing loop
	return 0;
}