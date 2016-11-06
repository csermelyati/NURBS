#include <GL/glut.h>
#include <math.h>
#include <iostream>
#include <vector>

class myPoint{

public:

  GLdouble x, y;
  myPoint(){
    x = 0.0;
    y = 0.0;
  }
  myPoint(GLdouble X, GLdouble Y){
    x = X;
    y = Y;
  }
  myPoint(const myPoint& a){
    x = a.x;
    y = a.y;
  }

  void set(const myPoint& a){
    x = a.x;
    y = a.y;
  }

private:

};
myPoint operator*(const myPoint& p, GLdouble a){
  myPoint ret(p.x*a,p.y*a);
  return ret;
}
myPoint operator*(GLdouble a, const myPoint& p){
  myPoint ret(p.x*a,p.y*a);
  return ret;
}
myPoint operator+(const myPoint& a, const myPoint& b){
  myPoint ret(a.x+b.x,a.y+b.y);
  return ret;
}
myPoint operator-(const myPoint& a, const myPoint& b){
  myPoint ret(a.x-b.x,a.y-b.y);
  return ret;
}


class myMatrix{
public:
  std::vector<std::vector<GLdouble>> M;
  myMatrix(){
  }
  myMatrix(GLdouble matrix[4][4]){
    M.push_back( std::vector<GLdouble>(&(matrix[0][0]),&(matrix[0][0])+4) );
    M.push_back( std::vector<GLdouble>(&(matrix[1][0]),&(matrix[1][0])+4) );
    M.push_back( std::vector<GLdouble>(&(matrix[2][0]),&(matrix[2][0])+4) );
    M.push_back( std::vector<GLdouble>(&(matrix[3][0]),&(matrix[3][0])+4) );
  }

  std::vector<myPoint> operator*(const std::vector<myPoint> points){
    std::vector<myPoint> ret;
    for (int i = 0; i < 4; i++) {
      myPoint akt;
      for (int j = 0; j < 4; j++) {
        akt = akt + (points[j]*M[i][j]);
      }
      ret.push_back(akt);
    }
    return ret;
  }
private:
};

std::vector<myPoint> operator*(GLdouble M[4][4], const std::vector<myPoint> points){
  std::vector<myPoint> ret;
  for (int i = 0; i < 4; i++) {
    myPoint akt;
    for (int j = 0; j < 4; j++) {
      akt = akt + (points[j]*M[i][j]);
    }
    ret.push_back(akt);
  }
  return ret;
}

myPoint operator*(const std::vector<double> a, const std::vector<myPoint> b){
  myPoint ret;
  for (int i = 0; i < a.size(); i++) {
    ret = ret + (a[i]*b[i]);
  }

  return ret;
}

class CsomoErtekek{
private:
  std::vector<GLdouble> sizes;
  std::vector<GLdouble> ertekek;

  void calcErtekek(){
    int i = 0;
    for (GLdouble akt : sizes) {
      if (i == 0){
        i++;
        continue;
      }
      else{
        ertekek[i] = ertekek[i-1] + akt;
      }
    }
  }

public:
  CsomoErtekek(){
    ertekek.push_back(0.00);
    sizes.push_back(0.00);
  }
  void add(GLdouble a){
    sizes.push_back(a);
    ertekek.push_back(ertekek[ertekek.size()-1]+a);
  }
  GLdouble get(int index){
    return ertekek[index]==0.00?0.00:ertekek[index] / ertekek[ertekek.size()];
  }
  void set(int index, GLdouble value){
    sizes[index] = value;
    calcErtekek();
  }
  int size(){
    return ertekek.size();
  }
};







//std::vector<GLdouble> U;
CsomoErtekek U;
std::vector<myPoint> D;
GLint dragged = -1;
GLint K = 3;

GLdouble calcN(GLdouble j, GLdouble k, GLdouble u){
  if(k == 1.00){
    if ((u>=U.get(j)) && (u<= U.get(j+1))){
      return 1.00;
    }
    else {
      return 0.00;
    }

  }
  else{
    return (((u-U.get(j))/(U.get(j+k-1)))*calcN(j, k-1, u)) + ((U.get(j+k)-u)/(U.get(j+k)-U.get(j+1))*calcN(j+1, k-1, u));
  }

}
GLdouble alfa(GLdouble j, GLdouble l, GLdouble u, GLdouble k){
  return l == 0 ? 1.00 : (u - U.get(j)) / ( U.get(j+k-1)-U.get(j));
}

myPoint calcD(GLdouble j, GLdouble l, GLdouble u, GLdouble k){
  if ( l == 0.00){
    return D[j];
  }
  else {
    return (alfa(j, l, u, k)*calcD(j, l-1, u, k)) + ((1-alfa(j, l, u, k))*calcD(j-1, l-1, u, k));
  }

}

std::vector<myPoint> calcNURBSPoints(){ // main calsulations
  std::vector<myPoint> ret;
  if (D.size()<K+1) return ret;
  
  for (int i = 0; i < D.size()-K; i++) {
      myPoint pontok[K+1];
      for (int j=0; j<K+1;j++){
	pontok[j] = D[i+j];
      }
      
      std::vector<myPoint> ps(&pontok[0],&pontok[0]+K+1);
      for (float j = 0; j < 1; j+=0.1) {
	
        ret.push_back(calcN(i,K,j)*calcD(i,K,j,K));
      }

    }
}
void printUs(){
    for (int i = 0; i< U.size();i++){
      std::cout << U.get(i) << ", ";
      
    }
    std::cout << std::endl;
}


GLdouble rad(GLint i)
{
	return i*(3.1415 / 180.0);
}

GLfloat dist2(myPoint P1, myPoint P2) {
	GLfloat t1 = P1.x - P2.x;
	GLfloat t2 = P1.y - P2.y;
	return t1 * t1 + t2 * t2;
}

void init(void)
{
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0.0, 700, 0.0, 600);
	glShadeModel(GL_FLAT);
	glEnable(GL_LINE_STIPPLE);
	glPointSize(10.0);
  glLineWidth(5.0);
}
void lineSegment(void)
{
	printUs();

	glClear(GL_COLOR_BUFFER_BIT);

	glColor3f(0.5, 0.5, 0.5);

	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINE_STRIP);
	  for (myPoint akt : calcNURBSPoints()) {
	    glVertex2d(akt.x, akt.y);
	  }
	glEnd();

	glutSwapBuffers();
}

GLint getActivePoint1(std::vector<myPoint> p, GLint x, GLint y) {
	GLint s = 100;
	myPoint P(x, y);
	for (GLint i = 0; i < p.size(); i++){
		myPoint value(p[i].x, p[i].y );
		if (dist2(value, P) < s)
			return i;
		}
	return -1;
}

void processMouse(GLint button, GLint action, GLint xMouse, GLint yMouse) {
	GLint i=0;
	if (button == GLUT_LEFT_BUTTON && action == GLUT_DOWN){
		if ((i = getActivePoint1(D, xMouse, 600 - yMouse)) != -1){
			dragged = i;
    }
    else{
      D.push_back(myPoint(xMouse, 600 - yMouse));
      if(D.size() > 1){
        myPoint seged (D[D.size()-1]-D[D.size()-2]);
        U.add(pow(seged.x,2.00)+pow(seged.y, 2.00));
      }
      glutPostRedisplay();
    }
  }
	if (button == GLUT_LEFT_BUTTON && action == GLUT_UP)
		dragged = -1;
}

void processMouseActiveMotion(GLint xMouse, GLint yMouse) {
	if (dragged >= 0) {
    int index = dragged>0 ? dragged : 1;

    D[dragged].x = xMouse;
  	D[dragged].y = 600 - yMouse;

    myPoint seged (D[index]-D[index-1]);
    U.set(index-1,pow(seged.x,2.00)+pow(seged.y, 2.00));
    if (index < D.size()-1){
      myPoint seged (D[index]-D[index+1]);
      U.set(index, pow(seged.x,2.00)+pow(seged.y, 2.00));
    }
    glutPostRedisplay();
	}

}

int main (int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(700, 600);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("B-Spline");
	init();
	glutDisplayFunc(lineSegment);
	glutMouseFunc(processMouse);
	glutMotionFunc(processMouseActiveMotion);
	glutMainLoop();
	return 0;
}
