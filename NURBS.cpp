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






//std::vector<GLdouble> U;

std::vector<GLdouble> W;
std::vector<myPoint> D;
GLint dragged = -1;
GLint selectedPoint = -1;
GLdouble K = 4.00;

class CsomoErtekek{
private:
  std::vector<GLdouble> sizes;

public:
  CsomoErtekek(){
  }
  
  void set( std::vector<myPoint> a){
    int size=a.size() + K;
    std::vector<GLdouble> temp;
    for ( int i=0; i<K-1; i++ )
        temp.push_back ( 0 );

    for ( int i=0; i<size-2*(K-1); i++ )
        temp.push_back ( i );
    GLdouble last=1;
    if ( temp.at ( temp.size()-1 ) !=0 )
        last=temp.at ( temp.size()-1 );
    for ( int i=0; i<K-1; i++ )
        temp.push_back ( last );
    sizes = temp;
  }
  GLdouble get(int index){
      return sizes[index];
  }
  int size(){
    return sizes.size();
  }
};

CsomoErtekek U;

GLdouble myDiv(GLdouble a, GLdouble b){
  if (a == 0.00 || b == 0.00) return 0.00;
  return a/b;
}




//NURBS funcs----------------------------------------------------------------------------------------------------------------------

GLdouble calcN(GLdouble j, GLdouble k, GLdouble u){
  //std::cout << "calcN j: " << j << " k: " << k << " u: " << u << std::endl;
  if(k == 1.00){
    if ((u>=U.get(j)) && (u< U.get(j+1))){
      return 1.00;
    }
    else {
      return 0.00;
    }

  }
  else{
    return ((myDiv((u-U.get(j)),(U.get(j+k-1))))*calcN(j, k-1, u)) +
           (myDiv((U.get(j+k)-u) , (U.get(j+k)-U.get(j+1))) * calcN(j+1, k-1, u));
  }

}


GLdouble alfa(GLdouble j, GLdouble l, GLdouble u){
  //std::cout << "alfa j: " << j << " l: " << l << " u: " << u << std::endl;
  return l >0.00 ? myDiv((u - U.get(j)) , ( U.get(j+K-l)-U.get(j))) : 1.00;
}

GLdouble weight(GLdouble j, GLdouble l, GLdouble u){
  //std::cout << "weight j: " << j << " l: " << l << " u: " << u << std::endl;
  return l>0.00?alfa(j, l, u)*weight(j, l-1, u) + (1-alfa(j, l, u))*weight(j-1, l-1, u) : W[j];
}

myPoint calcD(GLdouble j, GLdouble l, GLdouble u){
  //std::cout << "calcD j: " << j << " l: " << l << " u: " << u << std::endl;
  if ( l > 0.00){
    return ((alfa(j, l, u)* weight(j, l-1, u)*calcD(j, l-1, u)) + ((1-alfa(j, l, u))* weight(j-1, l-1, u)*calcD(j-1, l-1, u)))* (myDiv(1.00,weight(j, l, u)));
  }
  else {
    return D[j];
    
  }

}

std::vector<myPoint> calcNURBSPoints(){ // main calsulations
  std::vector<myPoint> ret;
  if (D.size()<K+1) return ret;

  for (int i = 0; i < U.size()-1; i++) {
    for (GLdouble u = U.get(i); u < U.get(i+1); u+=0.01){
      ret.push_back(calcD(i,K-1,u));
      
    }
   

  }
  return ret;

}
void printUs(){
    for (int i = 0; i< U.size();i++){
      //std::cout << U.get(i) << ", ";

    }
    //std::cout << std::endl;
}
//END--------------------------------------------------------------------------------------------------------------------------

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
	glClear(GL_COLOR_BUFFER_BIT);

	glColor3f(0.5, 0.5, 0.5);
  glBegin(GL_POINTS);
    for (int i=0; i < D.size();i++) {
      if (i == selectedPoint)
	glColor3f(0, 0.5, 0);
      else
	glColor3f(0.5, 0.5, 0.5);
      glVertex2d(D[i].x, D[i].y);
    }
  glEnd();

	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINE_STRIP);
	  //std::cout << "\n";
	  for (myPoint akt : calcNURBSPoints()) {
	    glVertex2d(akt.x, akt.y);
	  }
	  //std::cout << std::endl;
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
			selectedPoint = i;
			
			glutPostRedisplay();
			
		}
		else{
		  D.push_back(myPoint(xMouse, 600 - yMouse));
		  W.push_back(1.00);
		  if(D.size() > 1){
		    myPoint seged (D[D.size()-2]-D[D.size()-1]);
		    U.set(D);
		  }
		  //printUs();

		  glutPostRedisplay();
		}
	      }
	if (button == GLUT_LEFT_BUTTON && action == GLUT_UP){
		dragged = -1;
	}
	
}

void processMouseActiveMotion(GLint xMouse, GLint yMouse) {
	if (dragged >= 0) {
	  int index = dragged>0 ? dragged : 1;

	  D[dragged].x = xMouse;
	  D[dragged].y = 600 - yMouse;

	  U.set(D);
	  
	  glutPostRedisplay();
	}

}

void keyPressed (unsigned char key, int x, int y) {
    if (selectedPoint>=0){
     switch (key){
	case 'w':
	 W[selectedPoint] += 0.1;
	 glutPostRedisplay();
	break;
	case 's':
	  W[selectedPoint] -= 0.1;
	  glutPostRedisplay();
	break;
	case 'q':
	 K = K <D.size()-1? K+1:K;
	 std::cout << K << std::endl;
	 U.set(D);
	 glutPostRedisplay();
	break;
	case 'a':
	 K = K>2? K-1:K;
	 std::cout << K << std::endl;
	 U.set(D);
	 glutPostRedisplay();
	break;
	 
    }
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
	glutKeyboardFunc(keyPressed);
	glutDisplayFunc(lineSegment);
	glutMouseFunc(processMouse);
	glutMotionFunc(processMouseActiveMotion);
	glutMainLoop();
	return 0;
}
