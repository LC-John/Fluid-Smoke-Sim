#include <GL/glut.h>
#include <vector>
#include <ctime>
#include <Windows.h>
#include <cstdio>
#include "CFD.h"
#include "texture.h"

using namespace std;

#define	Window_H	(512)
#define Window_W	(512)
#define Window_D	(512)
#define FPS			(5)
#define Screen_W	(glutGet(GLUT_SCREEN_WIDTH))
#define Screen_H	(glutGet(GLUT_SCREEN_HEIGHT))
#define PI			(3.1415926)

C_CFD grid = C_CFD();
P_TEXTURE texture[8];
GLuint NoForAll[8];
bool Pause = false;
float rhoAngle = 2.0f;
float phiAngle = 0;
int currFrame = 0;
const int totalFrame = 100;

#define BMP_Header_Length 54  
void grab(const char *fileName)
{
	FILE*    pDummyFile;  //指向另一bmp文件，用于复制它的文件头和信息头数据  
	FILE*    pWritingFile;  //指向要保存截图的bmp文件  
	GLubyte* pPixelData;    //指向新的空的内存，用于保存截图bmp文件数据  
	GLubyte  BMP_Header[BMP_Header_Length];
	GLint    i, j;
	GLint    PixelDataLength;   //BMP文件数据总长度  

								// 计算像素数据的实际长度  
	i = Window_W * 3;   // 得到每一行的像素数据长度  
	while (i % 4 != 0)      // 补充数据，直到i是的倍数  
		++i;
	PixelDataLength = i * Window_H;  //补齐后的总位数  

									 // 分配内存和打开文件  
	pPixelData = (GLubyte*)malloc(PixelDataLength);
	if (pPixelData == 0)
		exit(0);

	pDummyFile = fopen("dummy.bmp", "rb");//只读形式打开  
	if (pDummyFile == 0)
		exit(0);

	pWritingFile = fopen(fileName, "wb"); //只写形式打开  
	if (pWritingFile == 0)
		exit(0);

	//把读入的bmp文件的文件头和信息头数据复制，并修改宽高数据  
	fread(BMP_Header, sizeof(BMP_Header), 1, pDummyFile);  //读取文件头和信息头，占据54字节  
	fwrite(BMP_Header, sizeof(BMP_Header), 1, pWritingFile);
	fseek(pWritingFile, 0x0012, SEEK_SET); //移动到0X0012处，指向图像宽度所在内存  
	i = Window_W;
	j = Window_H;
	fwrite(&i, sizeof(i), 1, pWritingFile);
	fwrite(&j, sizeof(j), 1, pWritingFile);

	// 读取当前画板上图像的像素数据  
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);  //设置4位对齐方式  
	glReadBuffer(GL_FRONT);
	glReadPixels(0, 0, Window_W, Window_H,
		GL_BGR_EXT, GL_UNSIGNED_BYTE, pPixelData);

	// 写入像素数据  
	fseek(pWritingFile, 0, SEEK_END);
	//把完整的BMP文件数据写入pWritingFile  
	fwrite(pPixelData, PixelDataLength, 1, pWritingFile);

	// 释放内存和关闭文件  
	fclose(pDummyFile);
	fclose(pWritingFile);
	free(pPixelData);
}

void init(int argc, char ** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(Window_W, Window_H);
	glutInitWindowPosition((Screen_W - Window_W) / 2, (Screen_H - Window_H) / 2);
	glutCreateWindow("Smoke-Sim");
	glClearColor(0, 0, 0, 0);
	//glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();
	glOrtho(-Window_W / 2, Window_W / 2,
		-Window_H / 2, Window_H / 2,
		-Window_D / 2, Window_D / 2);
	//glMatrixMode(GL_MODELVIEW);
	//glLoadIdentity();
	glEnable(GL_BLEND);
	glReadBuffer(GL_BACK);
	glDrawBuffer(GL_BACK);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);
	
	float DiffuseLight1[4] = { 1, 1, 1, 1 };
	float DiffuseLight2[4] = { 1, 1, 1, 1 };
	float DiffuseLight3[4] = { 1, 1, 1, 1 };
	float DiffuseLight4[4] = { 1, 1, 1, 1 };
	float AmbientLight[4] = { 0, 0, 0, 0 };

	float Position0[4] = { 0, 0, 0, 1 };
	glLightfv(GL_LIGHT0, GL_AMBIENT, AmbientLight);
	glLightfv(GL_LIGHT0, GL_POSITION, Position0);
	glEnable(GL_LIGHT0);
	float Position1[4] = { -Window_W / 2, Window_H / 2, Window_D / 2, 1 };
	glLightfv(GL_LIGHT1, GL_DIFFUSE, DiffuseLight1);
	glLightfv(GL_LIGHT1, GL_POSITION, Position1);
	glEnable(GL_LIGHT1);
	float Position2[4] = { -Window_W / 2, Window_H / 2, -Window_D / 2, 1 };
	glLightfv(GL_LIGHT2, GL_DIFFUSE, DiffuseLight2);
	glLightfv(GL_LIGHT2, GL_POSITION, Position2);
	glEnable(GL_LIGHT2);
	float Position3[4] = { Window_W / 2, Window_H / 2, -Window_D / 2, 1 };
	glLightfv(GL_LIGHT3, GL_DIFFUSE, DiffuseLight3);
	glLightfv(GL_LIGHT3, GL_POSITION, Position3);
	glEnable(GL_LIGHT3);
	float Position4[4] = { Window_W / 2, Window_H / 2, Window_D / 2, 1 };
	glLightfv(GL_LIGHT4, GL_DIFFUSE, DiffuseLight4);
	glLightfv(GL_LIGHT4, GL_POSITION, Position4);
	glEnable(GL_LIGHT4);
	glEnable(GL_LIGHTING);

	glGenTextures(8, NoForAll);
	for (int i = 0; i < 8; i++)
	{
		texture[i].init();
		texture[i].load(NoForAll, i);
	}
	glEnable(GL_TEXTURE_2D);
	/*
	grid.addSource(vec3(1, 0, 0),			vec3(26, 2, 26), Ball, 2.f, 1.f, 0.01f, true, 5.f);
	grid.addSource(vec3(1, 122./255., 0),	vec3(28, 2, 28), Ball, 2.f, 1.f, 0.01f, true, 5.f);
	grid.addSource(vec3(1, 1, 0),			vec3(30, 2, 30), Ball, 2.f, 1.f, 0.01f, true, 5.f);
	grid.addSource(vec3(0, 1, 0),			vec3(32, 2, 32), Ball, 2.f, 1.f, 0.01f, true, 5.f);
	grid.addSource(vec3(0, 255./255., 1),	vec3(34, 2, 34), Ball, 2.f, 1.f, 0.01f, true, 5.f);
	grid.addSource(vec3(0, 0, 1),			vec3(36, 2, 36), Ball, 2.f, 1.f, 0.01f, true, 5.f);
	grid.addSource(vec3(255./255., 0, 1),	vec3(38, 2, 38), Ball, 2.f, 1.f, 0.01f, true, 5.f);*/
	/*
	grid.addSource(vec3(1, 1, 1), vec3(28, 1, 28), Ball, 3.f, 0.4f, 0.01f, true, 4.f);
	grid.addSource(vec3(1, 1, 1), vec3(36, 1, 36), Ball, 3.f, 0.4f, 0.01f, true, 4.f);
	grid.initBuf(totalFrame);
	for (int i = 0; i < totalFrame; i++)
	{
		time_t t0, t1, t2;
		t0 = clock();
		grid.render(1, i,
			vec3(-Window_W / 3, -Window_H / 3, -Window_D / 3),
			vec3(Window_W / 3, Window_H / 3, Window_D / 3),
			NoForAll, 8);
		t1 = clock();
		grid.update(i);
		t2 = clock();
		printf_s("frame %d, %dms, %d ms\n", i, (int)(t1 - t0), (int)(t2 - t1));
	}
	grid.saveBuf(".savefile");
	*/
	grid.loadBuf(".savefile");
}

void Draw_Box(float box_l, float box_r, float box_d,
	float box_u, float box_f, float box_b,
	float surf_r, float surf_g, float surf_b, float surf_a,
	float line_r, float line_g, float line_b, float line_a,
	GLuint* NoForAll, int *t_idx)
{
	glDisable(GL_TEXTURE_2D);

	float surfaceAmbient[4] = { 0., 0., 0., 0. };
	float surfaceDiffuse[4] = { 0.1, 0.1, 0.1, 0.1 };
	float surfaceSpecular[4] = { 0.1, 0.1, 0.1, 0.1 };
	float surfaceEmission[4] = { 0., 0., 0., 0. };
	float surfaceShiness = 1.;
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, surfaceAmbient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, surfaceDiffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, surfaceSpecular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, surfaceEmission);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, surfaceShiness);

	glColor4f(line_r, line_g, line_b, line_a);
	glBegin(GL_LINES);
	glVertex3f(box_l, box_u, box_b); glVertex3f(box_l, box_u, box_f);
	glVertex3f(box_l, box_d, box_b); glVertex3f(box_l, box_d, box_f);
	glVertex3f(box_l, box_u, box_b); glVertex3f(box_l, box_d, box_b);
	glVertex3f(box_l, box_u, box_f); glVertex3f(box_l, box_d, box_f);
	glVertex3f(box_r, box_u, box_b); glVertex3f(box_r, box_u, box_f);
	glVertex3f(box_r, box_d, box_b); glVertex3f(box_r, box_d, box_f);
	glVertex3f(box_r, box_u, box_b); glVertex3f(box_r, box_d, box_b);
	glVertex3f(box_r, box_u, box_f); glVertex3f(box_r, box_d, box_f);
	glVertex3f(box_l, box_u, box_b); glVertex3f(box_r, box_u, box_b);
	glVertex3f(box_l, box_d, box_b); glVertex3f(box_r, box_d, box_b);
	glVertex3f(box_l, box_u, box_f); glVertex3f(box_r, box_u, box_f);
	glVertex3f(box_l, box_d, box_f); glVertex3f(box_r, box_d, box_f);
	glEnd();

	glColor4f(0, 0, 0, 1);
	glBindTexture(GL_TEXTURE_2D, NoForAll[t_idx[0]]);
	glBegin(GL_POLYGON);				// left
	glNormal3f(-1, 0, 0);
	glTexCoord2f(0, 0);
	glVertex3f(box_l, box_d, box_f);
	glTexCoord2f(0, 1);
	glVertex3f(box_l, box_d, box_b);
	glTexCoord2f(1, 1);
	glVertex3f(box_l, box_u, box_b);
	glTexCoord2f(1, 0);
	glVertex3f(box_l, box_u, box_f);
	glEnd();
	glColor4f(surf_r, surf_g, surf_b, surf_a);
	glBindTexture(GL_TEXTURE_2D, NoForAll[t_idx[1]]);
	glBegin(GL_POLYGON);				// right
	glNormal3f(1, 0, 0);
	glTexCoord2f(0, 0);
	glVertex3f(box_r, box_d, box_f);
	glTexCoord2f(0, 1);
	glVertex3f(box_r, box_d, box_b);
	glTexCoord2f(1, 1);
	glVertex3f(box_r, box_u, box_b);
	glTexCoord2f(1, 0);
	glVertex3f(box_r, box_u, box_f);
	glEnd();
	glColor4f(surf_r, surf_g, surf_b, surf_a);
	glBindTexture(GL_TEXTURE_2D, NoForAll[t_idx[2]]);
	glBegin(GL_POLYGON);				// up
	glNormal3f(0, 1, 0);
	glTexCoord2f(0, 0);
	glVertex3f(box_l, box_u, box_f);
	glTexCoord2f(0, 1);
	glVertex3f(box_r, box_u, box_f);
	glTexCoord2f(1, 1);
	glVertex3f(box_r, box_u, box_b);
	glTexCoord2f(1, 0);
	glVertex3f(box_l, box_u, box_b);
	glEnd();
	glColor4f(0, 0, 0, 1);
	glBindTexture(GL_TEXTURE_2D, NoForAll[t_idx[3]]);
	glBegin(GL_POLYGON);				// down
	glNormal3f(0, -1, 0);
	glTexCoord2f(0, 0);
	glVertex3f(box_l, box_d, box_f);
	glTexCoord2f(0, 1);
	glVertex3f(box_r, box_d, box_f);
	glTexCoord2f(1, 1);
	glVertex3f(box_r, box_d, box_b);
	glTexCoord2f(1, 0);
	glVertex3f(box_l, box_d, box_b);
	glEnd();
	glColor4f(surf_r, surf_g, surf_b, surf_a);
	glBindTexture(GL_TEXTURE_2D, NoForAll[t_idx[4]]);
	glBegin(GL_POLYGON);				// front
	glNormal3f(0, 0, -1);
	glTexCoord2f(0, 0);
	glVertex3f(box_l, box_u, box_f);
	glTexCoord2f(0, 1);
	glVertex3f(box_r, box_u, box_f);
	glTexCoord2f(1, 1);
	glVertex3f(box_r, box_d, box_f);
	glTexCoord2f(1, 0);
	glVertex3f(box_l, box_d, box_f);
	glEnd();
	glColor4f(0, 0, 0, 1);
	glBindTexture(GL_TEXTURE_2D, NoForAll[t_idx[5]]);
	glBegin(GL_POLYGON);				// back
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0);
	glVertex3f(box_l, box_u, box_b);
	glTexCoord2f(0, 1);
	glVertex3f(box_r, box_u, box_b);
	glTexCoord2f(1, 1);
	glVertex3f(box_r, box_d, box_b);
	glTexCoord2f(1, 0);
	glVertex3f(box_l, box_d, box_b);
	glEnd();

	glEnable(GL_TEXTURE_2D);
}

bool if_grab = true;
void myDisplay()
{
	glClear(GL_COLOR_BUFFER_BIT);

	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_BLEND);
	glBindTexture(GL_TEXTURE_2D, NoForAll[0]);
	glBegin(GL_POLYGON);
	glTexCoord2d(0, 0); glVertex3f(-Window_W / 2, -Window_H / 2, 0);
	glTexCoord2d(0, 1); glVertex3f(-Window_W / 2, Window_H / 2, 0);
	glTexCoord2d(1, 1); glVertex3f(Window_W / 2, Window_H / 2, 0);
	glTexCoord2d(1, 0); glVertex3f(Window_W / 2, -Window_H / 2, 0);
	glEnd();

	/*
	glPushMatrix();

	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_BLEND);
	glRotated(phiAngle, 1, 0, 0);
	glRotatef(rhoAngle, 0, 1, 0);
	grid.render(2, currFrame,
		vec3(-Window_W / 4, -Window_H / 4, -Window_D / 4),
		vec3(Window_W / 4, Window_H / 4, Window_D / 4),
		NoForAll, 8);
	int idx[] = { 0, 0, 0, 0, 0, 0 };
	Draw_Box(-Window_W / 4, Window_W / 4,
		-Window_H / 3, Window_H / 3,
		-Window_D / 4, Window_D / 4,
		0.3, 0.3, 0.3, 0.3, 
		0., 0., 0., 1, NoForAll, idx);

	glPopMatrix();
	*/
	/*
	int pixelSize = Window_W * Window_H * 3;
	GLubyte *pixelBuf = new GLubyte[pixelSize+10];
	glReadPixels(0, 0, Window_W, Window_H, GL_RGB, GL_UNSIGNED_BYTE, pixelBuf);
	for (int i = 0; i < pixelSize; i+=3)
	{
		if (i/Window_W % 50 == 0 && i%Window_W % 50 == 0)
			printf_s("(%d, %d), (%d, %d, %d)\n",
				i / Window_W, i%Window_W, pixelBuf[i], pixelBuf[i + 1], pixelBuf[i + 2]);
	}
	for (int i = 0; i < pixelSize; i++)
		pixelBuf[i] = 255 - pixelBuf[i];
	for (int i = 0; i < pixelSize; i += 3)
	{
		if (i / Window_W % 50 == 0 && i%Window_W % 50 == 0)
			printf_s("\t(%d, %d), (%d, %d, %d)\n",
				i / Window_W, i%Window_W, pixelBuf[i], pixelBuf[i + 1], pixelBuf[i + 2]);
	}
	glDrawPixels(Window_W, Window_H, GL_RGB, GL_UNSIGNED_BYTE, pixelBuf);
	*/

	glutSwapBuffers();
	
	if (if_grab)
	{
		char *sNum = new char[10];
		_itoa(currFrame, sNum, 10);
		string s = string(sNum);
		s += ".bmp";
		grab(s.data());
		delete[] sNum;
	}
}

int add_or_sub = 1;
void myTimer(int val)
{
	clock_t t_begin, t_paint, t_update;
	t_begin = clock();
	/*
	if (currFrame == 0)
	{
		grid.addSource(vec3(1, 0, 0), vec3(14, 1, 16), Ball, 3.5f, 0.5f, 0.01f, false);
		grid.addSource(vec3(0, 1, 0), vec3(16, 1, 16), Ball, 3.5f, 0.5f, 0.01f, false);
		grid.addSource(vec3(0, 0, 1), vec3(18, 1, 16), Ball, 3.5f, 0.5f, 0.01f, false);
	}
	*/
	
	if (!Pause)
	{
		//grid.update(currFrame);
		currFrame++;
		if (currFrame >= totalFrame)
		{
			if (if_grab)
				if_grab = false;
			currFrame -= totalFrame;
		}
		rhoAngle += 0.5;
		phiAngle += add_or_sub * 0.2;
		if (rhoAngle >= 360)
			rhoAngle -= 360;
		if (phiAngle >= 20)
		{
			phiAngle = 20;
			add_or_sub = -add_or_sub;
		}
		if (phiAngle <= 0)
		{
			phiAngle = 0;
			add_or_sub = -add_or_sub;
		}
	}

	t_update = clock();
	myDisplay();
	t_paint = clock();
	printf_s("%d, %dms, %dms\n", currFrame, (int)(t_update - t_begin), (int)(t_paint-t_update));
	glutTimerFunc(1000.f / FPS, myTimer, 1);
}

void myKeyBoard(unsigned char key, int x, int y)
{
	if (key == 'p' || key == 'P')
	{
		Pause = !Pause;
	}
	if (key == 'd' || key == 'D')
	{
		rhoAngle += 5;
		if (rhoAngle >= 360)
			rhoAngle -= 360;
	}
	if (key == 'a' || key == 'A')
	{
		rhoAngle -= 5;
		if (rhoAngle < 0)
			rhoAngle += 360;
	}
	if (key == 'w' || key == 'W')
	{
		phiAngle += 5;
		if (phiAngle >= 90)
			phiAngle = 90;

	}
	if (key == 's' || key == 'S')
	{
		phiAngle -= 5;
		if (phiAngle <= 0)
			phiAngle = 0;
	}
	if (key == 'c' || key == 'C')
	{
		grab("grab.bmp");
	}
}

int main(int argc, char** argv)
{
	init(argc, argv);
	glutDisplayFunc(myDisplay);
	glutTimerFunc(1000./FPS,myTimer,1);
	glutKeyboardFunc(myKeyBoard);
	glutMainLoop();

	system("pause");
	return 0;
}