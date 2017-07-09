#pragma once
#ifndef __TEXTURE_H__
#define __TEXTURE_H__

#include <gl/glut.h>

class P_TEXTURE
{
	const float PI;
	
	float PerlinNoise2D(float x, float y);
	float CosInterpolate(float a, float b, float t);
	
public:
	const static int TEXTURE_W = 512;
	const static int TEXTURE_H = 512;
	const static int TEXTURE_C = 4;
	const int Max_Texture_Number;

	GLubyte image[TEXTURE_H][TEXTURE_W][TEXTURE_C];

	P_TEXTURE();
	~P_TEXTURE();
	void init();
	void load(GLuint NoForAll[], int id);
};

#endif