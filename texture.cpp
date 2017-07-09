#include "texture.h"
#include <cmath>
#include <cstdlib>

P_TEXTURE::P_TEXTURE() : PI(3.1415926f), Max_Texture_Number(8)
{}

P_TEXTURE::~P_TEXTURE() {}

float P_TEXTURE::PerlinNoise2D(float x, float y)
{
	int n = (int)(x + y * 57);
	n = (n << 13) ^ n;
	float r = (1.0f - ((n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0f);
	return r;
}

float P_TEXTURE::CosInterpolate(float a, float b, float t)
{
	float ft = t * PI;
	float f = (1 - cos(ft)) * .5;
	return a*(1 - f) + b*f;
}

void P_TEXTURE::init()
{
	const int freq = 16;
	float tmpImage[(TEXTURE_H) / freq + 1][(TEXTURE_W) / freq + 1];
	int rx = rand() % 100, ry = rand() % 100;
	for (int i = 0; i <= TEXTURE_H; i += freq)
		for (int j = 0; j <= TEXTURE_W; j += freq)
			tmpImage[i / freq][j / freq] = PerlinNoise2D(i / freq + rx, j / freq + ry);
	for (int i = 0; i < TEXTURE_H; i++)
		for (int j = 0; j < TEXTURE_W; j++)
		{
			if ((i % freq == 0) && (j % freq == 0))
			{
				float c = (tmpImage[i / freq][j / freq] * 255.f);
				image[i][j][0] = (GLubyte)(c);
				image[i][j][1] = (GLubyte)(c);
				image[i][j][2] = (GLubyte)(c);
				image[i][j][3] = (GLubyte)(c);
			}
			else
			{
				float a = CosInterpolate(tmpImage[i / freq][j / freq], tmpImage[i / freq + 1][j / freq], (float)(i%freq) / (float)freq);
				float b = CosInterpolate(tmpImage[i / freq][j / freq], tmpImage[i / freq + 1][j / freq], (float)(i%freq) / (float)freq);
				float c = CosInterpolate(a, b, (float)(j%freq) / (float)freq) * 255.f;
				image[i][j][0] = (GLubyte)(c);
				image[i][j][1] = (GLubyte)(c);
				image[i][j][2] = (GLubyte)(c);
				image[i][j][3] = (GLubyte)(c);
			}
		}
}

void P_TEXTURE::load(GLuint NoForAll[], int id)
{
	glBindTexture(GL_TEXTURE_2D, NoForAll[id]);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
		TEXTURE_W, TEXTURE_H, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
}