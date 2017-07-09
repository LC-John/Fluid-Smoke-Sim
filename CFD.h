#ifndef __CFD_H__
#define __CFD_H__

#include "FrameBuffer.h"

#define MAX_ITER 5
#define INFLOW_PAUSE 100
#define C_EPS 1e-06f

inline float Noise2D(float *p);
float PerlinNoise2D(float x, float y, float alpha, float beta, int n);

enum SourceShape
{
	Ball,
	Cuboid
};

class vec3
{
public:
	int x, y, z;
	vec3(int _x = 0, int _y = 0, int _z = 0) : x(_x), y(_y), z(_z) {}
};

class C_SmokeSource
{
public:
	int id;
	vec3 color;
	vec3 pos;
	float r;
	float initialD, dD;
	bool fixedDensity;
	float initialVC, initialTC;
	SourceShape shape;
	C_SmokeSource();
	C_SmokeSource::C_SmokeSource(int _id, vec3 &_color, vec3 &_pos, float _r, float _initialD, float _dD, bool _fixedDensity, float _initialVC, float _initialTC);
	virtual ~C_SmokeSource();
};

class C_RecSmokeSource : public C_SmokeSource
{
public:
	float w, h, d;
	C_RecSmokeSource();
	C_RecSmokeSource(int _id, vec3 &_color, vec3 &_pos, float _r, float _initialD, float _dD, bool _fixedDensity, float _initialVC, float _initialTC, float _w, float _h, float _d);
	~C_RecSmokeSource();
};


class C_CFD
{
	int NX, NY, NZ;
	int p0x, p0y, p0z, p1x, p1y, p1z;
	int size;
	float fDens, fVsc;
	float *fU0, *fV0, *fW0, *fU1, *fV1, *fW1, *fDiv, *fP;
	float *CU, *CV, *CW, *D0, *D1, *T0, *T1, *BU, *N2D_B;
	int *M;
	int sourceNum;
	int *inSourceNum;
	float **DC0, **DC1;
	C_SmokeSource **source;
	FrameBuffer buf;

public:
	C_CFD();
	~C_CFD();
	int addSource(vec3 &_color, vec3 &_pos, SourceShape shape = Ball, float _r = 2.5f, float _initialD = 1.0f, float _dD = 0.05f, bool _fixedDensity = false, float _initialVC = 5.0f, float _initialTC = 300.0f, float _w = 2.0f, float _h = 2.0f, float _d = 2.0f);

	void addForce(float *u0, float *v0, float *w0, float *su,
		float *sv, float *sw);
	void diffuse(float vsc, float *a0, float *a1);
	void advect(float *a0, float *a1, float *u, float *v,
float *w);
	void project(float *u0, float *u1, float *v0, float *v1,
		float *w0, float *w1,
		float *div, float *p, float matDens);
	void setBnd(float *f);
	void setInflow(float *fu0, float *fv0, float *fw0);
	vec3 getVorticity(int i, int j, int k);
	void updateCurls(float *cu, float *cv, float *cw);
	void updateDensity(int frame);
	void updateVelocity();
	void updateTemperature();
	void addExtForces();
	void update(int frame);
	void render(int op, int frameN, vec3 canvasL, vec3 canvasH,
		GLuint* NoForAll, int textureCnt);
		// op = 0, real time; op = 1, compute each frame;
		// op = 2, draw each frame
	void writeFrameBuffer(int frameN, vec3 canvasL, vec3 canvasH,
		GLuint* NoForAll, int textureCnt);
	void setObstacleBnd();
	void genN2DMap(float *n2D, int nW, int nH, float offset);
	void draw_box(float box_l, float box_r, float box_d,
		float box_u, float box_f, float box_b,
		float surf_r, float surf_g, float surf_b, float surf_a,
		GLuint* NoForAll, int* t_idx);
	void initBuf(int a_frameN);
	void freeBuf();
	void saveBuf(string filePath);
	void loadBuf(string filePath);
};

#endif __CFD_H__