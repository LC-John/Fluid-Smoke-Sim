#include <cmath>
#include <utility>
#include <GL/glut.h>
#include <cstdlib>
#include "CFD.h"

using namespace std;
//#define DEBUG

#define I(i, j, k) ((i) + ((NX) + 2) * (j) + ((NX) + 2) * ((NY) + 2) * (k))
#define SWAP(x0, x)      \
	{                    \
		float *tmp = x0; \
		x0 = x;          \
		x = tmp;         \
	}
#define LERP(a, l, h) ((l) + (((h) - (l)) * (a)))
#define FOR_EACH_CELL                 \
	for (i = 1; i <= NX; ++i)         \
	{                                 \
		for (j = 1; j <= NY; ++j)     \
		{                             \
			for (k = 1; k <= NZ; ++k) \
			{
#define END_FOR \
	}           \
	}           \
	}

#define MAX_ITER 5
#define INFLOW_PAUSE 100
#define C_EPS 1e-06f

bool CFD_useObstacles = false;
float CFD_obstX = 0.0f, CFD_obstY = 1.0f, CFD_obstZ = 1.0f, CFD_obstR = 1.0f;
static float dt = 0.001f;
const int MAXSOURCE = 8;
static const float INFLOW_VEL = 1.0f;
static float CFD_vortConf = 3.0f, CFD_buoy = 0.3f;
float nx = 2.f, ny = 2.f, nz = 2.f;
static int CELLS_IN_X = 64, CELLS_IN_Y = 64, CELLS_IN_Z = 64;
float CFD_airT = 1.0f, CFD_smokeT = 300.0f, CFD_coolConst = 1000.0f; //1000
int i, j, k;
float x, y, z;

inline int sgn(float f)
{
	if (f > C_EPS)
		return 1;
	else if (f < -C_EPS)
		return -1;
	else
		return 0;
}

float trilerp(float f000, float f100, float f010, float f110, float f001,
			  float f101, float f011, float f111, float fx, float fy, float fz)
{
	float fx00, fx01, fx10, fx11, fxy0, fxy1, fxyz;

	/* interpolate across x dimension */
	fx00 = LERP(fx, f000, f100);
	fx01 = LERP(fx, f001, f101);
	fx10 = LERP(fx, f010, f110);
	fx11 = LERP(fx, f011, f111);

	/* interpolate across y fimension */
	fxy0 = LERP(fy, fx00, fx10);
	fxy1 = LERP(fy, fx01, fx11);

	/* interpolate across z fimension */
	fxyz = LERP(fz, fxy0, fxy1);

	return fxyz;
}

inline float Noise2D(float *p)
{
	int n = (int)(p[0] + p[1] * 57);
	n = (n << 13) ^ n;
	float r = (1.0f - ((n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0f);
	return r;
}

float PerlinNoise2D(float x, float y, float alpha, float beta, int n)
{
	int i;
	float val, sum = 0;
	float p[2], scale = 1;

	p[0] = x;
	p[1] = y;
	for (i = 0; i < n; i++)
	{
		val = Noise2D(p);
		sum += val / scale;
		scale *= alpha;
		p[0] *= beta;
		p[1] *= beta;
	}
	return (sum);
}
C_SmokeSource::C_SmokeSource() {}
C_SmokeSource::C_SmokeSource(int _id, vec3 &_color, vec3 &_pos, float _r, float _initialD, float _dD, bool _fixedDensity, float _initialVC, float _initialTC) : shape(Ball), id(_id), color(_color), pos(_pos), r(_r), initialD(_initialD), dD(_dD), fixedDensity(_fixedDensity), initialVC(_initialVC), initialTC(_initialTC)
{
}
C_SmokeSource:: ~C_SmokeSource() {}
C_RecSmokeSource::C_RecSmokeSource() : C_SmokeSource()
{
	shape = Cuboid;
}
C_RecSmokeSource::C_RecSmokeSource(int _id, vec3 &_color, vec3 &_pos, float _r, float _initialD, float _dD, bool _fixedDensity, float _initialVC, float _initialTC, float _w, float _h, float _d) : C_SmokeSource(_id, _color, _pos, _r, _initialD, _dD, _fixedDensity, _initialVC, _initialTC), w(_w), h(_w), d(_d) 
{
	shape = Cuboid;
}
C_RecSmokeSource::~C_RecSmokeSource() {}

C_CFD::C_CFD() : NX(CELLS_IN_X), NY(CELLS_IN_Y), NZ(CELLS_IN_Z), buf()
{
	p1x = (NX + 2) / 2;
	p0x = -p1x;
	p1y = (NY + 2) / 2;
	p0y = -p1y;
	p1z = (NZ + 2) / 2;
	p0z = -p1z;
	size = (NX + 2) * (NY + 2) * (NZ + 2);
	fDens = 1.0f;
	fVsc = 0.0f;
	sourceNum = 0;
	fU0 = new float[size];
	fV0 = new float[size];
	fW0 = new float[size];
	fU1 = new float[size];
	fV1 = new float[size];
	fW1 = new float[size];
	fDiv = new float[size];
	fP = new float[size];
	CU = new float[size];
	CV = new float[size];
	CW = new float[size];
	D0 = new float[size];
	D1 = new float[size];
	T0 = new float[size];
	T1 = new float[size];
	BU = new float[size];
	M = new int[size];
	N2D_B = new float[NX * NZ * 3];

	inSourceNum = new int[size];
	DC0 = new float*[MAXSOURCE];
	DC1 = new float*[MAXSOURCE];
	source = new C_SmokeSource*[MAXSOURCE];

	memset(fU0, 0, size * sizeof(float));
	memset(fV0, 0, size * sizeof(float));
	memset(fW0, 0, size * sizeof(float));
	memset(fU1, 0, size * sizeof(float));
	memset(fV1, 0, size * sizeof(float));
	memset(fW1, 0, size * sizeof(float));
	memset(fDiv, 0, size * sizeof(float));
	memset(fP, 0, size * sizeof(float));
	memset(CU, 0, size * sizeof(float));
	memset(CV, 0, size * sizeof(float));
	memset(CW, 0, size * sizeof(float));
	memset(D0, 0, size * sizeof(float));
	memset(D1, 0, size * sizeof(float));
	memset(T0, 0, size * sizeof(float));
	memset(T1, 0, size * sizeof(float));
	memset(BU, 0, size * sizeof(float));
	memset(M, 0, size * sizeof(int));

	memset(inSourceNum, 0, size * sizeof(int));
	for (int i = 0; i < MAXSOURCE; ++i)
	{
		DC0[i] = NULL;
		DC1[i] = NULL;
		source[i] = NULL;
	}

	// if (CFD_useObstacles == true)
	// {
	// 	FOR_EACH_CELL
	// 	x = (float)(i - p1x) - CFD_obstX;
	// 	y = (float)(j - p1y) - CFD_obstY;
	// 	z = (float)(k - p1z) - CFD_obstZ;
	// 	if ((x - 1) * (x - 1) + (y + (CELLS_IN_Y / 4)) * (y + (CELLS_IN_Y / 4)) + z * z < CFD_obstR * CFD_obstR)
	// 	{
	// 		M[I(i, j, k)] = 1;
	// 	}
	// 	END_FOR
	// }
}

C_CFD::~C_CFD()
{
	delete[] fU0;
	delete[] fV0;
	delete[] fW0;
	delete[] fU1;
	delete[] fV1;
	delete[] fW1;
	delete[] fDiv;
	delete[] fP;
	delete[] CU;
	delete[] CV;
	delete[] CW;
	delete[] D0;
	delete[] D1;
	delete[] T0;
	delete[] T1;
	delete[] BU;
	delete[] M;
	delete[] N2D_B;

	delete[] inSourceNum;
	for (int i = 0; i < MAXSOURCE; ++i)
	{
		if (DC0[i] != NULL)
			delete[] DC0[i];
		if (DC1[i] != NULL)
			delete[] DC1[i];
		if (source[i] != NULL)
			delete source[i];
	}
}

int C_CFD::addSource(vec3 &_color, vec3 &_pos, SourceShape shape, float _r, float _initialD, float _dD, bool _fixedDensity, float _initialVC, float _initialTC, float _w , float _h , float _d)
{
	if (sourceNum >= MAXSOURCE)
		return -1;

	int newId = sourceNum;

	C_SmokeSource *newSource = NULL;

	switch (shape)
	{
	case Ball:
		newSource = new C_SmokeSource(newId, _color, _pos, _r, _initialD, _dD, _fixedDensity, _initialVC, _initialTC);
		source[newId] = newSource;
		break;

	case Cuboid:
		newSource = new C_RecSmokeSource(newId, _color, _pos, _r, _initialD, _dD, _fixedDensity, _initialVC, _initialTC, _w, _h, _d);
		source[newId] = (C_SmokeSource *)newSource;
		break;

	default:
		return -1;
		break;
	}

	sourceNum++;
	DC0[newId] = new float[size];
	DC1[newId] = new float[size];
	memset(DC0[newId], 0, size * sizeof(float));
	memset(DC1[newId], 0, size * sizeof(float));

	return newId;
}

void C_CFD::addForce(float *u0, float *v0, float *w0, float *su,
					 float *sv, float *sw)
{
	for (int i = 0; i < size; ++i)
	{
		u0[i] += dt * su[i];
		v0[i] += dt * sv[i];
		w0[i] += dt * sw[i];
	}
}
void C_CFD::diffuse(float vsc, float *a0, float *a1)
{
	float _sum, _invA = 1.0f / (dt * 6.0f + 1);
	setBnd(a1);
	FOR_EACH_CELL
	_sum = a1[I(i - 1, j, k)] + a1[I(i + 1, j, k)] + a1[I(i, j - 1, k)] +
		   a1[I(i, j + 1, k)] + a1[I(i, j, k - 1)] + a1[I(i, j, k + 1)];
	a0[I(i, j, k)] = a1[I(i, j, k)] + dt * vsc * (_sum - a1[I(i, j, k)] * 6);
	END_FOR
	setBnd(a0);
}
void C_CFD::advect(float *a0, float *a1, float *u, float *v, float *w)
{
	int x0, x1, y0, y1, z0, z1;
	FOR_EACH_CELL
	x = i - (dt * u[I(i, j, k)]);
	y = j - (dt * v[I(i, j, k)]);
	z = k - (dt * w[I(i, j, k)]);
	if (x < 0.5f)
	{
		x = 0.5f;
	}
	if (x > NX + 0.5f)
	{
		x = NX + 0.5f;
	}
	x0 = (int)x;
	x1 = x0 + 1; // 4
	if (y < 0.5f)
	{
		y = 0.5f;
	}
	if (y > NY + 0.5f)
	{
		y = NY + 0.5f;
	}
	y0 = (int)y;
	y1 = y0 + 1; // 4
	if (z < 0.5f)
	{
		z = 0.5f;
	}
	if (z > NZ + 0.5f)
	{
		z = NZ + 0.5f;
	}
	z0 = (int)z;
	z1 = z0 + 1; // 4
	a0[I(i, j, k)] = trilerp(a1[I(x0, y0, z0)], a1[I(x1, y0, z0)], a1[I(x0, y1, z0)], a1[I(x1, y1, z0)],
							 a1[I(x0, y0, z1)], a1[I(x1, y0, z1)], a1[I(x0, y1, z1)], a1[I(x1, y1, z1)], x - x0, y - y0, z - z0);
	END_FOR
}

void C_CFD::project(float *u0, float *u1, float *v0, float *v1,
					float *w0, float *w1,
					float *div, float *p, float matDens)
{
	FOR_EACH_CELL
	div[I(i, j, k)] = ((u1[I(i + 1, j, k)] - u1[I(i - 1, j, k)]) + (v1[I(i, j + 1, k)] -
																	v1[I(i, j - 1, k)]) +
					   (w1[I(i, j, k + 1)] - w1[I(i, j, k - 1)])) *
					  0.5f;
	p[I(i, j, k)] = 0.0f;
	END_FOR
	setBnd(div);
	setBnd(p);
	int _cf = 6;
	float lhs;
	for (int l = 0; l < MAX_ITER; ++l)
	{
		FOR_EACH_CELL
		if (CFD_useObstacles)
		{
			int _L = 1 - M[I(i - 1, j, k)];
			int _R = 1 - M[I(i + 1, j, k)];
			int _B = 1 -
					 M[I(i, j - 1, k)];
			int _T = 1 - M[I(i, j + 1, k)];
			int _N = 1 - M[I(i, j, k - 1)];
			int _F = 1 - M[I(i, j, k + 1)];
			_cf = (_L + _R) + (_B + _T) + (_N + _F);
			lhs = _L * p[I(i - 1, j, k)] + _R * p[I(i + 1, j, k)] + _B * p[I(i, j - 1, k)] + _T * p[I(i, j + 1, k)] + _N * p[I(i, j, k - 1)] + _F * p[I(i, j, k + 1)];
		}
		else
		{
			lhs = p[I(i - 1, j, k)] + p[I(i + 1, j, k)] + p[I(i, j - 1, k)] +
				  p[I(i, j + 1, k)] + p[I(i, j, k - 1)] + p[I(i, j, k + 1)];
		}
		if (_cf != 0)
		{
			p[I(i, j, k)] = (lhs - (div[I(i, j, k)] / dt)) / _cf;
		}
		END_FOR
		setBnd(p);
	}
	FOR_EACH_CELL
	u0[I(i, j, k)] = u1[I(i, j, k)] - dt * (((p[I(i + 1, j, k)] - p[I(i - 1, j, k)]) * 0.5f));
	v0[I(i, j, k)] = v1[I(i, j, k)] - dt * (((p[I(i, j + 1, k)] - p[I(i, j - 1, k)]) * 0.5f));
	w0[I(i, j, k)] = w1[I(i, j, k)] - dt * (((p[I(i, j, k + 1)] - p[I(i, j, k - 1)]) * 0.5f));
	END_FOR
	setBnd(u0);
	setBnd(v0);
	setBnd(w0);
}

void C_CFD::setBnd(float *f)
{
	for (i = 1; i <= NX; ++i)
		for (j = 1; j <= NY; ++j)
		{
			f[I(i, j, 0)] = f[I(i, j, 1)];
			f[I(i, j, NZ + 1)] = f[I(i, j, NZ)];
		}
	for (j = 1; j <= NY; ++j)
		for (k = 1; k <= NZ; ++k)
		{
			f[I(0, j, k)] = f[I(1, j, k)];
			f[I(NX + 1, j, k)] = f[I(NX, j, k)];
		}
	for (i = 1; i <= NX; ++i)
		for (k = 1; k <= NZ; ++k)
		{
			f[I(i, 0, k)] = f[I(i, 1, k)];
			f[I(i, NY + 1, k)] = f[I(i, NY, k)];
		}
	for (i = 1; i <= NX; ++i)
	{
		f[I(i, 0, 0)] = (f[I(i, 1, 0)] + f[I(i, 0, 1)]) * 0.5f;
		f[I(i, NY + 1, 0)] = (f[I(i, NY, 0)] + f[I(i, NY + 1, 1)]) * 0.5f;
		f[I(i, 0, NZ + 1)] = (f[I(i, 0, NZ)] + f[I(i, 1, NZ + 1)]) * 0.5f;
		f[I(i, NY + 1, NZ + 1)] = (f[I(i, NY, NZ + 1)] + f[I(i, NY + 1, NZ)]) * 0.5f;
	}
	for (j = 1; j <= NY; ++j)
	{
		f[I(0, j, 0)] = (f[I(1, j, 0)] + f[I(0, j, 1)]) * 0.5f;
		f[I(NX + 1, j, 0)] = (f[I(NX, j, 0)] + f[I(NX + 1, j, 1)]) * 0.5f;
		f[I(0, j, NZ + 1)] = (f[I(0, j, NZ)] + f[I(1, j, NZ + 1)]) * 0.5f;
		f[I(NX + 1, j, NZ + 1)] = (f[I(NX, j, NZ + 1)] + f[I(NX + 1, j, NZ)]) * 0.5f;
	}
	for (k = 1; k <= NZ; ++k)
	{
		f[I(0, 0, k)] = (f[I(0, 1, k)] + f[I(1, 0, k)]) * 0.5f;
		f[I(0, NY + 1, k)] = (f[I(0, NY, k)] + f[I(1, NY + 1, k)]) * 0.5f;
		f[I(NX + 1, 0, k)] = (f[I(NX, 0, k)] + f[I(NX + 1, 1, k)]) * 0.5f;
		f[I(NX + 1, NY + 1, k)] = (f[I(NX + 1, NY, k)] + f[I(NX, NY + 1, k)]) * 0.5f;
	}
	f[I(0, 0, 0)] = (f[I(1, 0, 0)] + f[I(0, 1, 0)] + f[I(0, 0, 1)]) * 0.3333333f;
	f[I(NX + 1, 0, 0)] = (f[I(NX, 0, 0)] + f[I(NX + 1, 1, 0)] + f[I(NX + 1, 0, 1)]) * 0.3333333f;
	f[I(0, NY + 1, 0)] = (f[I(1, NY + 1, 0)] + f[I(0, NY, 0)] + f[I(0, NY + 1, 1)]) * 0.3333333f;
	f[I(0, 0, NZ + 1)] = (f[I(1, 0, NZ + 1)] + f[I(0, 1, NZ + 1)] + f[I(0, 0, NZ)]) * 0.3333333f;
	f[I(NX + 1, NY + 1, 0)] = (f[I(NX, NY + 1, 0)] + f[I(NX + 1, NY, 0)] + f[I(NX + 1, NY + 1, 1)]) * 0.3333333f;
	f[I(0, NY + 1, NZ + 1)] = (f[I(0, NY, NZ + 1)] + f[I(0, NY + 1, NZ)] + f[I(1, NY + 1, NZ + 1)]) * 0.3333333f;
	f[I(NX + 1, 0, NZ + 1)] = (f[I(NX, 0, NZ + 1)] + f[I(NX + 1, 0, NZ)] + f[I(NX + 1, 1, NZ + 1)]) * 0.3333333f;
	f[I(NX + 1, NY + 1, NZ + 1)] = (f[I(NX, NY + 1, NZ + 1)] +
									f[I(NX + 1, NY, NZ + 1)] + f[I(NX + 1, NY + 1, NZ)]) *
								   0.3333333f;
}

void C_CFD::setInflow(float *fu0, float *fv0, float *fw0)
{
	this->genN2DMap(N2D_B, NX, NZ, 0.05f);

	memset(inSourceNum, 0, size * sizeof(int));
	for (int id = 0; id < sourceNum; ++id)
	{
		vec3 sourcePos, recSourcePos;
		float sourceR = 0.f, initialVC = 0.f;
		int sX = 0,eX = 0,sY = 0,eY = 0,sZ = 0,eZ = 0;
		float shW = 0.f, shH = 0.f, shD = 0.f;
		C_RecSmokeSource * recSource = NULL;
		switch (source[id]->shape)
		{
		case Ball:
			sourcePos = source[id]->pos;
			sourceR = source[id]->r;
			initialVC = source[id]->initialVC;
			sX = (int)(sourcePos.x - sourceR), eX = (int)(sourcePos.x + sourceR + 1.f),
				sY = (int)(sourcePos.y - sourceR), eY = (int)(sourcePos.y + sourceR + 1.f),
				sZ = (int)(sourcePos.z - sourceR), eZ = (int)(sourcePos.z + sourceR + 1.f);
			if (sX < 1)
			{
				sX = 1;
			}
			if (eX > NX)
			{
				eX = NX;
			}
			if (sY < 1)
			{
				sY = 1;
			}
			if (eY > NY)
			{
				eY = NY;
			}
			if (sZ < 1)
			{
				sZ = 1;
			}
			if (eZ > NZ)
			{
				eZ = NZ;
			}

			for (j = sY; j <= eY; ++j)
			{
				y = float(j - sourcePos.y);
				for (i = sX; i <= eX; ++i)
				{
					x = float(i - sourcePos.x);
					for (k = sZ; k <= eZ; ++k)
					{
						z = float(k - sourcePos.z);

						if ((x * x + y * y + z * z) < sourceR * sourceR)
						{
							fu0[I(i, j, k)] = (N2D_B[k * NX * 3 + i * 3]) * INFLOW_VEL * 0.5f;
							fv0[I(i, j, k)] = (fv0[I(i, j, k)] * (float)inSourceNum[I(i, j, k)] + ((N2D_B[k * NX * 3 + i * 3 + 1] + 1.0f) * initialVC + INFLOW_VEL)) / (float)(inSourceNum[I(i, j, k)] + 1);
							fw0[I(i, j, k)] = (N2D_B[k * NX * 3 + i * 3 + 2]) * INFLOW_VEL * 0.5f;

							inSourceNum[I(i, j, k)]++;
						}
					}
				}
				break;
		case Cuboid:
			recSource = (C_RecSmokeSource *)source[id];
			recSourcePos = recSource->pos;
			initialVC = recSource->initialVC;
			shW = recSource->w / 2.f, shH = recSource->h / 2.f, shD = recSource->d / 2.f;
			sX = (int)(recSourcePos.x - shW), eX = (int)(recSourcePos.x + shW),
				sY = (int)(recSourcePos.y - shH), eY = (int)(recSourcePos.y + shH),
				sZ = (int)(recSourcePos.z - shD), eZ = (int)(recSourcePos.z + shD);
			if (sX < 1)
			{
				sX = 1;
			}
			if (eX > NX)
			{
				eX = NX;
			}
			if (sY < 1)
			{
				sY = 1;
			}
			if (eY > NY)
			{
				eY = NY;
			}
			if (sZ < 1)
			{
				sZ = 1;
			}
			if (eZ > NZ)
			{
				eZ = NZ;
			}

			for (j = sY; j <= eY; ++j)
				for (i = sX; i <= eX; ++i)
					for (k = sZ; k <= eZ; ++k)
					{
						fu0[I(i, j, k)] = (N2D_B[k * NX * 3 + i * 3]) * INFLOW_VEL * 0.5f;
						fv0[I(i, j, k)] = (fv0[I(i, j, k)] * (float)inSourceNum[I(i, j, k)] + ((N2D_B[k * NX * 3 + i * 3 + 1] + 1.0f) * initialVC + INFLOW_VEL)) / (float)(inSourceNum[I(i, j, k)] + 1);
						fw0[I(i, j, k)] = (N2D_B[k * NX * 3 + i * 3 + 2]) * INFLOW_VEL * 0.5f;

						inSourceNum[I(i, j, k)]++;
					}
			break;
		default:
			break;
			}
		}

		// for (i = 1; i <= NX; ++i)
		// 	for (k = 1; k <= NZ; ++k)
		// 	{
		// 		x = (float)(i - p1x);
		// 		z = (float)(k - p1z);
		// 		if ((x * x + z * z) < eR * eR)
		// 		{
		// 			fu0[I(i, 1, k)] = (N2D_B[k * NX * 3 + i * 3]) * INFLOW_VEL * 0.5f;
		// 			fv0[I(i, 1, k)] = ((N2D_B[k * NX * 3 + i * 3 + 1] + 1.0f) * 5.0f + INFLOW_VEL);
		// 			fw0[I(i, 1, k)] = (N2D_B[k * NX * 3 + i * 3 + 2]) * INFLOW_VEL * 0.5f;
		// 		}
		// 	}

		// memset(M, 0, size * sizeof(int));
		// if (CFD_useObstacles == true)
		// {
		// 	FOR_EACH_CELL
		// 	x = (float)(i - p1x) - CFD_obstX;
		// 	y = (float)(j - p1y) - CFD_obstY;
		// 	z = (float)(k - p1z) - CFD_obstZ;
		// 	if ((x - 1) * (x - 1) + (y + (CELLS_IN_Y / 4)) * (y + (CELLS_IN_Y / 4)) + z * z < CFD_obstR * CFD_obstR)
		// 	{
		// 		M[I(i, j, k)] = 1;
		// 	}
		// 	END_FOR
		// }
	}
}


vec3 C_CFD::getVorticity(int i, int j, int k)
{
	float DwDy, DvDz, DuDz, DwDx, DvDx, DuDy;
	DwDy = DvDz = DuDz = DwDx = DvDx = DuDy = 0.0f;
	DwDy = (fW0[I(i, j + 1, k)] - fW0[I(i, j - 1, k)]) * 0.5f;
	DvDz = (fV0[I(i, j, k + 1)] - fV0[I(i, j, k - 1)]) * 0.5f;
	DuDz = (fU0[I(i, j, k + 1)] - fU0[I(i, j, k - 1)]) * 0.5f;
	DwDx = (fW0[I(i + 1, j, k)] - fW0[I(i - 1, j, k)]) * 0.5f;
	DvDx = (fV0[I(i + 1, j, k)] - fV0[I(i - 1, j, k)]) * 0.5f;
	DuDy = (fU0[I(i, j + 1, k)] - fU0[I(i, j - 1, k)]) * 0.5f;
	return vec3(DwDy - DvDz, DuDz - DwDx, DvDx - DuDy);
}

inline float nrm1(const vec3 &v)
{
	return (fabs(v.x) + fabs(v.y) + fabs(v.z));
}

inline float nrm2(const vec3 &v)
{
	return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

void C_CFD::updateCurls(float *cu, float *cv, float *cw)
{
	memset(cu, 0, size * sizeof(float));
	memset(cv, 0, size * sizeof(float));
	memset(cw, 0, size * sizeof(float));
	float DoDx = 0.0f, DoDy = 0.0f, DoDz = 0.0f, len, invL;
	vec3 _v;
	FOR_EACH_CELL
	if (i >= 2 && i <= NX - 2)
		DoDx = (nrm2(getVorticity(i + 1, j, k)) - nrm2(getVorticity(i - 1, j, k))) * 0.5f;
	if (j >= 2 && j <= NY - 2)
		DoDy = (nrm2(getVorticity(i, j + 1, k)) - nrm2(getVorticity(i, j - 1, k))) * 0.5f;
	if (k >= 2 && k <= NZ - 2)
		DoDz = (nrm2(getVorticity(i, j, k + 1)) - nrm2(getVorticity(i, j, k - 1))) * 0.5f;
	len = sqrt(DoDx * DoDx + DoDy * DoDy + DoDz * DoDz);

	if (len >= C_EPS)
	{
		invL = 1.0f / len;
		DoDx *= invL;
		DoDy *= invL;
		DoDz *= invL;
	}
	else
	{
		DoDx = DoDy = DoDz = 0.0f;
	}

	_v = getVorticity(i, j, k);
	cu[I(i, j, k)] = (DoDy * _v.z - _v.y * DoDz) * CFD_vortConf;
	cv[I(i, j, k)] = (DoDz * _v.x - _v.z * DoDx) * CFD_vortConf;
	cw[I(i, j, k)] = (DoDx * _v.y - _v.x * DoDy) * CFD_vortConf;
	END_FOR
}

void C_CFD::updateDensity(int frame)
{

	for (int id = 0; id < sourceNum; ++id)
	{

		vec3 sourcePos, recSourcePos;
		float sourceR = 0.f;
		int sX = 0, eX = 0, sY = 0, eY = 0, sZ = 0, eZ = 0;
		float shW = 0.f, shH = 0.f, shD = 0.f;
		C_RecSmokeSource * recSource = NULL;
		switch (source[id]->shape)
		{
		case Ball:
			sourcePos = source[id]->pos;
			sourceR = source[id]->r;
			sX = (int)(sourcePos.x - sourceR), eX = (int)(sourcePos.x + sourceR + 1.f),
				sY = (int)(sourcePos.y - sourceR), eY = (int)(sourcePos.y + sourceR + 1.f),
				sZ = (int)(sourcePos.z - sourceR), eZ = (int)(sourcePos.z + sourceR + 1.f);
			if (sX < 1)
			{
				sX = 1;
			}
			if (eX > NX)
			{
				eX = NX;
			}
			if (sY < 1)
			{
				sY = 1;
			}
			if (eY > NY)
			{
				eY = NY;
			}
			if (sZ < 1)
			{
				sZ = 1;
			}
			if (eZ > NZ)
			{
				eZ = NZ;
			}

			for (j = sY; j <= eY; ++j)
			{
				y = float(j - sourcePos.y);
				for (i = sX; i <= eX; ++i)
				{
					x = float(i - sourcePos.x);
					for (k = sZ; k <= eZ; ++k)
					{
						z = float(k - sourcePos.z);

						if ((x * x + y * y + z * z) < sourceR * sourceR)
						{
							if (source[id]->fixedDensity)
							{
								DC0[id][I(i, j, k)] = source[id]->initialD;
							}
							else
							{
								if (frame == 0)
								{
									DC0[id][I(i, j, k)] = source[id]->initialD;
								}
								else
								{
									DC0[id][I(i, j, k)] += source[id]->dD;
								}
							}
							if (DC0[id][I(i, j, k)] < C_EPS)
							{
								DC0[id][I(i, j, k)] = 0.f;
							}
						}
					}
				}
				break;
		case Cuboid:
			recSource = (C_RecSmokeSource *)source[id];
			recSourcePos = recSource->pos;
			shW = recSource->w / 2.f, shH = recSource->h / 2.f, shD = recSource->d / 2.f;
			sX = (int)(recSourcePos.x - shW), eX = (int)(recSourcePos.x + shW),
				sY = (int)(recSourcePos.y - shH), eY = (int)(recSourcePos.y + shH),
				sZ = (int)(recSourcePos.z - shD), eZ = (int)(recSourcePos.z + shD);
			if (sX < 1)
			{
				sX = 1;
			}
			if (eX > NX)
			{
				eX = NX;
			}
			if (sY < 1)
			{
				sY = 1;
			}
			if (eY > NY)
			{
				eY = NY;
			}
			if (sZ < 1)
			{
				sZ = 1;
			}
			if (eZ > NZ)
			{
				eZ = NZ;
			}

			for (j = sY; j <= eY; ++j)
				for (i = sX; i <= eX; ++i)
					for (k = sZ; k <= eZ; ++k)
					{
						if (recSource->fixedDensity)
						{
							DC0[id][I(i, j, k)] = recSource->initialD;
						}
						else
						{
							if (frame == 0)
							{
								DC0[id][I(i, j, k)] = recSource->initialD;
							}
							else
							{
								DC0[id][I(i, j, k)] += recSource->dD;
							}
						}
						if (DC0[id][I(i, j, k)] < C_EPS)
						{
							DC0[id][I(i, j, k)] = 0.f;
						}
					}
			break;
		default:
			break;
			}
			setBnd(DC0[id]);
			SWAP(DC0[id], DC1[id]);
			advect(DC0[id], DC1[id], fU0, fV0, fW0);
		}

		// for (i = 1; i <= NX; ++i)
		// 	for (k = 1; k <= NZ; ++k)
		// 	{
		// 		x = (float)(i - p1x);
		// 		z = (float)(k - p1z);
		// 		if ((x * x + z * z) < eR * eR)
		// 		{
		// 			D0[I(i, 1, k)] = 1.0f;
		// 		}
		// 	}
		// if (CFD_useObstacles)
		// {
		// FOR_EACH_CELL
		// 	x = (i - p1x) - CFD_obstX;
		// 	y = (j - p1y) - CFD_obstY;
		// 	z = (k - p1z) - CFD_obstZ;
		// 	if ((x - 1) * (x - 1) + (y + (CELLS_IN_Y / 4)) * (y + (CELLS_IN_Y / 4)) + z * z < CFD_obstR * CFD_obstR)
		// 	{
		// 		D0[I(i, j, k)] = 1.0f;
		// 	}
		// END_FOR
		// }
		// setBnd(D0);
		// SWAP(D0, D1);
		// advect(D0, D1, fU0, fV0, fW0);
	}
	memset(D0, 0, size * sizeof(float));
	memset(D1, 0, size * sizeof(float));

	for (int id = 0; id < sourceNum; ++id)
	{
		FOR_EACH_CELL
			int ite = I(i, j, k);
			D0[I(i, j, k)] += DC0[id][I(i, j, k)];
			//printf_s("%f", D0[ite]);
		END_FOR
	}
	for (int id = 0; id < sourceNum; ++id)
	{
		FOR_EACH_CELL
			D1[I(i, j, k)] += DC1[id][I(i, j, k)];
		END_FOR
	}
}


void C_CFD::updateVelocity()
{
	this->addExtForces();
	//if (CFD_useObstacles)
	//{
	//	FOR_EACH_CELL
	//	x = (i - p1x) - CFD_obstX;
	//	y = (j - p1y) - CFD_obstY;
	//	z = (k - p1z) - CFD_obstZ;
	//	if ((x - 1) * (x - 1) + (y + (CELLS_IN_Y / 4)) * (y + (CELLS_IN_Y / 4)) + z * z < CFD_obstR * CFD_obstR)
	//	{
	//		float _x = x - 1;
	//		float _y = y + (CELLS_IN_Y / 4);
	//		float _z = z;
	//		float _l = sqrt(_x * _x + _y * _y + _z * _z);
	//		if (fabs(_l) <= C_EPS)
	//		{
	//			_l = C_EPS * sgn(_l);
	//		}
	//		float _invL = 5.0f / _l;
	//		fU0[I(i, j, k)] = _x * _invL;
	//		fV0[I(i, j, k)] = _y * _invL;
	//		fW0[I(i, j, k)] = _z * _invL;
	//	}
	//	END_FOR
	//}
	if (fVsc >= C_EPS)
	{
		SWAP(fU0, fU1);
		SWAP(fV0, fV1);
		SWAP(fW0, fW1);
		diffuse(fVsc, fU0, fU1);
		diffuse(fVsc, fV0, fV1);
		diffuse(fVsc, fW0, fW1);
		SWAP(fU0, fU1);
		SWAP(fV0, fV1);
		SWAP(fW0, fW1);
		project(fU0, fU1, fV0, fV1, fW0, fW1, fDiv, fP, fDens);
		//setObstacleBnd();
	}
	SWAP(fU0, fU1);
	SWAP(fV0, fV1);
	SWAP(fW0, fW1);
	advect(fU0, fU1, fU1, fV1, fW1);
	advect(fV0, fV1, fU1, fV1, fW1);
	advect(fW0, fW1, fU1, fV1, fW1);
	//setObstacleBnd();
	SWAP(fU0, fU1);
	SWAP(fV0, fV1);
	SWAP(fW0, fW1);
	project(fU0, fU1, fV0, fV1, fW0, fW1, fDiv, fP, fDens);
	//setObstacleBnd();
}

void C_CFD::updateTemperature()
{
	static float offset = 0.0f;
	offset += 0.0025f;

	memset(inSourceNum, 0, size * sizeof(int));
	for (int id = 0; id < sourceNum; ++id)
	{
		vec3 sourcePos, recSourcePos;
		float sourceR = 0.f, initialTC = 0.f;
		int sX = 0, eX = 0, sY = 0, eY = 0, sZ = 0, eZ = 0;
		float shW = 0.f, shH = 0.f, shD = 0.f;
		C_RecSmokeSource * recSource = NULL;
		switch (source[id]->shape)
		{
		case Ball:
			sourcePos = source[id]->pos;
			sourceR = source[id]->r;
			initialTC = source[id]->initialTC;
			sX = (int)(sourcePos.x - sourceR), eX = (int)(sourcePos.x + sourceR + 1.f),
				sY = (int)(sourcePos.y - sourceR), eY = (int)(sourcePos.y + sourceR + 1.f),
				sZ = (int)(sourcePos.z - sourceR), eZ = (int)(sourcePos.z + sourceR + 1.f);
			if (sX < 1)
			{
				sX = 1;
			}
			if (eX > NX)
			{
				eX = NX;
			}
			if (sY < 1)
			{
				sY = 1;
			}
			if (eY > NY)
			{
				eY = NY;
			}
			if (sZ < 1)
			{
				sZ = 1;
			}
			if (eZ > NZ)
			{
				eZ = NZ;
			}

			for (j = sY; j <= eY; ++j)
			{
				y = float(j - sourcePos.y);
				for (i = sX; i <= eX; ++i)
				{
					x = float(i - sourcePos.x);
					for (k = sZ; k <= eZ; ++k)
					{
						z = float(k - sourcePos.z);

						if ((x * x + y * y + z * z) < sourceR * sourceR)
						{
							if ((x * x + y * y + z * z) < sourceR * sourceR * 0.6f)
							{
								int _x = (i + 1) * (1.0f / NX) * 2.0f;
								int _y = (k + 1) * (1.0f / NZ) * 2.0f;
								float val = PerlinNoise2D(_x + offset, _y + offset, 1.5f, 1.0f, 4);
								T0[I(i, j, k)] = (T0[I(i, j, k)] * (float)inSourceNum[I(i, j, k)] + (fabs(val) + 0.95f) * initialTC * 0.45f) / (float)(inSourceNum[I(i, j, k)] + 1);
							}
							else
							{
								T0[I(i, j, k)] = (T0[I(i, j, k)] * (float)inSourceNum[I(i, j, k)] + initialTC * 0.45f) / (float)(inSourceNum[I(i, j, k)] + 1);
							}
							inSourceNum[I(i, j, k)]++;
						}
					}
				}
				break;
		case Cuboid:
			recSource = (C_RecSmokeSource *)source[id];
			recSourcePos = recSource->pos;
			initialTC = recSource->initialTC;
			shW = recSource->w / 2.f, shH = recSource->h / 2.f, shD = recSource->d / 2.f;
			sX = (int)(recSourcePos.x - shW), eX = (int)(recSourcePos.x + shW),
				sY = (int)(recSourcePos.y - shH), eY = (int)(recSourcePos.y + shH),
				sZ = (int)(recSourcePos.z - shD), eZ = (int)(recSourcePos.z + shD);
			if (sX < 1)
			{
				sX = 1;
			}
			if (eX > NX)
			{
				eX = NX;
			}
			if (sY < 1)
			{
				sY = 1;
			}
			if (eY > NY)
			{
				eY = NY;
			}
			if (sZ < 1)
			{
				sZ = 1;
			}
			if (eZ > NZ)
			{
				eZ = NZ;
			}

			for (j = sY; j <= eY; ++j)
				for (i = sX; i <= eX; ++i)
					for (k = sZ; k <= eZ; ++k)
					{
						if (i >= (int)(recSourcePos.x - shW*0.8f) && i < (int)(recSourcePos.x + shW*0.8f) && j >= (int)(recSourcePos.y - shH*0.8f) && j < (int)(recSourcePos.y + shH*0.8f) && k >= (int)(recSourcePos.z - shD*0.8f) && k < (int)(recSourcePos.z + shD*0.8f))
						{
							int _x = (i + 1) * (1.0f / NX) * 2.0f;
							int _y = (k + 1) * (1.0f / NZ) * 2.0f;
							float val = PerlinNoise2D(_x + offset, _y + offset, 1.5f, 1.0f, 4);
							T0[I(i, j, k)] = (T0[I(i, j, k)] * (float)inSourceNum[I(i, j, k)] + (fabs(val) + 0.95f) * initialTC * 0.45f) / (float)(inSourceNum[I(i, j, k)] + 1);
						}
						else
						{
							T0[I(i, j, k)] = (T0[I(i, j, k)] * (float)inSourceNum[I(i, j, k)] + initialTC * 0.45f) / (float)(inSourceNum[I(i, j, k)] + 1);
						}

						inSourceNum[I(i, j, k)]++;
					}
			break;
		default:
			break;
			}
		}

		// for (i = 1; i <= NX; ++i)
		// 	for (k = 1; k <= NZ; ++k)
		// 	{
		// 		x = (float)(i - p1x);
		// 		z = (float)(k - p1z);
		// 		float xSqr = x * x;
		// 		float zSqr = z * z;
		// 		float radiusSqr = eR * eR;
		// 		if (xSqr + zSqr < radiusSqr)
		// 		{
		// 			if (xSqr + zSqr < radiusSqr * 0.6f)
		// 			{
		// 				int _x = (i + 1) * (1.0f / NX) * 2.0f;
		// 				int _y = (k + 1) * (1.0f / NZ) * 2.0f;
		// 				float val = PerlinNoise2D(_x + offset, _y + offset, 1.5f, 1.0f, 4);
		// 				T0[I(i, 1, k)] = (fabs(val) + 0.95f) * CFD_smokeT * 0.45f;
		// 			}
		// 			else
		// 			{
		// 				T0[I(i, 1, k)] = CFD_smokeT * 0.45f;
		// 			}
		// 		}
		// 	}
		// FOR_EACH_CELL
		// if (CFD_useObstacles)
		// {
		// 	x = (i - p1x) - CFD_obstX;
		// 	y = (j - p1y) - CFD_obstY;
		// 	z = (k - p1z) - CFD_obstZ;
		// 	if ((x - 1) * (x - 1) + (y + (CELLS_IN_Y / 4)) * (y + (CELLS_IN_Y / 4)) + z * z < CFD_obstR * CFD_obstR)
		// 	{
		// 		T0[I(i, j, k)] = CFD_smokeT;
		// 	}
		// }
		// END_FOR
		setBnd(T0);
		SWAP(T0, T1);
		advect(T0, T1, fU0, fV0, fW0);
		FOR_EACH_CELL
			if (T0[I(i, j, k)] > CFD_airT)
				T0[I(i, j, k)] -= CFD_coolConst * 0.2f * pow((T0[I(i, j, k)] - CFD_airT) / (CFD_smokeT - CFD_airT), 5);
		if (T0[I(i, j, k)] < (CFD_smokeT - CFD_airT) * 0.1f)
			T0[I(i, j, k)] *= 0.75f;
		END_FOR setBnd(T0);
	}
}


void C_CFD::addExtForces()
{
	this->updateCurls(CU, CV, CW);
	addForce(fU0, fV0, fW0, CU, CV, CW);
	FOR_EACH_CELL
	float _t = T0[I(i, j, k)];
	BU[I(i, j, k)] = CFD_buoy * (T0[I(i, j, k)]);
	fV0[I(i, j, k)] += dt * BU[I(i, j, k)];
	END_FOR
}

void C_CFD::update(int frame)
{
	this->setInflow(fU0, fV0, fW0);
	float _maxFVel = 0.0f;
	FOR_EACH_CELL
	float _FVel = fabs(fU0[I(i, j, k)]) + fabs(fV0[I(i, j, k)]) +
				  fabs(fW0[I(i, j, k)]);
	if (_maxFVel < _FVel)
		_maxFVel = _FVel;
	END_FOR
	dt = 1.0f / (_maxFVel + (0.5f + 0.5f + 0.5f));
	this->updateVelocity();
	this->updateDensity(frame);
	this->updateTemperature();
}

void C_CFD::setObstacleBnd()
{
	if (CFD_useObstacles)
	{
		FOR_EACH_CELL
		if (M[I(i, j, k)] == 1)
		{
			int _emptyNeighbors = 0;
			if (M[I(i - 1, j, k)] == 0)
			{
				++_emptyNeighbors;
				fU0[I(i, j, k)] = -fU0[I(i - 1, j, k)];
				fV0[I(i, j, k)] = -fV0[I(i -
											 1,
										 j, k)];
				fW0[I(i, j, k)] = -fW0[I(i - 1, j, k)];
			}
			else if (M[I(i + 1, j, k)] == 0)
			{
				++_emptyNeighbors;
				fU0[I(i, j, k)] = -fU0[I(i + 1, j, k)];
				fV0[I(i, j, k)] = -fV0[I(i + 1, j, k)];
				fW0[I(i, j, k)] = -fW0[I(i + 1, j, k)];
			}
			if (M[I(i, j - 1, k)] == 0)
			{
				++_emptyNeighbors;
				fU0[I(i, j, k)] = -fU0[I(i, j - 1, k)];
				fV0[I(i, j, k)] = -fV0[I(i, j - 1, k)];
				fW0[I(i, j, k)] = -fW0[I(i, j - 1, k)];
			}
			else if (M[I(i, j + 1, k)] == 0)
			{
				++_emptyNeighbors;
				fU0[I(i, j, k)] = -fU0[I(i, j + 1, k)];
				fV0[I(i, j, k)] = -fV0[I(i, j + 1, k)];
				fW0[I(i, j, k)] = -fW0[I(i, j + 1, k)];
			}
			if (M[I(i, j, k - 1)] == 0)
			{
				++_emptyNeighbors;
				fU0[I(i, j, k)] = -fU0[I(i, j, k - 1)];
				fV0[I(i, j, k)] = -fV0[I(i, j, k - 1)];
				fW0[I(i, j, k)] = -fW0[I(i, j, k - 1)];
			}
			else if (M[I(i, j, k + 1)] == 0)
			{
				++_emptyNeighbors;
				fU0[I(i, j, k)] = -fU0[I(i, j, k + 1)];
				fV0[I(i, j, k)] = -fV0[I(i, j, k + 1)];
				fW0[I(i, j, k)] = -fW0[I(i, j, k + 1)];
			}
			if (_emptyNeighbors == 0 || _emptyNeighbors > 1)
			{
				fU0[I(i, j, k)] = 0.0f;
				fV0[I(i, j, k)] = 0.0f;
				fW0[I(i, j, k)] = 0.0f;
			}
		}
		END_FOR
	}
}

void C_CFD::genN2DMap(float *n2D, int nW, int nH, float offset)
{
	static float _o = 0.0f;
	_o += offset;
	for (j = 0; j < nH; ++j)
		for (i = 0; i < nW; ++i)
		{
			x = i * (1.0f / nW) * 2.0f;
			y = j * (1.0f / nH) * 2.0f;
			n2D[j * 3 * nW + i * 3] = PerlinNoise2D(x + _o, y, 1.5f, 1.0f, 4);
			n2D[j * 3 * nW + i * 3 + 1] = PerlinNoise2D(x, y + _o, 1.5f, 1.0f, 0.4);
			n2D[j * 3 * nW + i * 3 + 2] = PerlinNoise2D(x - _o, y, 1.5f, 1.0f, 40);
		}
}

void C_CFD::draw_box(float box_l, float box_r, float box_d,
	float box_u, float box_f, float box_b,
	float surf_r, float surf_g, float surf_b, float surf_a,
	GLuint* NoForAll, int *t_idx)
{
	glColor4f(surf_r, surf_g, surf_b, surf_a);
	float surfaceAmbient[4] = { 0., 0., 0., 0. };
	float surfaceDiffuse[4] = { surf_r, surf_g, surf_b, surf_a };
	float surfaceSpecular[4] = { 0., 0., 0., 0. };
	float surfaceEmission[4] = { 0., 0., 0., 0. };
	float surfaceShiness = 0.;
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glBindTexture(GL_TEXTURE_2D, NoForAll[t_idx[0]]);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, surfaceAmbient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, surfaceDiffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, surfaceSpecular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, surfaceEmission);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, surfaceShiness);
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
	glBindTexture(GL_TEXTURE_2D, NoForAll[t_idx[1]]);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, surfaceAmbient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, surfaceDiffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, surfaceSpecular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, surfaceEmission);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, surfaceShiness);
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
	glBindTexture(GL_TEXTURE_2D, NoForAll[t_idx[2]]);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, surfaceAmbient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, surfaceDiffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, surfaceSpecular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, surfaceEmission);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, surfaceShiness);
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
	glBindTexture(GL_TEXTURE_2D, NoForAll[t_idx[3]]);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, surfaceAmbient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, surfaceDiffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, surfaceSpecular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, surfaceEmission);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, surfaceShiness);
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
	glBindTexture(GL_TEXTURE_2D, NoForAll[t_idx[4]]);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, surfaceAmbient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, surfaceDiffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, surfaceSpecular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, surfaceEmission);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, surfaceShiness);
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
	glBindTexture(GL_TEXTURE_2D, NoForAll[t_idx[5]]);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, surfaceAmbient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, surfaceDiffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, surfaceSpecular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, surfaceEmission);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, surfaceShiness);
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
}

void C_CFD::writeFrameBuffer(int frameN, vec3 canvasL, vec3 canvasH,
	GLuint* NoForAll, int textureCnt)
{
	
}

void C_CFD::render(int op, int frameN, vec3 canvasL, vec3 canvasH,
	GLuint* NoForAll, int textureCnt)
{
	switch (op)
	{
	case 0:
	case 1:
		for (int a_s = 0; a_s < sourceNum; a_s++)
		{
			float *tmpD = DC0[a_s];
			vec3 color = source[a_s]->color;
			FOR_EACH_CELL
				float tmpx = (canvasH.x - canvasL.x) * (float)i / (float)CELLS_IN_X + canvasL.x;
			float tmpy = (canvasH.y - canvasL.y) * (float)j / (float)CELLS_IN_Y + canvasL.y;
			float tmpz = (canvasH.z - canvasL.z) * (float)k / (float)CELLS_IN_Z + canvasL.z;
			float dx = (canvasH.x - canvasL.x) / (float)CELLS_IN_X / nx;
			float dy = (canvasH.y - canvasL.y) / (float)CELLS_IN_Y / ny;
			float dz = (canvasH.z - canvasL.z) / (float)CELLS_IN_Z / nz;
			float tmpgs000 = ((tmpD[I(i, j, k)] / 2 > 1) ? 1 : tmpD[I(i, j, k)] / 2);
			float tmpgs001 = ((tmpD[I(i, j, k + 1)] / 2 > 1) ? 1 : tmpD[I(i, j, k + 1)] / 2);
			float tmpgs010 = ((tmpD[I(i, j + 1, k)] / 2 > 1) ? 1 : tmpD[I(i, j + 1, k)] / 2);
			float tmpgs100 = ((tmpD[I(i + 1, j, k)] / 2 > 1) ? 1 : tmpD[I(i + 1, j, k)] / 2);
			float tmpgs011 = ((tmpD[I(i, j + 1, k + 1)] / 2 > 1) ? 1 : tmpD[I(i, j + 1, k + 1)] / 2);
			float tmpgs101 = ((tmpD[I(i + 1, j, k + 1)] / 2 > 1) ? 1 : tmpD[I(i + 1, j, k + 1)] / 2);
			float tmpgs110 = ((tmpD[I(i + 1, j + 1, k)] / 2 > 1) ? 1 : tmpD[I(i + 1, j + 1, k)] / 2);
			float tmpgs111 = ((tmpD[I(i + 1, j + 1, k + 1)] / 2 > 1) ? 1 : tmpD[I(i + 1, j + 1, k + 1)] / 2);
			const float gs_epslion = 1e-3;
			int cnt = (tmpgs000 < gs_epslion) + (tmpgs001 < gs_epslion) + (tmpgs010 < gs_epslion)
				+ (tmpgs100 < gs_epslion) + (tmpgs110 < gs_epslion) + (tmpgs101 < gs_epslion)
				+ (tmpgs011 < gs_epslion) + (tmpgs111 < gs_epslion);
			if (i >= CELLS_IN_X || j >= CELLS_IN_Y || k >= CELLS_IN_Z || cnt > 6)
				continue;
			for (int ii = 0; ii < nx; ii++)
				for (int jj = 0; jj < ny; jj++)
					for (int kk = 0; kk < nz; kk++)
					{
						float tmpgs = trilerp(tmpgs000, tmpgs100, tmpgs010, tmpgs110,
							tmpgs001, tmpgs101, tmpgs011, tmpgs111,
							(float)ii / nx, (float)jj / ny, (float)kk / nz);
						//if (tmpgs < gs_epslion)
						//	continue;
						if (op == 0)
						{
							int *idx = new int[6];
							for (int i = 0; i < 6; i++)
								idx[i] = textureCnt;
							draw_box(tmpx + ii*dx, tmpx + ii*dx + dx,
								tmpy + jj*dy, tmpy + jj*dy + dy,
								tmpz + kk*dz, tmpz + kk*dz + dz,
								color.x, color.y, color.z, tmpgs,
								NoForAll, idx);
							delete[] idx;
						}
						else
							buf.insertFrame(frameN, tmpx + ii*dx, tmpx + ii*dx + dx,
								tmpy + jj*dy, tmpy + jj*dy + dy,
								tmpz + kk*dz, tmpz + kk*dz + dz,
								color.x, color.y, color.z, tmpgs);
					}
			END_FOR
		}
		break;
	case 2:
		for (int i = 0; i < buf.buf[frameN].r.size(); i++)
		{
			draw_box(buf.buf[frameN].box_l[i], buf.buf[frameN].box_r[i],
				buf.buf[frameN].box_d[i], buf.buf[frameN].box_u[i],
				buf.buf[frameN].box_f[i], buf.buf[frameN].box_b[i],
				buf.buf[frameN].r[i], buf.buf[frameN].g[i],
				buf.buf[frameN].b[i], buf.buf[frameN].a[i],
				NoForAll, buf.buf[frameN].textureIdx[i]);
		}
		break;
	}
}

void C_CFD::initBuf(int a_frameN)
{
	buf.init(a_frameN);
}
void C_CFD::freeBuf()
{
	buf.free();
}

void C_CFD::saveBuf(string filePath)
{
	buf.saveFrame(filePath);
}
void C_CFD::loadBuf(string filePath)
{
	buf.loadFrame(filePath);
}