#pragma once

#ifndef __FRAMEBUFFER_H__
#define __FRAMEBUFFER_H__

#include <vector>
#include <string>
using namespace std;

class Frame
{
public:
	vector<float> box_l, box_r, box_u, box_d, box_f, box_b;
	vector<float> r, g, b, a;
	vector<int*> textureIdx;

	Frame();
	~Frame();
	void insert(float a_box_l, float a_box_r, float a_box_u,
		float a_box_d, float a_box_f, float a_box_b,
		float a_r, float a_g, float a_b, float a_a, int *a_idx);
};

class FrameBuffer
{
public:
	Frame *buf;
	int frameNum;

	FrameBuffer();
	~FrameBuffer();
	void init(int a_n);
	void free();
	void insertFrame(int a_f,
		float a_box_l, float a_box_r, float a_box_u,
		float a_box_d, float a_box_f, float a_box_b,
		float a_r, float a_g, float a_b, float a_a);
	void saveFrame(string filePath);
	bool loadFrame(string filePath);
};

#endif

