#include "FrameBuffer.h" 
#include <fstream> 

Frame::Frame(){}
Frame::~Frame(){}

void Frame::insert(float a_box_l, float a_box_r,
	float a_box_u, float a_box_d, float a_box_f,
	float a_box_b, float a_r, float a_g, float a_b,
	float a_a, int *a_idx)
{
	box_l.push_back(a_box_l);
	box_r.push_back(a_box_r);
	box_u.push_back(a_box_u);
	box_d.push_back(a_box_d);
	box_f.push_back(a_box_f);
	box_b.push_back(a_box_b);
	r.push_back(a_r);
	g.push_back(a_g);
	b.push_back(a_b);
	a.push_back(a_a);
	textureIdx.push_back(a_idx);
}

FrameBuffer::FrameBuffer():buf(NULL), frameNum(-1) {}
FrameBuffer::~FrameBuffer()
{
	if (buf)
		delete[] buf;
}

void FrameBuffer::init(int a_n)
{
	buf = new Frame[a_n];
	frameNum = a_n;
}

void FrameBuffer::free()
{
	for (int i = 0; i < frameNum; i++)
		for (vector<int*>::iterator ii = buf[i].textureIdx.begin();
			ii != buf[i].textureIdx.end(); ii++)
			delete (*ii);
	if (buf)
		delete[] buf;
	frameNum = -1;
}

void FrameBuffer::insertFrame(int a_f,
	float a_box_l, float a_box_r, float a_box_u,
	float a_box_d, float a_box_f, float a_box_b,
	float a_r, float a_g, float a_b, float a_a)
{
	int *idx = new int[6];
	for (int i = 0; i < 6; i++)
		idx[i] = rand() % 8;
	buf[a_f].insert(a_box_l, a_box_r, a_box_u,
		a_box_d, a_box_f, a_box_b,
		a_r, a_g, a_b, a_a, idx);
}

void FrameBuffer::saveFrame(string filePath)
{
	int size;
	ofstream out(filePath, ios::out | ios::binary);
	out.write((char*)&frameNum, sizeof(frameNum));
	for (int i = 0; i < frameNum; i++)
	{
		size = (buf[i].r.size());
		out.write((char*)(&size), sizeof(size));
		for (int j = 0; j < buf[i].r.size(); j++)
		{
			out.write((char*)(&buf[i].box_l[j]), sizeof(float));
			out.write((char*)(&buf[i].box_r[j]), sizeof(float));
			out.write((char*)(&buf[i].box_d[j]), sizeof(float));
			out.write((char*)(&buf[i].box_u[j]), sizeof(float));
			out.write((char*)(&buf[i].box_f[j]), sizeof(float));
			out.write((char*)(&buf[i].box_b[j]), sizeof(float));
			out.write((char*)(&buf[i].r[j]), sizeof(float));
			out.write((char*)(&buf[i].g[j]), sizeof(float));
			out.write((char*)(&buf[i].b[j]), sizeof(float));
			out.write((char*)(&buf[i].a[j]), sizeof(float));
		}
	}
	out.close();
}

bool FrameBuffer::loadFrame(string filePath)
{
	int size;
	ifstream file(filePath, ios::in | ios::binary);
	if (!file.is_open())
	{
		return false;
	}

	file.read((char*)(&size), sizeof(size));
	init(size);
	for (int i = 0; i < frameNum; i++)
	{
		file.read((char*)(&size), sizeof(size));
		for (int j = 0; j < size; j++)
		{
			float tmp;
			file.read((char*)(&tmp), sizeof(float));
			buf[i].box_l.push_back(tmp);
			file.read((char*)(&tmp), sizeof(float));
			buf[i].box_r.push_back(tmp);
			file.read((char*)(&tmp), sizeof(float));
			buf[i].box_d.push_back(tmp);
			file.read((char*)(&tmp), sizeof(float));
			buf[i].box_u.push_back(tmp);
			file.read((char*)(&tmp), sizeof(float));
			buf[i].box_f.push_back(tmp);
			file.read((char*)(&tmp), sizeof(float));
			buf[i].box_b.push_back(tmp);
			file.read((char*)(&tmp), sizeof(float));
			buf[i].r.push_back(tmp);
			file.read((char*)(&tmp), sizeof(float));
			buf[i].g.push_back(tmp);
			file.read((char*)(&tmp), sizeof(float));
			buf[i].b.push_back(tmp);
			file.read((char*)(&tmp), sizeof(float));
			buf[i].a.push_back(tmp);
			int *tmpp;
			tmpp = new int[6];
			for (int i = 0; i < 6; i++)
				tmpp[i] = rand() % 8;
			buf[i].textureIdx.push_back(tmpp);
		}
	}
	file.close();
	return true;
}