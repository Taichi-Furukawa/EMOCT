#include <fstream>
#include <iostream>
#include "DensityDistribution.h"

DensityDistribution::DensityDistribution()
: width(0), height(0), depth(0)
{
}

DensityDistribution::DensityDistribution(const std::string& name)
: width(0), height(0), depth(0)
{
	std::ifstream file(name, std::ios::in | std::ios::binary);
	if (!file)
	{
		// ファイルを開けなかった場合はエラー
		return;
	}

	int n;
	file.read(reinterpret_cast<char*>(&width), sizeof(unsigned int));
	file.read(reinterpret_cast<char*>(&height), sizeof(unsigned int));
	file.read(reinterpret_cast<char*>(&depth), sizeof(unsigned int));
	file.read(reinterpret_cast<char*>(&n), sizeof(int));

	data.resize(width * height * depth, 0.0f);

	for (auto& density : data)
	{
		file.read(reinterpret_cast<char*>(&density), sizeof(float));
		density -= 1.0;
	}
}

DensityDistribution::~DensityDistribution()
{
}

bool DensityDistribution::save(const std::string& name) const
{
	int n = 1;
	std::ofstream file(name, std::ios::out | std::ios::binary | std::ios::trunc);
	if (!file)
	{
		// ファイルを開けなかった場合はエラー
		return false;
	}

	file.write(reinterpret_cast<const char*>(&width), sizeof(unsigned int));
	file.write(reinterpret_cast<const char*>(&height), sizeof(unsigned int));
	file.write(reinterpret_cast<const char*>(&depth), sizeof(unsigned int));
	file.write(reinterpret_cast<const char*>(&n), sizeof(int));
	
	for (auto density : data)
	{
		float number = density + 1.0f;
		file.write(reinterpret_cast<const char*>(&number), sizeof(float));
	}

	return true;
}

void DensityDistribution::clear()
{
	width = 0;
	height = 0;
	depth = 0;
	data.clear();
}

bool DensityDistribution::empty() const
{
	return data.empty();
}

void DensityDistribution::resize(unsigned int width, unsigned int height, unsigned int depth)
{
	this->width = width;
	this->height = height;
	this->depth = depth;
	data.resize(width * height * depth, 0.0f);
}

float DensityDistribution::at(unsigned int x, unsigned int y, unsigned int z) const
{
	return data[z * width * height + y * width + x];
}

float& DensityDistribution::at(unsigned int x, unsigned int y, unsigned int z)
{
	return data[z * width * height + y * width + x];
}

float& DensityDistribution::operator()(unsigned int x, unsigned int y, unsigned int z)
{
	return data[z * width * height + y * width + x];
}

float DensityDistribution::operator()(unsigned int x, unsigned int y, unsigned int z) const
{
	return data[z * width * height + y * width + x];
}