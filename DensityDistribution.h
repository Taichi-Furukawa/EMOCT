#pragma once

#include <string>
#include <vector>

class DensityDistribution
{
public:
	DensityDistribution();

	/// <summary>Flow-field dataをファイルから読み込みます。</summary>
	/// <param name="name">ファイル名</param>
	DensityDistribution(const std::string& name);
	~DensityDistribution();

	/// <summary>再構成データをファイルに保存します。</summary>
	/// <param name="filePath">ファイル名</param>
	/// <returns>成功した場合はtrue, 失敗した場合はfalse</returns>
	bool save(const std::string& name) const;

	void clear();
	bool empty() const;
	void resize(unsigned int width, unsigned int height, unsigned int depth);

	float at(unsigned int x, unsigned int y, unsigned int z) const;
	float& at(unsigned int x, unsigned int y, unsigned int z);
	float operator()(unsigned int x, unsigned int y, unsigned int z) const;
	float& operator()(unsigned int x, unsigned int y, unsigned int z);

	unsigned int width;
	unsigned int height;
	unsigned int depth;
	std::vector<float> data;
};