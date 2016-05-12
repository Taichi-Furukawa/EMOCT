#pragma once

#include <string>
#include <vector>

/// <summary>
/// PRJ(.prj)形式の投影データを扱います。
/// </summary>
class ProjectionData
{
public:
	struct rotation_angles
	{
		float roll;		// rotation angles of z-axis
		float pitch;	// rotation angles of x-axis
		float yaw;		// rotation angles of y-axis
	};

public:
	ProjectionData();

	/// <summary>投影データをファイルから読み込みます。</summary>
	/// <param name="name">ファイル名</param>
	ProjectionData(const std::string& name);
	~ProjectionData();

	/// <summary>投影データをファイルに保存します。</summary>
	/// <param name="name">ファイル名</param>
	/// <returns>成功した場合はtrue, 失敗した場合はfalse</returns>
	bool save(const std::string& name) const;

	bool empty() const;
	void resize(unsigned int width, unsigned int height, unsigned int count);

	float& at(unsigned int x, unsigned int y, unsigned int i);
	float at(unsigned int x, unsigned int y, unsigned int i) const;
	int& id(unsigned int x, unsigned int y, unsigned int i);
	int id(unsigned int x, unsigned int y, unsigned int i) const;

	float& operator()(unsigned int x, unsigned int y, unsigned int i);
	float operator()(unsigned int x, unsigned int y, unsigned int i) const;

	unsigned int width;
	unsigned int height;
	unsigned int count;
	std::vector<rotation_angles> angles;
	std::vector<float> data;
	std::vector<int> ids;
};