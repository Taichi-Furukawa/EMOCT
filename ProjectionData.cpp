#include <fstream>
#include <sstream>
#include "ProjectionData.h"

ProjectionData::ProjectionData()
: width(0), height(0), count(0)
{
}

ProjectionData::ProjectionData(const std::string& name)
{
	std::ifstream file(name);
	if (!file)
	{
		// ファイルを開けなかった場合はエラー
		return;
	}

	// ヘッダー読み取り
	{
		std::string string;
		if (!std::getline(file, string))
		{
			// 読み込めなかった場合はエラー
			return;
		}

		std::stringstream stream(string);
		stream >> width >> height >> count;
	}

	angles.resize(count);
	data.resize(count * height * width);
	ids.resize(data.size());

	// データの読み取り
	{
		std::string string;
		for (unsigned int i = 0; i < count; i++)
		{
			if (!std::getline(file, string))
			{
				// 読み込めなかった場合はエラー
				return;
			}
			angles[i].roll = std::stof(string);
			angles[i].pitch = 0.0f;
			angles[i].yaw = 0.0f;

			for (unsigned int y = 0; y < height; y++)
			{
				for (unsigned int x = 0; x < width; x++)
				{
					unsigned int index = i * height * width + y * width + x;

					if (!std::getline(file, string))
					{
						// 読み込めなかった場合はエラー
						return;
					}

					std::stringstream stream(string);
					stream >> data[index] >> ids[index];
				}
			}
		}
	}
}

ProjectionData::~ProjectionData()
{
}

bool ProjectionData::save(const std::string& name) const
{
	std::ofstream file(name, std::ios::out | std::ios::trunc);
	if (!file)
	{
		// ファイルを開けなかった場合はエラー
		return false;
	}

	std::stringstream stream;
	stream << width << " " << height << " " << count << std::endl;
	
	for (unsigned int i = 0; i < count; i++)
	{
		stream << angles[i].roll << std::endl;

		for (unsigned int y = 0; y < height; y++)
		{
			for (unsigned int x = 0; x < width; x++)
			{
				unsigned int index = i * height * width + y * width + x;

				stream << data[index] << " " << ids[index] << std::endl;
			}
		}
	}

	file.write(stream.str().c_str(), stream.str().size());

	return true;
}

bool ProjectionData::empty() const
{
	return data.empty();
}

void ProjectionData::resize(unsigned int width, unsigned int height, unsigned int count)
{
	angles.resize(count);
	data.resize(count * height * width, 0.0f);
	ids.resize(data.size(), 1);

	this->width = width;
	this->height = height;
	this->count = count;
}

float& ProjectionData::at(unsigned int x, unsigned int y, unsigned int i)
{
	return data[x * count * height + i * height + y];
}

float ProjectionData::at(unsigned int x, unsigned int y, unsigned int i) const
{
	return data[x * count * height + i * height + y];
}

int& ProjectionData::id(unsigned int x, unsigned int y, unsigned int i)
{
	return ids[x * count * height + i * height + y];
}

int ProjectionData::id(unsigned int x, unsigned int y, unsigned int i) const
{
	return ids[x * count * height + i * height + y];
}

float& ProjectionData::operator()(unsigned int x, unsigned int y, unsigned int i)
{
	return data[x * count * height + i * height + y];
}

float ProjectionData::operator()(unsigned int x, unsigned int y, unsigned int i) const
{
	return data[x * count * height + i * height + y];
}