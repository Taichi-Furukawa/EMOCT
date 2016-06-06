//=================================================================================//
//                                                                                 //
//  ilab library                                                                   //
//                                                                                 //
//  Copyright (C) 2011-2016 Terry                                                  //
//                                                                                 //
//  This file is a portion of the ilab library. It is distributed under the MIT    //
//  License, available in the root of this distribution and at the following URL.  //
//  http://opensource.org/licenses/mit-license.php                                 //
//                                                                                 //
//=================================================================================//

#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <fstream>

namespace ilab
{

/// <summary>
/// 投影データファイル(*.prj)を扱います。
/// </summary>
class projection final
{
public:
	projection()
		: m_width(0), m_height(0), m_counts(0)
	{
	}

	projection(size_t _width, size_t _height, size_t _counts)
		: m_width(_width), m_height(_height), m_counts(_counts), m_quantities(_counts * _height * _width, 0.0f), m_identities(_counts * _height * _width, blank_type::quantity), m_angles(_counts, 0.0f)
	{
	}

	/// <summary>投影データをファイルから読み込みます。</summary>
	/// <param name="_name">ファイル名</param>
	projection(const std::string& _name)
		: m_width(0), m_height(0), m_counts(0)
	{
		std::ifstream file(_name);
		if (!file)
		{
			/// ファイルを開けなかった場合はエラー
			return;
		}

		/// ヘッダー読み取り
		std::string string;
		if (!std::getline(file, string))
		{
			return;
		}

		std::stringstream header(string);
		uint32_t width, height, counts;
		header >> width >> height >> counts;

		m_width = static_cast<size_t>(width);
		m_height = static_cast<size_t>(height);
		m_counts = static_cast<size_t>(counts);
		m_angles.resize(m_counts);
		m_quantities.resize(m_counts * m_height * m_width);
		m_identities.resize(m_quantities.size());

		/// データ読み取り
		for (size_t i = 0; i < m_counts; i++)
		{
			if (!std::getline(file, string))
			{
				return;
			}
			m_angles[i] = std::stof(string);

			for (size_t y = 0; y < m_height; y++)
			{
				for (size_t x = 0; x < m_width; x++)
				{
					size_t index = i * m_height * m_width + y * m_width + x;
					if (!std::getline(file, string))
					{
						return;
					}

					std::stringstream block(string);
					float quantity;
					int32_t identity;
					block >> quantity >> identity;

					m_quantities[index] = quantity;
					m_identities[index] = static_cast<blank_type>(identity);
				}
			}
		}
	}

	~projection() = default;

	/// <summary>投影データをファイルに保存します。</summary>
	/// <param name="_name">ファイル名</param>
	/// <returns>成功した場合はtrue, 失敗した場合はfalse</returns>
	bool save(const std::string& _name) const
	{
		std::ofstream file(_name, std::ios::out | std::ios::trunc);
		if (!file)
		{
			return false;
		}

		std::stringstream stream;
		stream << static_cast<uint32_t>(m_width) << " " << static_cast<uint32_t>(m_height) << " " << static_cast<uint32_t>(m_counts) << std::endl;

		for (size_t i = 0; i < m_counts; i++)
		{
			stream << m_angles[i] << std::endl;
			for (size_t y = 0; y < m_height; y++)
			{
				for (size_t x = 0; x < m_width; x++)
				{
					size_t index = i * m_height * m_width + y * m_width + x;
					stream << m_quantities[index] << " " << static_cast<int32_t>(m_identities[index]) << std::endl;
				}
			}
		}
		file.write(stream.str().c_str(), stream.str().size());
		return true;
	}

	bool empty() const
	{
		return m_quantities.empty() && m_angles.empty() && m_identities.empty();
	}

	void resize(size_t _width, size_t _height, size_t _count)
	{
		m_width = _width;
		m_height = _height;
		m_counts = _count;
		m_angles.resize(_count);
		m_quantities.resize(_count * _height * _width);
		m_identities.resize(m_identities.size());
	}

	size_t width() const
	{
		return m_width;
	}

	size_t height() const
	{
		return m_height;
	}

	size_t counts() const
	{
		return m_counts;
	}

	std::vector<float>& angles()
	{
		return m_angles;
	}

	const std::vector<float>& angles() const
	{
		return m_angles;
	}

	std::vector<float>& quantities()
	{
		return m_quantities;
	}

	const std::vector<float>& quantities() const
	{
		return m_quantities;
	}

	std::vector<blank_type>& identities()
	{
		return m_identities;
	}

	const std::vector<blank_type>& identities() const
	{
		return m_identities;
	}

	float& quantity(size_t _x, size_t _y, size_t _i)
	{
		return m_quantities.at(_x * m_counts * m_height + _i * m_height + _y);
	}

	const float& quantity(size_t _x, size_t _y, size_t _i) const
	{
		return m_quantities.at(_x * m_counts * m_height + _i * m_height + _y);
	}

	blank_type& identity(size_t _x, size_t _y, size_t _i)
	{
		return m_identities.at(_x * m_counts * m_height + _i * m_height + _y);
	}

	const blank_type& identity(size_t _x, size_t _y, size_t _i) const
	{
		return m_identities.at(_x * m_counts * m_height + _i * m_height + _y);
	}

	float& angle(size_t _i)
	{
		return m_angles.at(_i);
	}

	const float& angle(size_t _i) const
	{
		return m_angles.at(_i);
	}

public:
	operator bool() const
	{
		return !empty();
	}

	bool operator!() const
	{
		return empty();
	}

private:
	size_t m_width;
	size_t m_height;
	size_t m_counts;
	std::vector<float> m_angles;
	std::vector<float> m_quantities;
	std::vector<blank_type> m_identities;
};

}