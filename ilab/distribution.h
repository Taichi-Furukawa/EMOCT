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

#include <algorithm>
#include <vector>
#include <string>
#include <fstream>

namespace ilab
{

enum class blank_type : int32_t
{
	quantity = 0,
	blank = 1,
	outside = 2,
};

/// <summary>分布データ(*.cfd)ファイルを扱います</summary>
class distribution final
{
public:
	distribution()
		: m_width(0), m_height(0), m_depth(0)
	{
	}

	distribution(size_t _width, size_t _height, size_t _depth)
		: m_width(_width), m_height(_height), m_depth(_depth), m_quantities(_width * _height * _depth), m_identities(_width * _height * _depth, blank_type::quantity)
	{
	}
	
	distribution(size_t _width, size_t _height, size_t _depth, float _value)
		: m_width(_width), m_height(_height), m_depth(_depth), m_quantities(_width * _height * _depth, _value), m_identities(_width * _height * _depth, blank_type::quantity)
	{
	}

	distribution(const std::string& _path, bool _is_binary = true)
		: m_width(0), m_height(0), m_depth(0)
	{
		std::ifstream file(_path, _is_binary ? std::ios::in | std::ios::binary : std::ios::in);
		if (!file)
		{
			// ファイルを開けなかった場合はエラー
			return;
		}

		uint32_t width, height, depth, n;
		file.read(reinterpret_cast<char*>(&width), sizeof(uint32_t));
		file.read(reinterpret_cast<char*>(&height), sizeof(uint32_t));
		file.read(reinterpret_cast<char*>(&depth), sizeof(uint32_t));
		file.read(reinterpret_cast<char*>(&n), sizeof(uint32_t));

		m_width = static_cast<size_t>(width);
		m_height = static_cast<size_t>(height);
		m_depth = static_cast<size_t>(depth);
		m_quantities.resize(m_width * m_height * m_depth);
		m_identities.resize(m_width * m_height * m_depth, blank_type::quantity);

		for (auto& quantity : m_quantities)
		{
			file.read(reinterpret_cast<char*>(&quantity), sizeof(float));
		}

		if (file.eof())
		{
			return;
		}

		for (auto& identity : m_identities)
		{
			file.read(reinterpret_cast<char*>(&identity), sizeof(int32_t));
		}
	}

	~distribution() = default;

	bool save(const std::string& _path, bool _is_binary = true) const
	{
		std::ofstream file(_path, _is_binary ? std::ios::out | std::ios::trunc | std::ios::binary : std::ios::out | std::ios::trunc);
		if (!file)
		{
			// ファイルを開けなかった場合はエラー
			return false;
		}

		const uint32_t width = static_cast<uint32_t>(m_width);
		const uint32_t height = static_cast<uint32_t>(m_height);
		const uint32_t depth = static_cast<uint32_t>(m_depth);
		const uint32_t n = 1;

		file.write(reinterpret_cast<const char*>(&width), sizeof(uint32_t));
		file.write(reinterpret_cast<const char*>(&height), sizeof(uint32_t));
		file.write(reinterpret_cast<const char*>(&depth), sizeof(uint32_t));
		file.write(reinterpret_cast<const char*>(&n), sizeof(uint32_t));
	
		for (auto& quantity : m_quantities)
		{
			file.write(reinterpret_cast<const char*>(&quantity), sizeof(float));
		}

		for (auto& identity : m_identities)
		{
			file.write(reinterpret_cast<const char*>(&identity), sizeof(int32_t));
		}

		return true;
	}

	distribution copy() const
	{
		distribution result;
		result.resize(*this);
		std::copy(m_quantities.begin(), m_quantities.end(), result.m_quantities.begin());
		std::copy(m_identities.begin(), m_identities.end(), result.m_identities.begin());
		return std::move(result);
	}

	void fill(float _value)
	{
		std::fill(m_quantities.begin(), m_quantities.end(), _value);
	}

	void resize(size_t _width, size_t _height, size_t _depth)
	{
		m_width = _width;
		m_height = _height;
		m_depth = _depth;
		m_quantities.resize(_width * _height * _depth);
		m_identities.resize(_width * _height * _depth);
	}

	void resize(const distribution& _distribution)
	{
		m_width = _distribution.m_width;
		m_height = _distribution.m_height;
		m_depth = _distribution.m_depth;
		m_quantities.resize(_distribution.m_quantities.size());
		m_identities.resize(_distribution.m_identities.size());
	}

	size_t width() const
	{
		return m_width;
	}

	size_t height() const
	{
		return m_height;
	}

	size_t depth() const
	{
		return m_depth;
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

	bool empty() const
	{
		return m_quantities.empty() && m_identities.empty();
	}

	float& quantity(size_t _x, size_t _y, size_t _z)
	{
		return m_quantities.at(_z * m_width * m_height + _y * m_width + _x);
	}

	const float& quantity(size_t _x, size_t _y, size_t _z) const
	{
		return m_quantities.at(_z * m_width * m_height + _y * m_width + _x);
	}

	blank_type& identity(size_t _x, size_t _y, size_t _z)
	{
		return m_identities.at(_z * m_width * m_height + _y * m_width + _x);
	}

	const blank_type& identity(size_t _x, size_t _y, size_t _z) const
	{
		return m_identities.at(_z * m_width * m_height + _y * m_width + _x);
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
	size_t m_depth;
	std::vector<float> m_quantities;
	std::vector<blank_type> m_identities;
};

}