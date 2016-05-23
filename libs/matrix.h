//=================================================================================//
//                                                                                 //
//  ArchEngine                                                                     //
//                                                                                 //
//  Copyright (C) 2011-2016 Terry                                                  //
//                                                                                 //
//  This file is a portion of the ArchEngine. It is distributed under the MIT      //
//  License, available in the root of this distribution and at the following URL.  //
//  http://opensource.org/licenses/mit-license.php                                 //
//                                                                                 //
//=================================================================================//

#pragma once

#include <array>
#include <initializer_list>

namespace arch
{

template<class T, size_t row_size, size_t column_size>
class matrix
{
public:
	typedef typename std::array<T, row_size * column_size>::value_type value_type;
	typedef typename std::array<T, row_size * column_size>::size_type size_type;
	typedef typename std::array<T, row_size * column_size>::pointer pointer;
	typedef typename std::array<T, row_size * column_size>::const_pointer const_pointer;
	typedef typename std::array<T, row_size * column_size>::reference reference;
	typedef typename std::array<T, row_size * column_size>::const_reference const_reference;
	typedef typename std::array<T, row_size * column_size>::iterator iterator;
	typedef typename std::array<T, row_size * column_size>::const_iterator const_iterator;
	typedef typename std::array<T, row_size * column_size>::reverse_iterator reverse_iterator;
	typedef typename std::array<T, row_size * column_size>::const_reverse_iterator const_reverse_iterator;

public:
	matrix() = default;
	
	constexpr matrix(const std::initializer_list<value_type>& list)
	: m_elements()
	{
		auto itr = list.begin();
		for (auto& value : m_elements)
		{
			if (itr == list.end())
			{
				break;
			}
			value = *itr;
			itr++;
		}
	}
	~matrix() = default;

	iterator begin()
	{
		return m_elements.begin();
	}

	const_iterator begin() const
	{
		return m_elements.begin();
	}

	iterator end()
	{
		return m_elements.end();
	}

	const_iterator end() const
	{
		return m_elements.end();
	}

	value_type& at(size_type y, size_type x)
	{
		return m_elements[y * column_size + x];
	}

	const value_type& at(size_type y, size_type x) const
	{
		return m_elements[y * column_size + x];

	}

	void fill(const value_type& value)
	{
		m_elements.fill(value);
	}

	///	<summary>行数を取得します。</summary>
	constexpr size_type rows() const
	{
		return row_size;
	}
	
	///	<summary>列数を取得します。</summary>
	constexpr size_type columns() const
	{
		return column_size;
	}

	constexpr size_type size() const
	{
		return m_elements.size();
	}

	value_type* data()
	{
		return m_elements.data();
	}

	const value_type* data() const
	{
		return m_elements.data();
	}

	const value_type* operator[](size_type y) const
	{
		return &m_elements[y * column_size];
	}

	value_type* operator[](size_type y)
	{
		return &m_elements[y * column_size];
	}

	bool is_symmetric() const
	{
		for (size_type y = 0; y < row_size; y++)
		{
			for (size_type x = 0; x < column_size; x++)
			{
				if (at(y, x) != at(x, y))
				{
					return false;
				}
			}
		}
		return true;
	}

public:
	static constexpr matrix zero()
	{
		matrix<value_type, row_size, column_size> mat;
		for (size_type y = 0; y < row_size; y++)
		{
			for (size_type x = 0; x < column_size; x++)
			{
				mat[y][x] = static_cast<value_type>(0.0);
			}
		}
		return std::move(mat);
	}

	static constexpr matrix identity()
	{
		static_assert(row_size != column_size, "row size is not equal column size.");
		matrix<value_type, row_size, column_size> mat;
		for (size_type y = 0; y < row_size; y++)
		{
			for (size_type x = 0; x < column_size; x++)
			{
				mat[y][x] = static_cast<value_type>(y == x ? 1.0 : 0.0);
			}
		}
		return std::move(m);
	}
	
private:
	std::array<value_type, row_size * column_size> m_elements;
};

}