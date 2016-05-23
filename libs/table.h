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

#include <vector>
#include <initializer_list>

namespace arch
{

template<class T>
class table
{
public:
	typedef typename std::vector<T>::value_type value_type;
	typedef typename std::vector<T>::size_type size_type;
	typedef typename std::vector<T>::pointer pointer;
	typedef typename std::vector<T>::const_pointer const_pointer;
	typedef typename std::vector<T>::reference reference;
	typedef typename std::vector<T>::const_reference const_reference;
	typedef typename std::vector<T>::iterator iterator;
	typedef typename std::vector<T>::const_iterator const_iterator;
	typedef typename std::vector<T>::reverse_iterator reverse_iterator;
	typedef typename std::vector<T>::const_reverse_iterator const_reverse_iterator;

public:
	table()
	: m_row_size(0), m_column_size(0)
	{
	}
	
	table(size_type row_size, size_type column_size)
	: m_row_size(row_size), m_column_size(column_size), m_elements(row_size * column_size)
	{
	}
	
	table(size_type row_size, size_type column_size, value_type value)
	: m_row_size(row_size), m_column_size(column_size), m_elements(row_size * column_size, value)
	{
	}
	
	~table() = default;

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
		return m_elements[y * m_column_size + x];
	}

	const value_type& at(size_type y, size_type x) const
	{
		return m_elements[y * m_column_size + x];

	}

	void fill(const value_type& _value)
	{
		std::fill(m_elements.begin(), m_elements.end(), _value);
	}

	void assign(size_type row_size, size_type column_size)
	{
		m_elements.assign(row_size * column_size);
		m_row_size = row_size;
		m_column_size = column_size;
	}

	void assign(size_type row_size, size_type column_size, value_type value)
	{
		m_elements.assign(row_size * column_size, value);
		m_row_size = row_size;
		m_column_size = column_size;
	}

	void resize(size_type row_size, size_type column_size)
	{
		std::vector<value_type> elements(row_size * column_size);
		for (size_type y = 0; y < m_row_size; y++)
		{
			for (size_type x = 0; x < m_column_size; x++)
			{
				size_type index = y * m_column_size + x;
				m_elements[index] = elements[index];
			}
		}

		m_elements = std::move(elements);
		m_row_size = row_size;
		m_column_size = column_size;
	}

	void resize(size_type row_size, size_type column_size, value_type value)
	{
		std::vector<value_type> elements(row_size * column_size, value);
		for (size_type y = 0; y < m_row_size; y++)
		{
			for (size_type x = 0; x < m_column_size; x++)
			{
				size_type index = y * m_column_size + x;
				m_elements[index] = elements[index];
			}
		}

		m_elements = std::move(elements);
		m_row_size = row_size;
		m_column_size = column_size;
	}


	///	<summary>行数を取得します。</summary>
	size_type rows() const
	{
		return m_row_size;
	}
	
	///	<summary>列数を取得します。</summary>
	size_type columns() const
	{
		return m_column_size;
	}

	size_type size() const
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
		return &m_elements[y * m_column_size];
	}

	value_type* operator[](size_type y)
	{
		return &m_elements[y * m_column_size];
	}

private:
	std::vector<value_type> m_elements;
	size_t m_row_size;
	size_t m_column_size;
};

template<typename value_type, typename char_type> inline std::basic_ostream<char_type>& operator <<(std::basic_ostream<char_type>& _os, const table<value_type>& _table)
{
	_os << '(';
	for (typename table<value_type>::size_type y = 0; y < _table.rows(); y++)
	{
		_os << '(';
		for (typename table<value_type>::size_type x = 0; x < _table.columns(); x++)
		{
			_os << _table.at(y, x);
			if (x != _table.columns() - 1)
			{
				_os << ',';
			}
		}
		_os << ')';
		if (y != _table.rows() - 1)
		{
			_os << ','<<std::endl;
		}
	}
	return _os << ')';
}

}