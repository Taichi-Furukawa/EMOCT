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

#include <cmath>
#include <algorithm>
#include <functional>
#include "distribution.h"
#include "projection.h"

namespace ilab
{

class projector final
{
public:
	projector()
		: m_weight([](float _length){ return 1.0f - _length; })
	{
	}

	projector(const std::function<float(float)>& _weight) noexcept
		: m_weight(_weight)
	{
	}

	std::vector<distribution> calculate_projected_points(const distribution& _distribution, const std::vector<float>& _angles) const
	{
		const float center_x = static_cast<float>(_distribution.width()) / 2.0f;
		const float center_y = static_cast<float>(_distribution.height()) / 2.0f;

		const auto intersect = [](float _x, float _y, float _r)
		{
			return _x * _x + _y * _y <= _r * _r;
		};

		std::vector<distribution> projected_points(_angles.size());
		for (size_t i = 0; i < projected_points.size(); i++)
		{
			const float cos_angle = cos(_angles[i]);
			const float sin_angle = sin(_angles[i]);

			projected_points[i].resize(_distribution);

			for (size_t z = 0; z < _distribution.depth(); z++)
			{
				for (size_t y = 0; y < _distribution.height(); y++)
				{
					for (size_t x = 0; x < _distribution.width(); x++)
					{
						const float direction_x = (static_cast<float>(x) - center_x) * cos_angle;
						const float direction_y = (static_cast<float>(y) - center_y) * sin_angle;
						const float real_point = direction_x + direction_y + center_y;
						const size_t near_point = static_cast<size_t>(real_point);

						if (near_point >= projected_points[i].height() - 1)
						{
							projected_points[i].quantity(x, y, z) = -1.0f;
							projected_points[i].identity(x, y, z) = blank_type::outside;
							continue;
						}

						if (_distribution.identity(x, y, z) != blank_type::quantity || !intersect(static_cast<float>(x) - center_x, static_cast<float>(y) - center_y, center_y))
						{
							projected_points[i].quantity(x, y, z) = -1.0f;
							projected_points[i].identity(x, y, z) = blank_type::outside;
						}
						else
						{
							projected_points[i].quantity(x, y, z) = real_point;
							projected_points[i].identity(x, y, z) = blank_type::quantity;
						}
					}
				}
			}
		}
		return std::move(projected_points);
	}
	
	projection project(const distribution& _distribution, const std::vector<float>& _angles) const
	{
		auto projected_points = calculate_projected_points(_distribution, _angles);
		return std::move(project(_distribution, _angles, projected_points));
	}

	projection project(const distribution& _distribution, const std::vector<float>& _angles, const std::vector<distribution>& _projected_points) const
	{
		projection results(_distribution.depth(), _distribution.height(), _projected_points.size());

		for (size_t i = 0; i < results.counts(); i++)
		{
			results.angle(i) = _angles[i];
			for (size_t z = 0; z < _distribution.depth(); z++)
			{
				for (size_t y = 0; y < _distribution.height(); y++)
				{
					for (size_t x = 0; x < _distribution.width(); x++)
					{
						if (_projected_points[i].identity(x, y, z) != blank_type::quantity)
						{
							continue;
						}

						const float real_point = _projected_points[i].quantity(x, y, z);
						const size_t near_point = static_cast<size_t>(real_point);
						const float distance = 1.0f - (real_point - static_cast<float>(near_point));

						results.quantity(z, near_point + 0, i) += _distribution.quantity(x, y, z) * m_weight(distance);
						results.quantity(z, near_point + 1, i) += _distribution.quantity(x, y, z) * m_weight(1.0f - distance);
					}
				}
			}
		}

		for (size_t i = 0; i < results.counts(); i++)
		{
			for (size_t z = 0; z < _distribution.depth(); z++)
			{
				for (size_t y = 0; y < _distribution.height(); y++)
				{
					for (size_t x = 0; x < _distribution.width(); x++)
					{
						if (_projected_points[i].identity(x, y, z) != blank_type::quantity)
						{
							continue;
						}

						const float real_point = _projected_points[i].quantity(x, y, z);
						const size_t near_point = static_cast<size_t>(real_point);
						const float distance = 1.0f - (real_point - static_cast<float>(near_point));

						if (_distribution.identity(x, y, z) == blank_type::blank)
						{
							results.quantity(z, near_point, i) = 0.0f;
							results.identity(z, near_point, i) = blank_type::blank;
							break;
						}
						results.identity(z, near_point, i) = blank_type::quantity;
					}
				}
			}
		}
		return std::move(results);
	}

	void set_weight(const std::function<float(float)>& _weight)
	{
		m_weight = _weight;
	}

	const std::function<float(float)>& get_weight() const
	{
		return m_weight;
	}

private:
	std::function<float(float)> m_weight;
};

}