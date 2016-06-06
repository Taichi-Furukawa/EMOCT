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

#include "projection.h"
#include "projector.h"
#include "distribution.h"

namespace ilab
{

/// <summary>投影データから分布データを再構成します。</summary>
class art final
{
public:
	art() = default;
	~art() = default;

	/// <summary>投影データからデータをZ軸周りで再構成します。</summary>
	distribution reconstruct(const projection& _projection, const projector& _projector, size_t _steps)
	{
		distribution result(_projection.height(), _projection.height(), _projection.width(), 0.0f);

		// Calculate the projected points;
		const auto projected_points = _projector.calculate_projected_points(result, _projection.angles());

		projection correction(_projection.width(), _projection.height(), _projection.counts());

		for (size_t step = 0; step < _steps; step++)
		{
			for (size_t z = 0; z < result.depth(); z++)
			{
				// Reproject.
				projection reprojection = _projector.project(result, _projection.angles(), projected_points);
				
				// Correct.
				for (size_t i = 0; i < reprojection.quantities().size(); i++)
				{
					if (reprojection.identities()[i] == blank_type::quantity)
					{
						correction.quantities()[i] = (_projection.quantities()[i] - reprojection.quantities()[i]) / static_cast<float>(correction.height());
					}
				}

				// Backproject.
				result.fill(0.0f);
				for (size_t y = 0; y < result.height(); y++)
				{
					for (size_t x = 0; x < result.width(); x++)
					{
						float total = 0.0f;
						for (size_t i = 0; i < _projection.counts(); i++)
						{
							if (projected_points[i].identity(x, y, z) != blank_type::quantity)
							{
								continue;
							}

							const float real_point = projected_points[i].quantity(x, y, z);
							const size_t near_point = static_cast<size_t>(real_point);
							const float distance = 1.0f - (real_point - static_cast<float>(near_point));

							if (_projection.identity(z, near_point, i) == blank_type::quantity)
							{
								result.quantity(x, y, z) += correction.quantity(z, near_point + 0, i) * _projector.get_weight()(distance);
								result.quantity(x, y, z) += correction.quantity(z, near_point + 1, i) * _projector.get_weight()(1.0f - distance);
								total++;
							}
						}
						result.quantity(x, y, z) /= total;
					}
				}
			}
		}
		return std::move(result);
	}
};

}