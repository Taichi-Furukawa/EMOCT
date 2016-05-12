#include "CTUtility.h"
#include "CUDAMath.h"
#include <iterator>

std::vector<float> CreateFittingPositions(const ProjectionData& projections)
{
	std::vector<float> fittingPositions(projections.height * projections.height * projections.count);

	// 除外領域と投影位置を算出
	const int center			= projections.height / 2;		// 中心位置
	const int maxRadius			= projections.height - center;	// 最大半径
	const int squaredMaxRadius	= maxRadius * maxRadius;		// 最大半径の2乗(斜辺の長さの2乗)

	for (unsigned int y = 0; y < projections.height; y++)
	{
		for (unsigned int x = 0; x < projections.height; x++)
		{
			int distanceX = x - center;												// X軸方向における現在位置と中心との距離(半径)
			int distanceY = y - center;												// Y軸方向における現在位置と中心との距離(半径)
			int squaredDistance = distanceX * distanceX + distanceY * distanceY;	// 中心位置から現在の位置までの距離の2乗

			// 再構成領域は円形なので現在位置が領域内に収まっているかを判断
			if (squaredDistance <= squaredMaxRadius)
			{
				for (unsigned int i = 0; i < projections.count; i++)
				{
					int fittingIndex = i * projections.height * projections.height + y * projections.height + x;
					float positionX = static_cast<float>(distanceX) * cos(projections.angles[i].roll);					// ある投影角における位置Xから再構成領域における位置Xを算出
					float positionY = static_cast<float>(distanceY) * sin(projections.angles[i].roll);					// ある投影角における位置Yから再構成領域における位置Yを算出
					fittingPositions[fittingIndex] = positionX + positionY + static_cast<float>(center);				// 投影データのフィッティング位置

					// フィッティング位置が領域外(0.0〜m_DimensionY-1が1つの投影データの領域)の場合は範囲外に設定
					if (fittingPositions[fittingIndex] < 0.0f || fittingPositions[fittingIndex] > static_cast<float>(projections.height - 1))	
					{
						fittingPositions[fittingIndex] = -1.0f;
					}
				}
			}
			else
			{
				for (unsigned int i = 0; i < projections.count; i++)
				{
					int fittingIndex = i * projections.height * projections.height + y * projections.height + x;
					fittingPositions[fittingIndex] = -1.0f;
				}
			}
		}
	}

	return std::move(fittingPositions);
}

ProjectionData CreateProjectionDataAtX(const DensityDistribution& data, const std::vector<float>& angles)
{
	ProjectionData projections;
	projections.resize(data.depth, data.width, angles.size());

	for (unsigned int i = 0; i < projections.angles.size(); i++)
	{
		projections.angles[i] = projections.angles[i];
	}

	for (auto& projection : projections.data)
	{
		projection = 0.0;
	}

	auto fittingPositions = CreateFittingPositions(projections);

	for (unsigned int z = 0; z < projections.width; z++)
	{
		for (unsigned int y = 0; y < projections.height; y++)
		{
			for (unsigned int x = 0; x < projections.height; x++)
			{
				for (unsigned int i = 0; i < projections.angles.size(); i++)
				{
					float position		= fittingPositions[i * projections.height * projections.height + y * projections.height + x];
					int fittingIndex	= static_cast<int>(position);
					float distance		= position - static_cast<float>(fittingIndex);
					int projectionIndex	= i * projections.height + fittingIndex;
					float beginWeight	= Math::WeightFactor(distance);
					float endWeight		= Math::WeightFactor(1.0f - distance);

					projections.data[projectionIndex] += data.at(x, y, z) * beginWeight;
					projections.data[projectionIndex + 1] += data.at(x, y, z) * endWeight;
				}
			}
		}
	}

	return std::move(projections);
}

