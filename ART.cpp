#include <iostream>
#include <array>
#include <ctime>
#include "ART.h"
#include "Xorshift.h"
#include "MersenneTwister.h"
#include "CTUtility.h"

ART::ART()
: m_DimensionX(0),
m_DimensionY(0),
m_DimensionZ(0),
m_DimensionI(0),
m_Steps(30),
m_ElapsedTime(0.0)
{
}

ART::~ART()
{
}

void ART::setSteps(unsigned int steps)
{
	m_Steps = steps;
}

bool ART::reconstruct(const ProjectionData& projections, DensityDistribution& reconstructions)
{
	bool result = false;
	m_DimensionX = projections.height;
	m_DimensionY = projections.height;
	m_DimensionZ = projections.width;
	m_DimensionI = projections.count;

	// 必要なバッファを確保
	reconstructions.resize(m_DimensionX, m_DimensionY, m_DimensionZ);
	m_SmoothingReconstructions.resize(m_DimensionX * m_DimensionY);
	m_Reprojections.resize(m_DimensionY * m_DimensionI);

	// 再投影位置を算出
	m_FittingPositions = CreateFittingPositions(projections);

	// 再構成開始時間を取得
	time_t beginTime;
	time_t endTime;
	time(&beginTime);

	// 再構成開始
	for (unsigned int z = 0; z < m_DimensionZ; z++)
	{
		for (unsigned int step = 0; step < m_Steps; step++)
		{
			reprojectionPass(projections, reconstructions, z);		// 再投影
			feedbackPass(projections, z);							// フィルタリング
			backProjectionPass(projections, reconstructions, z);	// 逆投影
			gaussian(reconstructions, z);							// ガウスフィルタで平滑化
		}
		time(&endTime);
		m_ElapsedTime = difftime(endTime, beginTime);

		std::cout << "z=" << z << " exit (elapsed time " << m_ElapsedTime << "[sec])" << std::endl;
	}
	return true;
}

double ART::elapsedTime() const
{
	return m_ElapsedTime;
}

float ART::getFittingPosition(unsigned int x, unsigned int y, unsigned int index) const
{
	return m_FittingPositions[index * m_DimensionY * m_DimensionX + y * m_DimensionX + x];
}

bool ART::isExcluded(const ProjectionData& projections, unsigned int x, unsigned int y, unsigned int z) const
{
	int total = 0;
	for (unsigned int i = 0; i < projections.count; i++)
	{
		float position = getFittingPosition(x, y, i);
		int fittingIndex = static_cast<int>(position);

		total += (fittingIndex < 0) ? 1 : projections.id(z, fittingIndex, i);
	}

	return static_cast<unsigned int>(total) >= projections.count;
}

void ART::setReprojection(unsigned int y, unsigned int index, float reprojection)
{
	m_Reprojections[index * m_DimensionY + y] = reprojection;
}

float ART::getReprojection(unsigned int y, unsigned int index) const
{
	return m_Reprojections[index * m_DimensionY + y];
}

void ART::addReprojection(unsigned int y, unsigned int index, float reprojection)
{
	m_Reprojections[index * m_DimensionY + y] += reprojection;
}

void ART::reprojectionPass(const ProjectionData& projections, const DensityDistribution& reconstructions, unsigned int z)
{
	std::fill(m_Reprojections.begin(), m_Reprojections.end(), 0.0f);

	for (unsigned int i = 0; i < m_DimensionI; i++)
	{
		for (unsigned int y = 0; y < m_DimensionY; y++)
		{
			for (unsigned int x = 0; x < m_DimensionX; x++)
			{
				if (!isExcluded(projections, x, y, z))
				{
					float projection		= getFittingPosition(x, y, i);
					int fittingIndex		= static_cast<int>(projection);
					float distance			= projection - static_cast<float>(fittingIndex);
					float reconstruction	= reconstructions.at(x, y, z);

					addReprojection(fittingIndex, i, reconstruction * Math::WeightFactor(distance));
					addReprojection(fittingIndex + 1, i, reconstruction * Math::WeightFactor(1.0f - distance));
				}
			}
		}	
	}
}

void ART::feedbackPass(const ProjectionData& projections, unsigned int z)
{
	for (unsigned int i = 0; i < m_DimensionI; i++)
	{
		for (unsigned int y = 0; y < m_DimensionY; y++)
		{
			setReprojection(y, i, projections.at(z, y, i) - getReprojection(y, i));
		}
	}
}

void ART::backProjectionPass(const ProjectionData& projections, DensityDistribution& reconstructions, unsigned int z)
{
	float factor = static_cast<float>(m_DimensionX * m_DimensionI) / 2.0f;

	for (unsigned int i = 0; i < m_DimensionI; i++)
	{
		for (unsigned int y = 0; y < m_DimensionY; y++)
		{
			for (unsigned int x = 0; x < m_DimensionX; x++)
			{
				float position			= getFittingPosition(x, y, i);
				int fittingIndex		= static_cast<int>(position);
				if(fittingIndex==-1){
					continue;
				}
				float distance			= position - static_cast<float>(fittingIndex);
				float beginWeight		= Math::WeightFactor(distance);
				float endWeight			= Math::WeightFactor(1.0f - distance);
				int currentID			= projections.id(z, fittingIndex, i);
				int nextID				= projections.id(z, fittingIndex + 1, i);

				if (!(currentID || nextID))
				{
					// 現在のIDと次のIDが両方0(除外しない)の場合は現在位置の投影データと次の投影データで補間する
					reconstructions.at(x, y, z) += (getReprojection(fittingIndex, i) * beginWeight + getReprojection(fittingIndex + 1, i) * endWeight) / factor;
				}
				else if (!currentID && nextID)
				{
					// 現在のIDが0, 次のIDが1の場合は現在位置の投影データと前の投影データで次のデータを補間する
					float gradient = getReprojection(fittingIndex, i) - getReprojection(fittingIndex - 1, i);
					reconstructions.at(x, y, z) += (getReprojection(fittingIndex, i) * beginWeight + (getReprojection(fittingIndex, i) + gradient) * endWeight) / factor;
				}
				else if (currentID && !nextID)
				{
					// 現在のIDが1で次のIDが0の場合
					float gradient = getReprojection(fittingIndex + 2, i) - getReprojection(fittingIndex + 1, i);
					reconstructions.at(x, y, z) += ((getReprojection(fittingIndex + 1, i) - gradient) * beginWeight + getReprojection(fittingIndex + 1, i) * endWeight) / factor;
				}
			}
		}
	}
}

void ART::gaussian(DensityDistribution& reconstructions, unsigned int z)
{
	const std::array<float, 9> weights = 
	{
		1.0f, 2.0f, 1.0f,
		2.0f, 4.0f, 2.0f,
		1.0f, 2.0f, 1.0f,
	};
	const float total = 16.0f;

	std::fill(m_SmoothingReconstructions.begin(), m_SmoothingReconstructions.end(), 0.0f);

	for (unsigned int y = 1; y < m_DimensionY - 1; y++)
	{
		for (unsigned int x = 1; x < m_DimensionX - 1; x++)
		{
			int index = y * m_DimensionX + x;

			for (unsigned int dy = 0; dy < 3; dy++)
			{
				for (unsigned int dx = 0; dx < 3; dx++)
				{
					m_SmoothingReconstructions.at(index) += reconstructions.at(x + dx - 1, y + dy - 1, z) * weights[dy * 3 + dx] / total;
				}
			}
		}
	}

	for (unsigned int y = 1; y < m_DimensionY - 1; y++)
	{
		for (unsigned int x = 1; x < m_DimensionX - 1; x++)
		{
			int index = y * m_DimensionX + x;
			reconstructions.at(x, y, z) = m_SmoothingReconstructions.at(index);
		}
	}
}
