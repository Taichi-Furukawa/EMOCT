#pragma once

#include <string>
#include <vector>
#include <arch.h>
#include "ProjectionData.h"
#include "DensityDistribution.h"

/// <summary>投影データからデータを再構成します。</summary>
class art
{
public:
	art()
		: m_dimension(arch::vector3<size_t>::zero()), m_cameras(0), m_steps(30)
	{
	}

	~art() = default;

	/// <summary>投影データからデータをZ軸周りで再構成します。</summary>
	/// <returns>成功した場合はtrue, 失敗した場合はfalse</returns>
	bool reconstruct(const arch::table<float>& _projection, arch::grid<float>& _reconstruction)
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

	void set_steps(size_t _steps)
	{
		m_steps = _steps;
	}

	size_t get_steps() const
	{
		return m_steps;
	}

private:
	bool is_excluded(const ProjectionData& projections, unsigned int x, unsigned int y, unsigned int z) const
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

	float getFittingPosition(unsigned int x, unsigned int y, unsigned int index) const;
	void setReprojection(unsigned int y, unsigned int index, float reprojection);
	void addReprojection(unsigned int y, unsigned int index, float reprojection);
	float getReprojection(unsigned int y, unsigned int index) const;

	void reprojectionPass(const ProjectionData& projections, const DensityDistribution& reconstructions, unsigned int z);
	void feedbackPass(const ProjectionData& projections, unsigned int z);
	void backProjectionPass(const ProjectionData& projections, DensityDistribution& reconstructions, unsigned int z);
	void gaussian(DensityDistribution& reconstructions, unsigned int z);

	arch::vector3<size_t> m_dimension;
	size_t m_cameras;
	size_t m_steps;

	arch::grid<float> m_fitting_positions;
	arch::table<float> m_Reprojections;
	std::vector<float> m_SmoothingReconstructions;
};