#pragma once

#include <string>
#include <vector>
#include <arch.h>
#include "ProjectionData.h"
#include "DensityDistribution.h"

/// <summary>���e�f�[�^����f�[�^���č\�����܂��B</summary>
class art
{
public:
	art()
		: m_dimension(arch::vector3<size_t>::zero()), m_cameras(0), m_steps(30)
	{
	}

	~art() = default;

	/// <summary>���e�f�[�^����f�[�^��Z������ōč\�����܂��B</summary>
	/// <returns>���������ꍇ��true, ���s�����ꍇ��false</returns>
	bool reconstruct(const arch::table<float>& _projection, arch::grid<float>& _reconstruction)
	{
		bool result = false;
		m_DimensionX = projections.height;
		m_DimensionY = projections.height;
		m_DimensionZ = projections.width;
		m_DimensionI = projections.count;

		// �K�v�ȃo�b�t�@���m��
		reconstructions.resize(m_DimensionX, m_DimensionY, m_DimensionZ);
		m_SmoothingReconstructions.resize(m_DimensionX * m_DimensionY);
		m_Reprojections.resize(m_DimensionY * m_DimensionI);

		// �ē��e�ʒu���Z�o
		m_FittingPositions = CreateFittingPositions(projections);

		// �č\���J�n���Ԃ��擾
		time_t beginTime;
		time_t endTime;
		time(&beginTime);

		// �č\���J�n
		for (unsigned int z = 0; z < m_DimensionZ; z++)
		{
			for (unsigned int step = 0; step < m_Steps; step++)
			{
				reprojectionPass(projections, reconstructions, z);		// �ē��e
				feedbackPass(projections, z);							// �t�B���^�����O
				backProjectionPass(projections, reconstructions, z);	// �t���e
				gaussian(reconstructions, z);							// �K�E�X�t�B���^�ŕ�����
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