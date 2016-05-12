#pragma once

#include <string>
#include <vector>
#include "CUDAMath.h"
#include "ProjectionData.h"
#include "DensityDistribution.h"

class ART
{
public:
	ART();
	~ART();

	/// <summary>���e�f�[�^�����f�[�^���č\�����܂��B</summary>
	/// <returns>���������ꍇ��true, ���s�����ꍇ��false</returns>
	bool reconstruct(const ProjectionData& projections, DensityDistribution& reconstructions);

	void setSteps(unsigned int steps);
	double elapsedTime() const;

private:
	bool isExcluded(const ProjectionData& projections, unsigned int x, unsigned int y, unsigned int z) const;
	float getFittingPosition(unsigned int x, unsigned int y, unsigned int index) const;
	void setReprojection(unsigned int y, unsigned int index, float reprojection);
	void addReprojection(unsigned int y, unsigned int index, float reprojection);
	float getReprojection(unsigned int y, unsigned int index) const;

	void reprojectionPass(const ProjectionData& projections, const DensityDistribution& reconstructions, unsigned int z);
	void feedbackPass(const ProjectionData& projections, unsigned int z);
	void backProjectionPass(const ProjectionData& projections, DensityDistribution& reconstructions, unsigned int z);
	void gaussian(DensityDistribution& reconstructions, unsigned int z);

	unsigned int		m_DimensionX;
	unsigned int		m_DimensionY;
	unsigned int		m_DimensionZ;
	unsigned int		m_DimensionI;
	unsigned int		m_Steps;
	double				m_ElapsedTime;
	std::vector<float>	m_FittingPositions;
	std::vector<float>	m_Reprojections;
	std::vector<float>	m_SmoothingReconstructions;
};
