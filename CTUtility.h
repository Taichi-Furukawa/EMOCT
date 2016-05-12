#pragma once
#include "DensityDistribution.h"
#include "ProjectionData.h"

std::vector<float> CreateFittingPositions(const ProjectionData& projections);
ProjectionData CreateProjectionDataAtX(const DensityDistribution& data, const std::vector<float>& angles);