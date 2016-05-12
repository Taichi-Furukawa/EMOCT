#pragma once
#include <cmath>


//#include <cuda_runtime.h>

namespace Math
{

const float Pi	= 3.14159265358979323846f;

template<class Type>  inline Type Max(Type a, Type b)
{
	return a > b ? a : b;
}

template<class Type> inline Type Min(Type a, Type b)
{
	return a < b ? a : b;
}

template<class Type>  inline Type Clamp(Type value, Type minimum, Type maximum)
{
	return Min(Max(value, minimum), maximum);
}

 inline float Hermite(float beginPosition, float beginVelocity, float endPosition, float endVelocity, float factor)
{
	// f(x) = (2Pb+Nb-2Pe+Ne)x^3+(-3Pb-2Nb+3Pe-Ne)x^2+Nbx+Pb
	float factor2 = factor * factor;
	float factor3 = factor * factor2;
	float a = 2.0f * beginPosition + beginVelocity - 2.0f * endPosition + endVelocity;
	float b = -3.0f * beginPosition - 2.0f * beginVelocity + 3.0f * endPosition - endVelocity;
	return a * factor3 + b * factor2 + beginVelocity * factor + beginPosition;
}

/// <summary>�d�݌W�����擾���܂��B</summary>
/// <param name="length">2�_�Ԃ̋���(0.0�`1.0)</param>
/// <returns>�d�݌W��</returns>
 inline float WeightFactor(float length)
{
	return 1.0f - 3.0f * length * length + 2.0f * length * length * length;
}

}
