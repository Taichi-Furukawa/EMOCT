#pragma once

#include <random>

namespace MersenneTwister
{

static std::mt19937 g_RandomEngine;

inline unsigned long Random()
{
	return g_RandomEngine();
}

template <class Type> inline Type RandomInt(Type Min, Type Max)
{
	unsigned long Range = Max - Min + 1;
	return (Random() % Range) + Min;
}

template <class Type> inline Type RandomReal(Type minimum, Type maximum)
{
	Type number	= static_cast<Type>(Random()) / static_cast<Type>(g_RandomEngine.max());
	return number * (maximum - minimum) + minimum;
}

}