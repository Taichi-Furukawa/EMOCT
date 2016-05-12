#pragma once

namespace Xorshift
{

static unsigned long g_RandomMax	= ULONG_MAX;

inline unsigned long Random()
{
	static unsigned long xors_x		= 123456789;
	static unsigned long xors_y		= 362436069;
	static unsigned long xors_z		= 521288629;
	static unsigned long xors_w		= 88675123;

	unsigned long t = (xors_x ^(xors_x << 11));
	xors_x	= xors_y;
	xors_y	= xors_z;
	xors_z	= xors_w;
	xors_w	= (xors_w ^(xors_w >> 19)) ^(t ^(t >> 8));
	return xors_w;
}

template <class Type> inline Type RandomInt(Type Min, Type Max)
{
	unsigned long Range = Max - Min + 1;
	return (Random() % Range) + Min;
}

template <class Type> inline Type RandomReal(Type minimum, Type maximum)
{
	Type number	= static_cast<Type>(Random()) / static_cast<Type>(g_RandomMax);
	return number * (maximum - minimum) + minimum;
}


}