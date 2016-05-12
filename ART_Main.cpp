/*
タイムアウトを回避して実行するにはレジストリの
HKEY_LOCAL_MACHINE/SYSTEM/CurrentControlSet/Control/GraphicsDrivers
にDWORD値を追加する必要があります。

http://msdn.microsoft.com/ja-jp/Library/Windows/Hardware/ff569918(v=vs.85).aspx

TdrLevel: REG_DWORD
0:タイムアウトを無効にします。
1:バグチェックのみ行います。
2:VGAへ回復します。(未実装)
3:タイムアウトした場合に回復します。
*/

#include <iostream>
#include <sstream>
#include <string>
#include "ART.h"

using namespace Math;

int main(int argc, char* argv[])
{
	const std::string sectionName = "Constant";
	
	if (argc < 2)
	{
		std::cout << "Arguments[1] : Projection data file path(.prj)" << std::endl;
		return 0;
	}

	std::string data;
	data = argv[1];

	std::cout << "Projection Data : " << data.c_str() << std::endl;
	std::cout << "Loading Projection Data" << std::endl;

	ProjectionData projection(data);
	if(projection.empty())
	{
		std::cout << "Can't load Projection Data" << std::endl;
		return 1;
	}
	std::cout << "Width : " << projection.width << std::endl;
	std::cout << "Height : " << projection.height << std::endl;
	std::cout << "Count : " << projection.count << std::endl;

	DensityDistribution reconstruction;
	unsigned int steps;

	steps = std::stoi(argv[2]);

	std::cout << "Steps : " << steps << std::endl;

	ART art;
	art.setSteps(steps);

	std::cout << "Running ART" << std::endl;
	if(!art.reconstruct(projection, reconstruction))
	{
		std::cout << "reconstruct failed" << std::endl;
		return 1;
	}

	std::cout << "Time : " << art.elapsedTime() << "[sec]" << std::endl;

	const std::string densityFileName("3D-DENSITY");
	if(!reconstruction.save(densityFileName))
	{
		std::cout << "Can't save 3D-Density" << std::endl;
		return 0;
	}

	std::cout << "Complete" << std::endl;
	return 0;
}