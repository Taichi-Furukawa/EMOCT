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
#include "ilab.h"
#include "EMO.h"

using namespace Math;
using namespace ilab;

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

	ilab::projection projection(data);
	if(projection.empty())
	{
		std::cout << "Can't load Projection Data" << std::endl;
		return 1;
	}
	std::cout << "Width : " << projection.width() << std::endl;
	std::cout << "Height : " << projection.height() << std::endl;
	std::cout << "Count : " << projection.counts() << std::endl;

	distribution reconstruction;
	unsigned int steps;

	std::cout << "Steps : " << steps << std::endl;

	std::string dist_data;
	if (argc == 3)
	{
		std::cout << "GOT Arguments[2] : Initial Reconstruct Image file path(distribution file)" << std::endl;
		dist_data = argv[2];
	}

	//GACT gact(projection,dist_data);//call GACT
	//gact.Evolution();
    EMO emo(projection);//call EMO

    /*
	const std::string densityFileName("3D-DENSITY");
	if(!reconstruction.save(densityFileName))
	{
		std::cout << "Can't save 3D-Density" << std::endl;
		return 0;
	}
    */

	std::cout << "Complete" << std::endl;
	return 0;
}