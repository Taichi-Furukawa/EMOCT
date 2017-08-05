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
#include <omp.h>
#include "ilab.h"
#include "EMO.h"

using namespace Math;
using namespace ilab;

/*
int get_irho() {
    constexpr float rho0 = 1.0f;
    constexpr size_t cameras = 4;
    constexpr float dtheta = 180.0f / static_cast<float>(cameras);

    const std::string model_path = "experiment_data/no_object(196,196,1).cfd";
    const std::string output_path = "pdata.prj";

    // Calculate projection angles.
    std::vector<float> angles(cameras);
    for (size_t i = 0; i < angles.size(); i++)
    {
        angles[i] = dtheta * static_cast<float>(i);
    }

    // Load the density distribution;
    ilab::distribution rho(model_path);
    if (rho.empty())
    {
        return -1;
    }

    // Calculate the dimensionless number.
    for (auto& r : rho.quantities())
    {
        r = (r - rho0) / rho0;
    }

    // Calculate the blanks.
    const auto intersect = [](float _x, float _y, float _r)
    {
        return _x * _x + _y * _y <= _r * _r;
    };

    const float center_x = static_cast<float>(rho.width()) / 2.0f;
    const float center_y = static_cast<float>(rho.height()) / 2.0f;
    for (size_t z = 0; z < rho.depth(); z++)
    {
        for (size_t y = 0; y < rho.height(); y++)
        {
            for (size_t x = 0; x < rho.width(); x++)
            {
                if (intersect(static_cast<float>(x) - center_x, static_cast<float>(y) - center_y, center_y))
                {
                    rho.identity(x, y, z) = ilab::blank_type::quantity;
                }
                else
                {
                    rho.identity(x, y, z) = ilab::blank_type::outside;
                }
            }
        }
    }

    // Calculate the projection data.
    const auto weight_function = [](float _length)
    {
        //return 1.0 - _length;
        return 1.0f - 3.0f * _length * _length + 2.0f * _length * _length * _length;
    };
    ilab::projector projector(weight_function);
    ilab::projection irho = projector.project(rho, angles);
    for(auto &r:irho.angles()){
        r = r* static_cast<float>(M_PI/180.0);
    }
    // Save the result at the PRJ format file.
    irho.save(output_path);
    InverseDomain inv(irho);
    inv.save_notshift("new_p_data.png");

    return 0;
}
 */

projection miss_pdata_center(projection p,int r,int saveflag){
    const int center_y = static_cast<int>(p.height()) / 2;
    for(int i=0;i<p.counts();i++){
        for(int t=center_y-r;t<center_y+r;t++){
            p.quantity(0, static_cast<size_t>(t), static_cast<size_t>(i)) = 0.0f;

        }
    }
    if(saveflag==1){
        p.save("p_data/no_object_missing("+to_string(p.height())+","+to_string(p.height())+","+to_string(p.width())+")_"+to_string(p.counts())+".prj");
    }
    return p;
}


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

#pragma omp parallel
    {
        // ここがコア数分だけ並列に実行される。１コアだと１つです。
        printf("Hello, World ! %d of %d\n", omp_get_thread_num(), omp_get_num_threads());
    }

    //get_irho();
    EMO emo(projection);//call EMO
    emo.evolution();
    //projection = miss_pdata_center(projection,20,0);
    //GACT gact(projection,dist_data);//call GACT
    //gact.Evolution();
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