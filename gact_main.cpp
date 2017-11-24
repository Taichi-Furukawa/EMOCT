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

void dist_normalize(distribution dist){
    float max=FLT_MIN;
    float min=FLT_MAX;
    int begin = 0;
    int end = 255;
    for(int i=0;i<dist.quantities().size();i++){
        if(max<dist.quantities()[i]){
            max = dist.quantities()[i];
        }
        if(min>dist.quantities()[i]){
            min = dist.quantities()[i];
        }
    }
    for(int i=0;i<dist.quantities().size();i++){
        dist.quantities()[i] = dist.quantities()[i]*(max - min)/ 255;
    }
    dist.save("normalize");

}

void save_newp_data(projection p_data){
    constexpr float rho0 = 1.0f;
    size_t cameras = p_data.counts();
    float dtheta = 180.0f / static_cast<float>(cameras);
    vector<float> angles(cameras);
    for (size_t i = 0; i < angles.size(); i++)
    {
        angles[i] = dtheta * static_cast<float>(i);
    }



    distribution rho("experiment_data/another_layer/no_object(196,196,1)layer35.cfd");
    const auto intersect = [](float _x, float _y, float _r)
    {
        return _x * _x + _y * _y <= _r * _r;
    };

    for (auto& r : rho.quantities())
    {
        r = (r - rho0) / rho0;
    }

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
    //vector<distribution> projected_points =  projector.calculate_projected_points(rho,angles);
    ilab::projection irho = projector.project(rho, angles);

    // Save the result at the PRJ format file.
    irho.save("p_data/layer_35/no_object(196,196,1)_"+to_string(dtheta)+"deg.prj");
}

void save_new_dist_another_layer(distribution dist){
    size_t layer_number = 34;
    ilab::distribution new_dist(196,196,1);
    for(int i=0;i<new_dist.width();i++){
        for(int j=0;j<new_dist.height();j++){
            new_dist.quantity(static_cast<size_t>(i),static_cast<size_t>(j),0)=dist.quantity(static_cast<size_t>(i),static_cast<size_t>(j),layer_number);
            new_dist.quantity(static_cast<size_t>(i),static_cast<size_t>(j),0) = new_dist.quantity(static_cast<size_t>(i),static_cast<size_t>(j),0);
        }
    }
    new_dist.save("new_distribution.cfd");

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

    save_newp_data(projection);
    //distribution baseDist("experiment_data/no_object(197,197,70).cfd");
    //save_new_dist_another_layer(baseDist);
    EMO emo(projection);//call EMO
    emo.evolution();
    //GACT ga(projection,"");
    //ga.Evolution();


	std::cout << "Complete" << std::endl;
	return 0;
}