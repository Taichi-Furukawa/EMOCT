//
// Created by Furukawa on 2017/05/19.
//

#ifndef EMOCT_INVERSEDOMAIN_H
#define EMOCT_INVERSEDOMAIN_H
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <complex.h>
#include <Eigen/Core>
#include <unsupported/Eigen/FFT>
#include <numeric>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/eigen.hpp>
#include "CUDAMath.h"
#include "ilab.h"

using namespace std;
using namespace ilab;
using namespace Eigen;

enum class info_type : int32_t
{
    unknown = 0,
    known = 1
};

class InverseDomain {
public:
    //initializer
    InverseDomain(projection p);
    InverseDomain(distribution dist);

    InverseDomain()
            : m_width(0), m_height(0)
    {
    }

    InverseDomain(size_t _width, size_t _height)
    : m_width(_width), m_height(_height), m_quantities(_width * _height), m_identities(_width * _height, blank_type::quantity)
    {
    }

    //functions
    distribution TwoDimFFT();

    size_t width() const
    {
        return m_width;
    }

    size_t height() const
    {
        return m_height;
    }

    std::vector<complex<float>>& quantities()
    {
        return m_quantities;
    }

    const std::vector<complex<float>>& quantities() const
    {
        return m_quantities;
    }

    std::vector<blank_type>& identities()
    {
        return m_identities;
    }

    const std::vector<blank_type>& identities() const
    {
        return m_identities;
    }

    std::vector<info_type>& informations()
    {
        return m_infomation;
    }


    void resize(size_t _width, size_t _height)
    {
        m_width = _width;
        m_height = _height;
        m_quantities.resize(_width * _height);
        m_identities.resize(_width * _height);
        m_infomation.resize(_width * _height);
    }

    bool save(const std::string& _path) const;
    bool save_notshift(const std::string& _path) const;
    //properties
    complex<float>& quantity(size_t _x, size_t _y){
        return m_quantities.at(_y * m_width + _x);
    }

    const complex<float>& quantity(size_t _x, size_t _y) const {
        return m_quantities.at(_y * m_width + _x);
    }

    blank_type& identity(size_t _x, size_t _y) {
        return m_identities.at(_y * m_width + _x);
    }

    const blank_type& identity(size_t _x, size_t _y) const {
        return m_identities.at(m_width * m_height + _y * m_width + _x);
    }

    info_type& infomation(size_t _x, size_t _y) {
        return m_infomation.at(_y * m_width + _x);
    }
    void visualize();

private:
    size_t m_width;
    size_t m_height;
    std::vector<complex<float>> m_quantities;
    std::vector<blank_type> m_identities;
    std::vector<info_type> m_infomation;

};


#endif //EMOCT_INVERSEDOMAIN_H
