//
// Created by Furukawa on 2017/05/19.
//

#include <iostream>
#include <cvaux.h>
#include "InverseDomain.h"

cv::Mat normalize(cv::Mat in) {

    //最大と最小の初期値を設定
    float max = 0;
    float min = FLT_MAX;

    //入力画像をclone
    cv::Mat out = in.clone();

    //最大と最小を取得
    for(int j = 0; j < in.rows; j++){
        for(int i = 0; i < in.cols; i++){
            if(max < in.at<float>(j,i)){
                max = in.at<float>(j,i);
            }
            if(min > in.at<float>(j,i)){
                min = in.at<float>(j,i);
            }
        }
    }

    //取得した最大値と最小値で正規化
    for(int j = 0; j < in.rows; j++){
        for(int i = 0; i < in.cols; i++){
            out.at<float>(j,i) = in.at<float>(j,i) * (max - min)/ 255;
        }
    }
    //0~255になったMatを返す
    return out;
}

//投影データから一次元フーリエ変換を用いて逆空間画像を構築する
InverseDomain::InverseDomain(ilab::projection p) {
    resize(p.height(),p.height());
    const auto intersect = [](float _x, float _y, float _r)
    {
        return _x * _x + _y * _y <= _r * _r;
    };
    const float center_x = static_cast<float>(this->width()) / 2.0f;
    const float center_y = static_cast<float>(this->height()) / 2.0f;
    for (size_t y = 0; y < this->height(); y++) {
        for (size_t x = 0; x < this->width(); x++){
            this->infomation(x, y) = info_type::unknown;
            if (intersect(static_cast<float>(x) - center_x, static_cast<float>(y) - center_y, center_y)){
                this->identity(x, y) = ilab::blank_type::quantity;
            }else{
                this->identity(x, y) = ilab::blank_type::outside;
            }
        }
    }
    FFT<float> FFT;
    vector<complex<float>> freq;
    for(int i=0;i<p.counts();i++){
        vector<float> d(p.height());
        for (int j = 0; j < p.height(); ++j) {
            d[j]=p.quantity(0, static_cast<size_t >(j), static_cast<size_t>(i));
        }
        //dはpのi番目の投影データが格納された１次元配列
        FFT.fwd(freq,d);//fftする

        for(int f=0;f<freq.size();f++){
            float r;
            float theta = p.angle(static_cast<size_t>(i));
            if (f >= freq.size() / 2) {
                r = freq.size() / 2 - (f - freq.size() / 2);
                int x = abs(static_cast<int>(r*sin(-theta-M_PI)+(freq.size()/2+1)));
                int y = abs(static_cast<int>(r*cos(-theta-M_PI)+(freq.size()/2)));
                this->quantity(static_cast<size_t>(x), static_cast<size_t>(y)) = freq[f];
                this->infomation(static_cast<size_t>(x), static_cast<size_t>(y)) = info_type::known;
            } else {
                r = f;
                int x = static_cast<int>(r*sin(-theta)+(freq.size()/2));
                int y = static_cast<int>(r*cos(-theta)+(freq.size()/2));
                this->quantity(static_cast<size_t>(x), static_cast<size_t>(y)) = freq[f];
                this->infomation(static_cast<size_t>(x), static_cast<size_t>(y)) = info_type::known;
            }
        }
    }
}
//実空間(再構成画像)から２次元フーリエ変換を用いて逆空間画像を構築する．
InverseDomain::InverseDomain(ilab::distribution dist) {
    resize(dist.width(),dist.height());
    this->identities() = dist.identities();

    InverseDomain temp(dist.width(),dist.height());
    FFT<float> FFT;
    vector<complex<float>> freq;

    for(int i=0;i<this->height();i++){
        vector<float> d(this->width());
        for(int j=0;j<this->width();j++){
            d[j]=dist.quantity(static_cast<size_t >(i), static_cast<size_t>(j),0);
        }
        FFT.fwd(freq,d);
        for(int j=0;j<this->width();j++){
            temp.quantity(static_cast<size_t>(i), static_cast<size_t>(j))=freq[j];
        }
    }
    for(int i=0;i<this->width();i++){
        vector<complex<float>> d(this->height());
        for(int j=0;j<this->height();j++){
            d[j]=temp.quantity(static_cast<size_t>(j), static_cast<size_t>(i));
        }
        FFT.fwd(freq,d);
        for(int j=0;j<this->height();j++){
            this->quantity(static_cast<size_t>(j), static_cast<size_t>(i))=freq[j];
        }
    }
    MatrixXcf mat = this->quantity2matrix();
    MatrixXcf tempmat(mat.rows(), mat.cols());
    int cx = (int)mat.cols() / 2;
    int cy = (int)mat.rows() / 2;// 追記:rowとcol逆
    tempmat.block(0, 0, cx, cy) = mat.block(cx, cy, cx, cy);
    tempmat.block(cx, cy, cx, cy) = mat.block(0, 0, cx, cy);
    tempmat.block(cx, 0, cx, cy) = mat.block(0, cy, cx, cy);
    tempmat.block(0, cy, cx, cy) = mat.block(cx, 0, cx, cy);
    mat = tempmat;
    this->matrix2quantity(mat);
    //Infomationは設定されないので注意

}

ilab::distribution InverseDomain::TwoDimFFT() {//自身の逆空間画像をフーリエ逆変換して再構成画像に直す．
    MatrixXcf mat = this->quantity2matrix();
    shift(mat);
    InverseDomain baseinv;
    baseinv.resize(this->width(),this->height());
    for(int i=0;i<this->height();i++){
        for(int j=0;j<this->width();j++){
            baseinv.quantity(static_cast<size_t>(j), static_cast<size_t>(i)) = mat(i,j);
        }
    }
    distribution newdist(this->width(),this->height(),1);
    newdist.identities() = this->identities();
    //フーリエ逆変換するフーリエ逆変換する
    distribution temp(this->width(),this->height(),1);
    FFT<float> FFT;
    vector<float> time;
    for(int i=0;i<this->height();i++){
        vector<complex<float>> d(this->width());
        for(int j=0;j<this->width();j++){
            d[j]=baseinv.quantity(static_cast<size_t >(i), static_cast<size_t>(j));
        }
        FFT.inv(time,d);
        for(int j=0;j<this->width();j++){
            temp.quantity(static_cast<size_t>(i), static_cast<size_t>(j),0)=time[j];
        }
    }
    for(int i=0;i<this->width();i++){
        vector<complex<float>> d(this->height());
        for(int j=0;j<this->height();j++){
            d[j]=temp.quantity(static_cast<size_t >(j), static_cast<size_t>(i),0);
        }
        FFT.inv(time,d);
        for(int j=0;j<this->height();j++){
            newdist.quantity(static_cast<size_t>(j), static_cast<size_t>(i),0)=time[j];
        }
    }
    return newdist;
}

bool InverseDomain::save(const std::string& _path) const{
    MatrixXf image(this->height(),this->width());
    for(int i=0;i<this->height();i++){
        for(int j=0;j<this->width();j++){
            image(i,j) = abs(this->quantity(static_cast<size_t>(i), static_cast<size_t>(j)));
        }
    }
    shift(image);
    cv::Mat resultimage;
    cv::eigen2cv(image, resultimage);
    resultimage = normalize(resultimage);
    cv::imwrite(_path,resultimage);
    return 1;
}

bool InverseDomain::save_notshift(const std::string& _path) const{
    MatrixXf image(this->height(),this->width());
    for(int i=0;i<this->height();i++){
        for(int j=0;j<this->width();j++){
            image(i,j) = abs(this->quantity(static_cast<size_t>(i), static_cast<size_t>(j)));
        }
    }
    cv::Mat resultimage;
    cv::eigen2cv(image, resultimage);
    resultimage = normalize(resultimage);
    return cv::imwrite(_path,resultimage);
}

void InverseDomain::visualize(){
    MatrixXf image(this->height(),this->width());
    for(int i=0;i<this->height();i++){
        for(int j=0;j<this->width();j++){
            image(i,j) = abs(this->quantity(static_cast<size_t>(i), static_cast<size_t>(j)));
        }
    }
    //shift(image);
    cv::Mat resultimage;
    cv::eigen2cv(image, resultimage);
    resultimage = normalize(resultimage);
    cv::imshow("Result", resultimage);
    cv::waitKey(0);
}

void const InverseDomain::shift(MatrixXf &freq) const{
    MatrixXf temp(freq.rows(), freq.cols());

    int cx = (int)freq.cols() / 2;
    int cy = (int)freq.rows() / 2;// 追記:rowとcol逆
    temp.block(0, 0, cx, cy) = freq.block(cx, cy, cx, cy);
    temp.block(cx, cy, cx, cy) = freq.block(0, 0, cx, cy);
    temp.block(cx, 0, cx, cy) = freq.block(0, cy, cx, cy);
    temp.block(0, cy, cx, cy) = freq.block(cx, 0, cx, cy);

    freq = temp;
}
void const InverseDomain::shift(MatrixXcf &freq) const{
    MatrixXcf temp(freq.rows(), freq.cols());

    int cx = (int)freq.cols() / 2;
    int cy = (int)freq.rows() / 2;// 追記:rowとcol逆
    temp.block(0, 0, cx, cy) = freq.block(cx, cy, cx, cy);
    temp.block(cx, cy, cx, cy) = freq.block(0, 0, cx, cy);
    temp.block(cx, 0, cx, cy) = freq.block(0, cy, cx, cy);
    temp.block(0, cy, cx, cy) = freq.block(cx, 0, cx, cy);
//    cout<<temp(temp.rows()-1,temp.cols()-1)<<endl;

    freq = temp;
}

MatrixXcf InverseDomain::quantity2matrix() {
    MatrixXcf image(this->height(),this->width());
    for(int i=0;i<this->height();i++){
        for(int j=0;j<this->width();j++){
            image(i,j) = this->quantity(static_cast<size_t>(j), static_cast<size_t>(i));
        }
    }
    return image;
}

void InverseDomain::matrix2quantity(MatrixXcf mat) {
    for(int i=0;i<this->height();i++){
        for(int j=0;j<this->width();j++){
            this->quantity(static_cast<size_t>(j), static_cast<size_t>(i)) = mat(i,j);
        }
    }
}





