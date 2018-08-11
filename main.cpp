#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

void cal(float panelParam[4])
{
    panelParam[0] = 1;
}

/*
int main()
{


    //vector<Point2f>===> vector<Point2i>
    vector<Point2f> v1;
    vector<Point2i> v2;
    Point2f p;

    p.x=1.1;
    v1.push_back(Point2f(1.1,2.4));
    v1.push_back(Point2f(2.1,3.5));
    v1.push_back(Point2f(-3.4,-4.6));
    Mat(v1).convertTo(v2,CV_32SC1);


    float K[]={0.0,1.0,2.0,3.0,
                0.0,1.0,2.0,3.0,
                0.0,1.0,2.0,3.0};


    cv::Mat origMat(3,4,CV_32FC1,K);

    cout<<origMat<<endl;



    vector<Point3f> p1,p2;
    p1 = Mat_<Point3f>(origMat);//说明：取出Mat的每一行，然后作为一个向量


    p2 = Mat_<Point3f>(origMat.t());

    cout<<origMat.t()<<endl;

    Mat senM = Mat(p2);
    senM.convertTo(senM,CV_32FC1);
    //cvtColor(senM,senM,CV_32FC1);
    cout<<senM<<endl;




    return 0;
}

 */

