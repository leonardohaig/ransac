#include "LineParamEstimator.h"
#include "PlaneParamEstimator.h"
#include "Ransac.h"
#include<iostream>
using namespace std;
/*
 * Example of using the Ransac class for robust parameter estimation.
 *
 * Author: Ziv Yaniv (zivy@cs.huji.ac.il)
 */
/*
int  main(int argc, char* argv[])
{
	std::vector<double> lineParameters;
	LineParamEstimator lpEstimator(0.5); //for a point to be on the line it has to be closer than 0.5 units from the line
	std::vector<Point2D> pointData;
	std::vector<Point2D *> pointDataPtr;
	int numForEstimate = 2;
	int numSamples = 20;
	int numOutliers = 80;
	double desiredProbabilityForNoOutliers = 0.999;
	double maximalOutlierPercentage = 0.1 + (double)numOutliers/(double)(numSamples + numOutliers);
	double noiseSpreadRadius = 0.4;
	double outlierSpreadRadius = 10;
	int i;
	double newX, newY, dx, dy, norm;

    //1.Create data with outliers

    //randomly select a direction [dx,dy] and create a line passing through the origin
	//for each point sampled on the line add random noise, finally add outlying 
	//points in the direction of the line normal.

	srand((unsigned)time(NULL)); //seed random number generator
   
	//get random direction
	dx = rand();
	dy = rand();
//	dx = rand()%100;
//	dy = rand()%100;
	norm = sqrt(dx*dx + dy*dy); 
	dx/= norm;
	dy/= norm;
	dx *= (rand() > RAND_MAX/2 ? 1 : -1);
	dy *= (rand() > RAND_MAX/2 ? 1 : -1);
//	dx *= ((rand()%100) > RAND_MAX/2 ? 1 : -1);
//	dy *= ((rand()%100) > RAND_MAX/2 ? 1 : -1);

    //add 'numSamples' points
	for(i=0; i<numSamples; i++)
	{
		newX = i*dx + noiseSpreadRadius*(double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1);
		newY = i*dy + noiseSpreadRadius*(double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1);
//		newX = i*dx + noiseSpreadRadius*(double)(rand()%100)/(double)RAND_MAX * ((rand()%100) > RAND_MAX/2 ? 1 : -1);
//		newY = i*dy + noiseSpreadRadius*(double)(rand()%100)/(double)RAND_MAX * ((rand()%100) > RAND_MAX/2 ? 1 : -1);
		
		pointDataPtr.push_back(new Point2D(newX,newY));
		pointData.push_back(*(pointDataPtr[i]));
	}

	//'numOutliers' points
	double centerX = -dy*100;
	double centerY = dx*100;
	for(i=0; i<numOutliers; i++) 
	{
		newX = centerX + outlierSpreadRadius * (double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1);
		newY = centerY + outlierSpreadRadius * (double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1);
//		newX = centerX + outlierSpreadRadius * (double)(rand()%100)/(double)RAND_MAX * ((rand()%100) > RAND_MAX/2 ? 1 : -1);
//		newY = centerY + outlierSpreadRadius * (double)(rand()%100)/(double)RAND_MAX * ((rand()%100) > RAND_MAX/2 ? 1 : -1);
		pointDataPtr.push_back(new Point2D(newX,newY));
		pointData.push_back(*(pointDataPtr[pointDataPtr.size()-1]));
	}
	
	double dotProd;
                        
	//2. Compare least squares approach to Ransac

	cout<<"Total number of points used: "<<pointData.size()<<endl;
	cout<<"Number of outliers: "<<numOutliers<<endl;
	            //The real line parameters
	cout<<"Real line parameters [nx,ny,ax,ay]\n\t [ "<<-dy<<", "<<dx<<", 0, 0 ]"<<endl;
//	cout<<"Real line parameters [-ny,nx,ax,ay]\n\t [ "<<-dy<<", "<<dx<<", 0, 0 ]"<<endl;

    //A least squares estimate of the line parameters
	lpEstimator.leastSquaresEstimate(pointDataPtr,lineParameters);
	cout<<"Least squares line parameters: [n_x,n_y,a_x,a_y]\n\t [ "<<lineParameters[0]<<", "<<lineParameters[1]<<", ";
	cout<<lineParameters[2]<<", "<<lineParameters[3]<<" ]"<<endl;
	dotProd = lineParameters[0]*(-dy) + lineParameters[1]*dx;
	cout<<"\tDot product of real and computed line normals[+-1=correct]: "<<dotProd<<endl;
    dotProd = (-dy)*(lineParameters[2]) + dx*lineParameters[3];
	cout<<"\tCheck if computed point is on real line [0=correct]: "<<dotProd<<endl;

    //A RANSAC estimate of the line parameters

//	double usedData = Ransac<Point2D,double>::compute(lineParameters, &lpEstimator,  pointData, numForEstimate); //for line, numForEstimate == 2
	                                                                                                             //for plane, numForEstimate == 3

    double usedData = Ransac<Point2D,double>::compute(lineParameters, &lpEstimator,  pointData, numForEstimate, desiredProbabilityForNoOutliers, maximalOutlierPercentage); 

	cout<<"RANSAC line parameters [n_x,n_y,a_x,a_y]\n\t [ "<<lineParameters[0]<<", "<<lineParameters[1]<<", ";
	cout<<lineParameters[2]<<", "<<lineParameters[3]<<" ]"<<endl;

	dotProd = lineParameters[0]*(-dy) + lineParameters[1]*dx;
	cout<<"\tDot product of real and computed line normals[+-1=correct]: "<<dotProd<<endl;

//	dotProd = lineParameters[0]*dx - lineParameters[1]*dy;
	dotProd = lineParameters[0]*dx + lineParameters[1]*dy;
	cout<<"\tDot product of real and computed line direction[0=correct]: "<<dotProd<<endl;

    dotProd = (-dy)*(lineParameters[2]) + dx*lineParameters[3];
	cout<<"\tCheck if computed point is on real line [0=correct]: "<<dotProd<<endl;

	dotProd = lineParameters[0]*(lineParameters[2]) + lineParameters[1]*lineParameters[3];
	cout<<"\talso can be Checked if computed point is on real line [0=correct]: "<<dotProd<<endl;

	cout<<"\tPercentage of points which were used for final estimate: "<<usedData<<endl;
	
	return 0;
}  */

int main(int argc, char* argv[])
{
    double array[130][3];
    double newX, newY, newZ;

    for (int i = 0; i < 12; i++) 
	{
		for (int j = 0; j < 3; j++)
		{
			array[i][j] = 0.0;
		}
    }
	
	srand((unsigned)time(NULL)); //seed random number generator
	double noiseSpreadRadius = 0.4;
	
    for (int i = 0; i < 12; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            array[i][j] = (double)rand();
        }
    }

	/*---test estimate() function begin---*/
	double arrReal[3][3];

	arrReal[0][0] = 1.0;
    arrReal[0][1] = 0.0;
	arrReal[0][2] = 0.0;

	arrReal[1][0] = 0.0;
    arrReal[1][1] = 1.0;
	arrReal[1][2] = 0.0;

	arrReal[2][0] = 0.0;
    arrReal[2][1] = 0.0;
	arrReal[2][2] = 1.0;

   	std::vector<double> realParameters;
	std::vector<cv::Point3d> pointReal;
	std::vector<cv::Point3d *> pointRealPtr;

	for (int i = 0; i < 3; i++) 
	{
		newX = arrReal[i][0];
		newY = arrReal[i][1];
		newZ = arrReal[i][2];
		pointRealPtr.push_back(new cv::Point3d(newX, newY, newZ));
		pointReal.push_back(*(pointRealPtr[pointRealPtr.size()-1]));
	}
    /*---test estimate() function end---*/

	/*---test leastSquaresEstimate() function begin--- */
    for (int i = 0; i < 40; i++)  //point 0~11, are on the plane
    {
        array[i][0] = 1.0 + noiseSpreadRadius*(double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1);
		array[i][1] = 0.0 + noiseSpreadRadius*(double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1);
		array[i][2] = 0.0 + noiseSpreadRadius*(double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1);
    }

	for (int i = 40; i < 80; i++)
    {
        array[i][0] = 0.0 + noiseSpreadRadius*(double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1);
		array[i][1] = 1.0 + noiseSpreadRadius*(double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1);
		array[i][2] = 0.0 + noiseSpreadRadius*(double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1);
    }

	for (int i = 80; i < 120; i++)
    {
        array[i][0] = 0.0 + noiseSpreadRadius*(double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1);
		array[i][1] = 0.0 + noiseSpreadRadius*(double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1);
		array[i][2] = 1.0 + noiseSpreadRadius*(double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1);
    }

	for (int i = 120; i < 130; i++)
    {
        array[i][0] = 100.0 + noiseSpreadRadius*(double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1);
		array[i][1] = 100.0 + noiseSpreadRadius*(double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1);
		array[i][2] = 100.0 + noiseSpreadRadius*(double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1);
    }

    std::vector<double> planeParameters;
	std::vector<cv::Point3d> pointData;
	std::vector<cv::Point3d *> pointDataPtr;

    /*---test leastSquaresEstimate() function end--- */
    /*
	PlaneParamEstimator lpEstimator(0.5); 
	lpEstimator.estimate(pointRealPtr, realParameters);
	lpEstimator.leastSquaresEstimate(pointDataPtr, planeParameters);

	double realA = realParameters[0];  //point-normal formZ: A(x-x0)+B(y-y0)+C(y-y0)=0
	double realB = realParameters[1];
	double realC = realParameters[2];
	double realX = realParameters[3];
	double realY = realParameters[4];
	double realZ = realParameters[5];

	double resultA = planeParameters[0];  //point-normal form
	double resultB = planeParameters[1];
	double resultC = planeParameters[2];
	double resultX = planeParameters[3]; 
	double resultY = planeParameters[4]; 
    double resultZ = planeParameters[5]; 
    */

	/*---test RANSAC begin--- */
    //point 12~15, are not on one plane
//  array[12][0] = 1.1;
//	array[12][1] = 0.3;
//	array[12][2] = 0.2;

//  array[13][0] = 0.1;
//	array[13][1] = 1.3;
//	array[13][2] = -0.2;

//  array[14][0] = 0.1;
//	array[14][1] = 0.23;
//	array[14][2] = 0.8;

//	array[15][0] = 0.1;
//	array[15][1] = 300.86;
//	array[15][2] = 0.04;

	for (int i = 0; i < 130; i++)
	{
		newX = array[i][0];
		newY = array[i][1];
		newZ = array[i][2];
		pointDataPtr.push_back(new cv::Point3d(newX, newY, newZ));
		pointData.push_back(*(pointDataPtr[pointDataPtr.size()-1]));
	}

	PlaneParamEstimator lpEstimator(0.5); 
	int numForEstimate = 3;

	double usedData = Ransac<cv::Point3d,double>::compute(planeParameters, &lpEstimator,  pointData, numForEstimate);

	double resultA = planeParameters[0];
	double resultB = planeParameters[1];
	double resultC = planeParameters[2];
	double resultX = planeParameters[3]; 
	double resultY = planeParameters[4]; 
    double resultZ = planeParameters[5]; 

	cout<<"RANSAC plane parameters [A, B, C, X0, Y0, Z0]\n\t [ "<<planeParameters[0]<<", "<<planeParameters[1]<<", "
	                                                        	<<planeParameters[2]<<", "<<planeParameters[3]<<", "
		                                                        <<planeParameters[4]<<", "<<planeParameters[5]<<", ";
	cout<<"\tPercentage of points which were used for final estimate: "<<usedData<<endl;
    /*---test RANSAC end--- */

	return 0;
}
