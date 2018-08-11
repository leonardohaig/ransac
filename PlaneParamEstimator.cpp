#include <math.h>
#include "PlaneParamEstimator.h"
#include "Matrix.h"

PlaneParamEstimator::PlaneParamEstimator(double delta) : m_deltaSquared(delta*delta) {}
/*****************************************************************************/
/*
 * Compute the plane parameters  [n_x, n_y, n_z, a_x,a_y, a_z]
 */
void PlaneParamEstimator::estimate(std::vector<cv::Point3d *> &data, std::vector<double> &parameters)
{
	parameters.clear();
	if(data.size() < 3)
		return;

	double AB_x = data[1]->x - data[0]->x;  
	double AB_y = data[1]->y - data[0]->y;
	double AB_z = data[1]->z - data[0]->z;

	double AC_x = data[2]->x - data[0]->x;
	double AC_y = data[2]->y - data[0]->y;
	double AC_z = data[2]->z - data[0]->z;

	double A = AB_y*AC_z - AB_z*AC_y;
	double B = AB_z*AC_x - AB_x*AC_z;
	double C = AB_x*AC_y - AB_y*AC_x;

	double norm = sqrt(A*A + B*B + C*C);
	
	parameters.push_back(A/norm);
	parameters.push_back(B/norm);
	parameters.push_back(C/norm);
	parameters.push_back(data[0]->x);
	parameters.push_back(data[0]->y);	
	parameters.push_back(data[0]->z);
}
/*****************************************************************************/
/*
 * Compute the plane parameters  [n_x, n_y, n_z, a_x, a_y, a_z]
 */
void PlaneParamEstimator::leastSquaresEstimate(std::vector<cv::Point3d *> &data, std::vector<double> &parameters)
{
	double *Matrix[3],*IMatrix[3], Y[3];
	double A, B, C;
	A = B = C = 0.0;

    for (int i = 0;i < 3;i++)
    {
        Matrix[i]  = new double[3];
        IMatrix[i] = new double[3];
    }

    for (int i = 0; i < 3;i++)
    {
        for (int j = 0;j < 3;j++)
        {
            *(Matrix[i] + j) = 0.0;
        }
    }

	for (int i = 0; i < 3; i++)
	{
		Y[i] = 0.0;
	}

    for (int i = 0; i < data.size(); i++)
    {
        *(Matrix[0]) += data[i]->x*data[i]->x;
        *(Matrix[1]) += data[i]->y*data[i]->x;
        *(Matrix[2]) += data[i]->z*data[i]->x;
        Y[0] -= data[i]->x;
    }
	for (int i = 0; i < data.size(); i++)
    {
        *(Matrix[0] + 1) += data[i]->x*data[i]->y;
        *(Matrix[1] + 1) += data[i]->y*data[i]->y;
        *(Matrix[2] + 1) += data[i]->z*data[i]->y;
        Y[1] -= data[i]->y;
    }

    for (int i = 0; i < data.size(); i++)
    {
        *(Matrix[0] + 2) += data[i]->x*data[i]->z;
        *(Matrix[1] + 2) += data[i]->y*data[i]->z;
        *(Matrix[2] + 2) += data[i]->z*data[i]->z;
        Y[2] -= data[i]->z;
    }

	double d = Determinant(Matrix, 3);

	if (abs(d) < 0.0001)
	{
		printf("\n ¾ØÕóÆæÒì");
//		getchar();

		return;
	}

	Inverse(Matrix, IMatrix, 3, d);
	for (int i = 0; i < 3; i++)
	{
		A += *(IMatrix[0] + i)*Y[i];
		B += *(IMatrix[1] + i)*Y[i];
		C += *(IMatrix[2] + i)*Y[i];
	}

	double norm = sqrt(A*A + B*B + C*C);

	double meanX, meanY, meanZ;
	meanX = meanY = meanZ = 0.0;
	for (int i = 0; i < data.size(); i++)
	{
		meanX += data[i]->x;
	    meanY += data[i]->y;
        meanZ += data[i]->z;
	}

	meanX /= data.size();
	meanY /= data.size();
	meanZ /= data.size();

	parameters.push_back(A/norm);
	parameters.push_back(B/norm);
	parameters.push_back(C/norm);
	parameters.push_back(meanX);
	parameters.push_back(meanY);
	parameters.push_back(meanZ);

//	printf("\n A = %5.3f, B = %5.3f, C = %5.3f", A, B, C);

	for (int i = 0; i < 3; i++)
	{
		delete [] Matrix[i];
		delete [] IMatrix[i];
	}
}
/*****************************************************************************/
/*
 * Given the plane parameters  [n_x,n_y,n_z,a_x,a_y,a_z] check if
 * [n_x, n_y, n_z] dot [data.x-a_x, data.y-a_y, data.z-a_z] < m_delta
 */
bool PlaneParamEstimator::agree(std::vector<double> &parameters, cv::Point3d &data)
{
	double signedDistance = parameters[0]*(data.x - parameters[3])
		                   +parameters[1]*(data.y - parameters[4])
						   +parameters[2]*(data.z - parameters[5]);

	return ((signedDistance*signedDistance) < m_deltaSquared);
}
/*****************************************************************************/
void PlaneParamEstimator::debugTest(ostream &out)
{

}
