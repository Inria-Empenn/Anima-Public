#pragma once

#include "animaTODEstimatorImageFilter.h"

#include <vtkVector.h>
#include <vtkCellArray.h>
#include <vtkCellArrayIterator.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataReader.h>

#include <animaODFSphericalHarmonicBasis.h>
#include <animaODFFunctions.h>
#include <animaVectorOperations.h>
#include <animaShapesReader.h>
#include <animaReadWriteFunctions.h>
#include <animaMatrixOperations.h>
#include <animaSphereOperations.h>
#include "animaLogExpMapsUnitSphere.h"

#include <random>

namespace anima
{

void TODEstimatorImageFilter::BeforeThreadedGenerateData()
{
    Superclass::BeforeThreadedGenerateData();

    TRefImage::Pointer refImg = anima::readImage<TRefImage>(m_RefFileName);

    unsigned int vectorLength = (m_LOrder + 1)*(m_LOrder + 2)/2;
    m_VectorLength = vectorLength;
    TOutputImage *output = this->GetOutput();
    output->SetVectorLength(vectorLength);
    output->SetRegions(refImg->GetLargestPossibleRegion());
    output->SetSpacing(refImg->GetSpacing());
    output->SetDirection(refImg->GetDirection());
    output->SetOrigin(refImg->GetOrigin());
    output->Allocate();

    VectorType tmp;
    tmp.SetSize((m_LOrder+1)*(m_LOrder+2)/2);
    tmp.Fill(0);
    output->FillBuffer(tmp);

    test = itk::Image<int, 3>::New();
    test->SetRegions(refImg->GetLargestPossibleRegion());
    test->SetSpacing(refImg->GetSpacing());
    test->SetDirection(refImg->GetDirection());
    test->SetOrigin(refImg->GetOrigin());
    test->Allocate();

    m_CstDir[0] = 0;
    m_CstDir[1] = 0;
    m_CstDir[2] = 1;

    m_NbSample = 200;
    anima::GetSphereEvenSampling(m_SphereSampl, m_NbSample);

    m_SpherHarm.set_size(m_NbSample, vectorLength);

    m_ODFSHBasis = new anima::ODFSphericalHarmonicBasis(m_LOrder);
    this->precomputeSH();

    const int Xmax = output->GetLargestPossibleRegion().GetSize()[0];
    const int Ymax = output->GetLargestPossibleRegion().GetSize()[1];
    const int Zmax = output->GetLargestPossibleRegion().GetSize()[2];

    const int Xmin = output->GetLargestPossibleRegion().GetIndex()[0];
    const int Ymin = output->GetLargestPossibleRegion().GetIndex()[1];
    const int Zmin = output->GetLargestPossibleRegion().GetIndex()[2];

    const int dirX = output->GetDirection()[0][0];
    const int dirY = output->GetDirection()[1][1];
    const int dirZ = output->GetDirection()[2][2];

    anima::ShapesReader trackReader;
    trackReader.SetFileName(m_InputFileName);
    trackReader.Update();

    this->ComputeCoefs();

    m_ImgDir.resize(Xmax * Ymax * Zmax);

    vtkSmartPointer<vtkPolyData> tracks = trackReader.GetOutput();
    vtkSmartPointer<vtkCellArray> cells = tracks->GetLines();
    vtkSmartPointer<vtkPoints> points = tracks->GetPoints();

    const vtkIdType *indices;
    vtkIdType numberOfPoints;

    float nbLines = tracks->GetNumberOfLines();
    float lineCount = 0;

    std::cout << "Number of fibers + test : " << nbLines << std::endl;

    anima::ODFSphericalHarmonicBasis basis(m_LOrder);

    auto iter = vtk::TakeSmartPointer(cells->NewIterator());
    for (iter->GoToFirstCell(); !iter->IsDoneWithTraversal(); iter->GoToNextCell())
    {
        iter->GetCurrentCell(numberOfPoints, indices);
        //        std::cout  << "Progress : " << (int) (((float)lineCount / nbLines)*100) <<  "%\r";
        //        fflush (stdout);
        lineCount++;

        FiberType fiber = readFiber(numberOfPoints, indices, points);
        int nbPoints = fiber.size();


        for(int i = 0; i < nbPoints; i++)
        {
            DirType dir = getFiberDirection(i, fiber);
            //        OutputImagePointer::DirectionType refDir = output->GetDirection();

            //            if(output->GetDirection()[0][0] > 0)
            //                dir[0] *= -1;
            //            if(output->GetDirection()[1][1] < 0)
            //                dir[1] *= -1;
            //            if(output->GetDirection()[2][2] < 0)
            //                dir[2] *= -1;
            //            dir[0] *= -1;
            //                        dir[1] *= -1;
            //            dir[2] *= -1;

            //            anima::Normalize(dir, dir);
            PointType point = getCenterVoxel(i, fiber);
            itk::Index<3> index, index2;

            output->TransformPhysicalPointToIndex(point, index);

            index[0] *= dirX;
            index[1] *= dirY;
            index[2] *= dirZ;

            index2[0] = std::floor(point[0]);
            index2[1] = std::floor(point[0]);
            index2[2] = std::floor(point[0]);

            m_ImgDir[index[0] + Xmax * (index[1] + Ymax * index[2])].push_back(dir);

            //    //    itk::ImageRegionIterator <itk::Image<int, 3>> testItr(test, test->GetLargestPossibleRegion());
            //    //    while(!testItr.IsAtEnd())
            //    //    {
            //    //        itk::Image<int, 3>::IndexType index;
            //    //        index = testItr.GetIndex();
            //    //        testItr.Set(imgDir[index[0] + Xmax * (index[1] + Ymax * index[2])].size());
            //    //        ++testItr;
            //    //    }
            //    //    anima::writeImage<itk::Image<int, 3>>("./test.nii.gz", test);


        }
    }
}

void TODEstimatorImageFilter::ComputeCoefs()
{
    m_GaussCoefs.SetSize(91);
    m_GaussCoefs.Fill(0.0);

    double l1 = 1.71e-4;
    double l2 = 2e-5;
    double d = l1 - l2;
    double spi = sqrt(M_PI);
    double pi = M_PI;
    double at = M_PI+2*atan((d-l2)/(2*sqrt(d*l2)));

    double c0 = 1/(2*sqrt(pi));
    double c1 = - sqrt(5) * (3*sqrt(l2)*(d+l2)*(pi + atan((d-l2)/(2*sqrt(l2*d)))) - 4*sqrt(d)*(2*d+3*l2)) / (16*sqrt(pi)*pow(d,1.5));
    double c2 = 3 * (4*sqrt(d)*(16*pow(d,2)+115*d*l2+105*pow(l2,2)) - 15*sqrt(l2)*(3*pow(d,2)+10*d*l2+7*pow(l2,2))*(pi+2*atan((d-l2)/(2*sqrt(d*l2))))) / (128*spi*pow(d,2.5));
    double c3 = -105*sqrt(13)*(at*sqrt(l2)*(5*pow(d,3)+35*l2*pow(d,2)+63*pow(l2,2)*d+33*pow(l2,3)) - 4*sqrt(d)*(33*pow(l2,3)+52*d*pow(l2,2)+103*pow(d,2)*l2/5+128*pow(d,3)/105)) / (1024*spi*pow(d,3.5));
    double c4 = -315*sqrt(17) * (at*sqrt(l2)*(35*pow(d,4)+420*pow(d,3)*l2+1386*pow(d,2)*pow(l2,2)+1716*pow(l2,3)*d+715*pow(l2,4)) - 4*sqrt(d)*(715*pow(l2,4)+4433*pow(l2,3)*d/3+957*pow(l2,2)*pow(d,2)+6967*l2*pow(d,3)/35+2048*pow(d,4)/315)) / (16384 * spi * pow(d,4.5));
    double c5 = -24255*sqrt(21) * (at * sqrt(l2)*(9*pow(d,5)+165*pow(d,4)*l2+858*pow(d,3)*pow(l2,2)+12870*pow(d,2)*pow(l2,3)/7+12155*d*pow(l2,4)/7+4199*pow(l2,5)/7) - 4*sqrt(d)*(4199*pow(l2,5)/7+32266*pow(l2,4)*d/21+20691*pow(l2,3)*pow(d,2)/15+24830*pow(l2,2)*pow(d,3)/49+28799*l2*pow(d,4)/441+32768*pow(d,5)/24255)) / (262144*pow(d,5.5)*spi);


    //    m_GaussCoefs[0] = this->m_ODFSHBasis->getNthSHValueAtPosition(0, 0, 0, 0);
    //    std::vector<double> Dir(3), sphDir(2);
    //    Dir[0] = 0;
    //    Dir[1] = 0;
    //    Dir[2] = 1;
    //    anima::TransformCartesianToSphericalCoordinates(Dir, sphDir);
    //    int k = 0;
    //    for(int l = 0; l <= m_LOrder; l+=2)
    //    {
    //        for(int m = -l; m <= l; m++)
    //        {
    //            m_GaussCoefs[k] = this->m_ODFSHBasis->getNthSHValueAtPosition(l, m, sphDir[0], sphDir[1]);
    //            k++;
    //        }
    //    }

    m_GaussCoefs[0] = c0;
    m_GaussCoefs[3] = c1;
    m_GaussCoefs[10] = c2;
    m_GaussCoefs[21] = c3;
    m_GaussCoefs[36] = c4;
    m_GaussCoefs[55] = c5;
}

void TODEstimatorImageFilter::DynamicThreadedGenerateData(const TODEstimatorImageFilter::OutputImageRegionType &outputRegionForThread)
{
    TOutputImage *output = this->GetOutput();
    itk::ImageRegionIterator <TOutputImage> outItr(output, outputRegionForThread);
    typename TOutputImage::IndexType index, index_out;

    VectorType nullVector;
    nullVector.SetSize((m_LOrder+1)*(m_LOrder+2)/2);
    nullVector.Fill(0);

    const int Xmax = output->GetLargestPossibleRegion().GetSize()[0];
    const int Ymax = output->GetLargestPossibleRegion().GetSize()[1];
    const int Zmax = output->GetLargestPossibleRegion().GetSize()[2];
    const int dirX = output->GetDirection()[0][0];
    const int dirY = output->GetDirection()[1][1];
    const int dirZ = output->GetDirection()[2][2];


    while(!outItr.IsAtEnd())
    {


        index = outItr.GetIndex();

        //        index[0] *= -1;
        //        index[1] *= -1;
        //        index[2] *= -1;


        //        index[0] *= dirX;
        //        index[1] *= dirY;
        //        index[2] *= dirZ;

        std::vector<VectorType> vecCoefs;
        if(!(m_ImgDir[index[0] + Xmax * (index[1] + Ymax * index[2])].size()) == 0)
        {
            int tmpSize = m_ImgDir[index[0] + Xmax * (index[1] + Ymax * index[2])].size();
            DirVectorType mainDirs;
            DirVectorType localDir(tmpSize);
            localDir = m_ImgDir[index[0] + Xmax * (index[1] + Ymax * index[2])];
            this->getMainDirections(localDir, mainDirs);

            //            std::cout << "----------------------------------------" << std::endl << std::endl;
            //            for(int i = 0; i < tmpSize; i++)
            //                std::cout << localDir[i] << ", " << std::endl;
            //            std::cout << std::endl << std::endl;
            //            for(int i = 0; i < mainDirs.size(); i++)
            //                std::cout << mainDirs[i] << ", " << std::endl;

            //            std::cout << std::endl << std::endl;

            VectorType tmpCoefs(m_VectorLength);
            tmpCoefs.Fill(0);
            for(int i = 0; i < mainDirs.size(); i++)
            {
                if(!isnan(mainDirs[i][0]))
                {
                    this->getSHCoef(mainDirs[i], tmpCoefs);
                    //                    std::vector<double> Dir(3), sphDir(2);
                    //                    for(int k = 0; k < 3; k++)
                    //                        Dir[k] = mainDirs[i][k];
                    //                    anima::TransformCartesianToSphericalCoordinates(Dir, sphDir);
                    //                    int j = 0;
                    //                    for(int l = 0; l <= m_LOrder; l+=2)
                    //                    {
                    //                        for(int m = -l; m <= l; m++)
                    //                        {
                    //                            tmpCoefs[j] = this->m_ODFSHBasis->getNthSHValueAtPosition(l, m, sphDir[0], sphDir[1]);
                    //                            j++;
                    //                        }
                    //                    }

                    vecCoefs.push_back(tmpCoefs);
                }
            }

            VectorType resCoefs(m_VectorLength);
            resCoefs.Fill(0);
            if(vecCoefs.size() != 0)
            {
                //                this->averageODFs(vecCoefs, resCoefs);
                //                for(int k = 0; k < vecCoefs.size(); k++)
                //                {
                //                    for(int i = 0; i < m_VectorLength; i++)
                //                        resCoefs[i] += (vecCoefs[k][i] / vecCoefs.size());
                //                }
                for(int k = 0; k < vecCoefs.size(); k++)
                {
                    for(int i = 0; i < m_VectorLength; i++)
                        resCoefs[i] += vecCoefs[k][i];
                }
                if(m_Normalize)
                {
                    for(int i = 0; i < m_VectorLength; i++)
                        resCoefs[i] /= vecCoefs.size();
                }
                outItr.Set(resCoefs);
            }

            else
                outItr.Set(nullVector);

            std::vector<double> testSh(m_VectorLength);
            for(int i = 0; i < m_VectorLength; i++)
            {
                testSh[i] = resCoefs[i];
            }
            localDir.clear();
            localDir.shrink_to_fit();
        }
        else
            outItr.Set(nullVector);

        ++outItr;

        this->IncrementNumberOfProcessedPoints();
    }

}

void TODEstimatorImageFilter::AfterThreadedGenerateData()
{

}

//void TODEstimatorImageFilter::GenerateData()
//{






//}




void TODEstimatorImageFilter::getMainDirections(DirVectorType inDirs, DirVectorType &mainDirs)
{

    if(inDirs.size() == 0)
    {
        mainDirs.resize(0);
        return;
    }
    const int numIt = 200;
    int numDirs = inDirs.size(), crit = 0, numClusters = 0;

    if(numDirs < 4)
    {
        numClusters = numDirs;
        mainDirs.resize(numClusters);
        for(int i = 0; i < numClusters; i++)
        {
            mainDirs[i] = inDirs[i];
        }
        return;
    }
    else
        numClusters = 4;

    DirVectorType clustersAverage(numClusters);
    std::vector<int> clusters(numDirs);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, numDirs);


    for(int i = 0; i < numClusters; i++)
    {
        int k = distrib(gen);
        clustersAverage[i] = inDirs[k];
    }

    int It = 0;
    bool stop = false;
    double dist = 0;
    while(It < numIt && !stop)
    {
        std::vector<int> clustersOld(numDirs);
        for(int i = 0; i < numDirs; i++)
            clustersOld[i] = clusters[i];

        for(int i = 0; i < numDirs; i++)
        {
            double minDist = 1e10;
            for(int j = 0; j < numClusters; j++)
            {
                //                double dist = getEuclideanDistance(inDirs[i], clustersAverage[j]);
                //                dist = std::sqrt(std::pow((clustersAverage[j][0] - inDirs[i][0]),2) + std::pow((clustersAverage[j][1] - inDirs[i][1]),2) + std::pow((clustersAverage[j][2] - inDirs[i][2]),2));
                dist = std::sqrt((clustersAverage[j][0] - inDirs[i][0]) * (clustersAverage[j][0] - inDirs[i][0]) + (clustersAverage[j][1] - inDirs[i][1]) * (clustersAverage[j][1] - inDirs[i][1]) + (clustersAverage[j][2] - inDirs[i][2]) * (clustersAverage[j][2] - inDirs[i][2]));
                if(dist < minDist)
                {
                    minDist = dist;
                    clusters[i] = j;
                }
            }
        }

        crit = 0;
        for(int i = 0; i < numClusters; i++)
        {
            clustersAverage[i] = getNewClusterAverage(i, inDirs, clusters);
        }

        for(int i = 0; i < numDirs; i++)
            crit += clusters[i] - clustersOld[i];

        if(crit == 0)
            stop = true;

        ++It;
    }

    mainDirs.resize(numClusters);
    for(int i = 0; i < numClusters; i++)
        mainDirs[i] = clustersAverage[i];


}


typename TODEstimatorImageFilter::DirType TODEstimatorImageFilter ::getNewClusterAverage(int numCluster, DirVectorType &indirs, std::vector<int> &cluster)
{
    DirType sum, sumNorm;
    DirVectorType dirs;
    dirs = indirs;
    for(int j = 0; j < 3; j++)
        sum[j] = 0;

    int size = 0;
    for(int i = 0 ; i < dirs.size(); i++)
    {
        if(cluster[i] == numCluster)
        {
            for(int j = 0; j < 3; j++)
                sum[j] = sum[j] + dirs[i][j];
            ++size;
        }
    }
    for(int j = 0; j < 3; j++)
        sum[j] = sum[j]/size;


    anima::Normalize(sum, sumNorm);
    return sumNorm;
}


double TODEstimatorImageFilter::getEuclideanDistance(DirType dir1, DirType dir2)
{
    return std::sqrt(std::pow((dir2[0] - dir1[0]),2) + std::pow((dir2[1] - dir1[1]),2) + std::pow((dir2[2] - dir1[2]),2));
}


void TODEstimatorImageFilter::getSHCoef(DirType tmpDir, VectorType &coefs)
{
    int pos = 0, T = this->GetOutput()->GetVectorLength();
    VectorType tmp, resSH, rotatedModel;
    tmp.SetSize(T);
    DirType dir;
    this->GetOutput()->TransformLocalVectorToPhysicalVector(tmpDir, dir);
    //    this->GetOutput()->TransformPhysicalVectorToLocalVector(tmpDir, dir);
    //    dir[0] *= -1;
    //        dir[1] *= -1;
    //    dir[2] *= -1;

    for(int k = 0; k <= m_LOrder; k+=2)
    {
        for(int m = -k; m <= k; m ++)
        {
            tmp[pos] = m_GaussCoefs[pos];
            pos++;
        }
    }

    //    if(this->GetOutput()->GetDirection()[0][0] > 0)
    //        dir[0] *= -1;
    //    if(this->GetOutput()->GetDirection()[1][1] < 0)
    //        dir[1] *= -1;
    //    if(this->GetOutput()->GetDirection()[2][2] < 0)
    //        dir[2] *= -1;
    //    itk::Matrix<double> itkRotationMatrix = anima::GetRotationMatrixFromVectors(m_CstDir, dir);
    itk::Matrix<double> itkRotationMatrix = anima::GetRotationMatrixFromVectors(dir, m_CstDir);

    vnl_matrix<double> rotationMatrix;
    rotationMatrix.set_size(3,3);
    for(int x = 0; x < 3; ++x)
        for(int y = 0; y < 3; ++y)
            rotationMatrix.put(x, y, itkRotationMatrix[y][x]);

    //        vnl_matrix<double> rotationMatrix = this->GetRotationMatrix(CstDir, dir);
    std::vector<double> eulerAngles;
    anima::GetEulerAnglesFromRotationMatrix(rotationMatrix, eulerAngles);

    rotatedModel = tmp;

    vnl_matrix<double> ODFRotationMatrix;
    for (unsigned int l = 0;l <= m_LOrder;l += 2)
    {
        anima::EstimateLocalODFRotationMatrix(ODFRotationMatrix, l, eulerAngles[0],eulerAngles[1],eulerAngles[2]);

        unsigned int mBaseInd = (l*l + l + 2)/2 - l - 1;
        for (unsigned int m = 0;m <= 2*l;++m)
        {
            rotatedModel[mBaseInd + m] = 0;

            for (unsigned int mp = 0;mp <= 2*l;++mp)
                rotatedModel[mBaseInd + m] += ODFRotationMatrix(m,mp)*tmp[mBaseInd + mp];
        }
    }
    coefs = rotatedModel;

}

void TODEstimatorImageFilter::averageODFs(std::vector<VectorType> &vecCoefs, VectorType &resOdf)
{
    int numImages = vecCoefs.size();

    std::vector<std::vector<double>> odfHistos(numImages);
    std::vector<VectorType> vecSqrtCoefs(numImages);
    std::vector<double> avgSqrtHisto(m_NbSample);
    VectorType avgSqrtCoefs(m_VectorLength);

    for(int i = 0; i < numImages; i++)
    {
        odfHistos[i].resize(m_NbSample);

        this->discretizeODF(vecCoefs[i], odfHistos[i]);
        vecSqrtCoefs[i] = getSquareRootODF(odfHistos[i]);
    }

    this->getAverageCoefs(vecSqrtCoefs, avgSqrtCoefs);
    this->discretizeODF(avgSqrtCoefs, avgSqrtHisto);
    resOdf = this->getSquareODF(avgSqrtHisto);
}

void TODEstimatorImageFilter::precomputeSH()
{
    int k = 0;

    std::vector<double> tmpDir(2), revDir(3);
    for(int i = 0; i < m_NbSample; i++)
    {
        anima::TransformCartesianToSphericalCoordinates(m_SphereSampl[i], tmpDir);
        int c = 0;
        for(double l = 0; l <= m_LOrder; l+=2)
        {
            for(double m = -l; m <= l; m++)
            {
                m_SpherHarm(k, c++) = m_ODFSHBasis->getNthSHValueAtPosition(l, m, tmpDir[0], tmpDir[1]);
            }
        }
        k++;
    }
}

void TODEstimatorImageFilter::discretizeODF(VectorType ODFCoefs, std::vector<double> &ODFDiscret)
{
    int k = 0;
    double resVal2 = 0;


    std::vector<double> tmpDir(2), revDir(3);
    for(int i = 0; i < m_NbSample; i++)
    {
        anima::TransformCartesianToSphericalCoordinates(m_SphereSampl[i], tmpDir);
        resVal2 = m_ODFSHBasis->getValueAtPosition(ODFCoefs, tmpDir[0], tmpDir[1]);

        ODFDiscret[i] = resVal2;
        k++;
    }
}



typename TODEstimatorImageFilter::VectorType TODEstimatorImageFilter::getSquareRootODF(std::vector<double> ODFDiscret)
{
    vnl_matrix<double> Odf(m_NbSample, 1);
    vnl_matrix<double> Coef(m_VectorLength, 1);

    for(int i = 0; i < m_NbSample; ++i)
        Odf(i, 0) = std::sqrt(ODFDiscret[i]);

    Coef = vnl_matrix_inverse<double>(m_SpherHarm.transpose() * m_SpherHarm).as_matrix() * m_SpherHarm.transpose()*Odf;

    VectorType ModelValue(m_VectorLength);
    for(int i = 0; i < m_VectorLength; ++i)
    {
        ModelValue[i] = Coef(i, 0);
    }
    return ModelValue;
}


typename TODEstimatorImageFilter::VectorType
TODEstimatorImageFilter::getSquareODF(std::vector<double> ODFDiscret)
{
    vnl_matrix<double> Odf(m_NbSample, 1);
    vnl_matrix<double> Coef(m_VectorLength, 1);

    for(int i = 0; i < m_NbSample; ++i)
        Odf(i, 0) = ODFDiscret[i] * ODFDiscret[i];

    Coef = vnl_matrix_inverse<double>(m_SpherHarm.transpose() * m_SpherHarm).as_matrix() * m_SpherHarm.transpose()*Odf;

    VectorType ModelValue(m_VectorLength);
    for(int i = 0; i < m_VectorLength; ++i)
    {
        ModelValue[i] = Coef(i, 0);
    }

    ModelValue[0] = 1/(2*sqrt(M_PI));
    return ModelValue;
}


void
TODEstimatorImageFilter::getAverageCoefs(std::vector<VectorType> &vecCoefs, VectorType &avgCoef)
{
    int numImage = vecCoefs.size();

    std::vector<std::vector<double>> arrayCoef(numImage), arrayLogMap(numImage);
    std::vector<double> expMap(m_VectorLength), mean(m_VectorLength), nextMean(m_VectorLength), tangent(m_VectorLength);

    const int T = 100;

    for(int i = 0; i < numImage; i++)
    {
        arrayCoef[i].resize(m_VectorLength);
        arrayLogMap[i].resize(m_VectorLength);

        for(int n = 0; n < m_VectorLength; n++)
        {
            arrayCoef[i][n] = vecCoefs[i][n];
            //            nextMean[n] += 1/(double)numImage * arrayCoef[i][n];
        }
    }

    nextMean = arrayCoef[0];

    for(int t = 0; t < T; t++)
    {
        mean = nextMean;

        for(int i = 0; i < numImage; i++)
            anima::sphere_log_map(arrayCoef[i], mean, arrayLogMap[i]);

        std::fill(tangent.begin(), tangent.end(), 0);
        for(int n = 0; n < m_VectorLength; n++)
            for(int i = 0; i < numImage; i++)
                tangent[n] += 1/(double)numImage * arrayLogMap[i][n];

        anima::sphere_exp_map(tangent, mean, nextMean);
    }

    //    nextMean[0] = 2.82094792e-01;
    for(int i = 0; i < m_VectorLength; ++i)
        avgCoef[i] = nextMean[i];
}



void TODEstimatorImageFilter::processFiber(FiberType &fiber, baseSH &basis)
{
    const int T = (m_LOrder+1)*(m_LOrder+2)/2;
    VectorType tmp, resSH, rotatedModel;
    tmp.SetSize(T);
    resSH.SetSize(T);
    rotatedModel.SetSize(T);
    //    itk::Vector<float, 3> tmp, resSH;
    DirType CstDir;
    CstDir[0] = 0;
    CstDir[1] = 0;
    CstDir[2] = 1;


    TOutputImage *output = this->GetOutput();
    int nbPoints = fiber.size();

    for(int i = 0; i < nbPoints; i++)
    {
        DirType dir = getFiberDirection(i, fiber);
        //        OutputImagePointer::DirectionType refDir = output->GetDirection();

        if(output->GetDirection()[0][0] > 0)
            dir[0] *= -1;
        if(output->GetDirection()[1][1] < 0)
            dir[1] *= -1;
        if(output->GetDirection()[2][2] < 0)
            dir[2] *= -1;



        int pos = 0;
        for(int k = 0; k <= m_LOrder; k+=2)
        {
            for(int m = -k; m <= k; m ++)
            {
                tmp[pos] = m_GaussCoefs[pos];
                pos++;
            }
        }

        itk::Matrix<double> itkRotationMatrix = anima::GetRotationMatrixFromVectors(CstDir, dir);
        vnl_matrix<double> rotationMatrix;
        rotationMatrix.set_size(3,3);
        for(int x = 0; x < 3; ++x)
            for(int y = 0; y < 3; ++y)
                rotationMatrix.put(x, y, itkRotationMatrix[y][x]);

        //        vnl_matrix<double> rotationMatrix = this->GetRotationMatrix(CstDir, dir);
        std::vector<double> eulerAngles;
        anima::GetEulerAnglesFromRotationMatrix(rotationMatrix, eulerAngles);

        rotatedModel = tmp;

        vnl_matrix<double> ODFRotationMatrix;
        for (unsigned int l = 0;l <= m_LOrder;l += 2)
        {
            anima::EstimateLocalODFRotationMatrix(ODFRotationMatrix, l, eulerAngles[0],eulerAngles[1],eulerAngles[2]);

            unsigned int mBaseInd = (l*l + l + 2)/2 - l - 1;
            for (unsigned int m = 0;m <= 2*l;++m)
            {
                rotatedModel[mBaseInd + m] = 0;

                for (unsigned int mp = 0;mp <= 2*l;++mp)
                    rotatedModel[mBaseInd + m] += ODFRotationMatrix(m,mp)*tmp[mBaseInd + mp];
            }
        }
        PointType point = getCenterVoxel(i, fiber);
        itk::Index<3> index;

        output->TransformPhysicalPointToIndex(point, index);
        //        index[0] *= -1;
        //        index[1] *= -1;
        //        resSH = (output->GetPixel(index) + tmp);
        //        resSH /= resSH[0];
        output->SetPixel(index, rotatedModel);

        resSH.Fill(0);
        tmp.Fill(0);
        rotatedModel.Fill(0);

    }
}


typename TODEstimatorImageFilter::DirType
TODEstimatorImageFilter::getFiberDirection(int index, FiberType &fiber)
{
    DirType resDir;
    if(index == 0)
    {
        for(int i = 0; i < 3; i++)
            resDir[i] = fiber[index+1][i] - fiber[index][i];
    }
    else if(index == fiber.size() - 1)
    {
        for(int i = 0; i < 3; i++)
            resDir[i] = fiber[index][i] - fiber[index-1][i];
    }
    else
    {
        for(int i = 0; i < 3; i++)
            resDir[i] = fiber[index-1][i] - fiber[index+1][i];
    }

    return resDir;
}


typename TODEstimatorImageFilter::PointType
TODEstimatorImageFilter::getCenterVoxel(int index, FiberType &fiber)
{
    PointType resPoint;
    for(int i = 0; i < 3; i++)
    {
        resPoint[i] = fiber[index][i];
    }

    return resPoint;
}


typename TODEstimatorImageFilter::FiberType TODEstimatorImageFilter::readFiber(vtkIdType numberOfPoints, const vtkIdType *indices, vtkPoints *points)
{
    FiberType fiber;
    //    fiber.resize(numberOfPoints);
    for (vtkIdType i = 0; i < numberOfPoints; i++)
    {
        double tmpPoint[3];
        points->GetPoint(indices[i], tmpPoint);
        PointType point;
        for(int i = 0; i < 3; i++)
            point[i] = tmpPoint[i];

        fiber.push_back(tmpPoint);
    }

    return fiber;
}


vnl_matrix <double> TODEstimatorImageFilter::GetRotationMatrix(DirType dir1, DirType dir2)
{
    DirType nu;
    anima::ComputeCrossProduct(dir2, dir1, nu);

    double s = anima::ComputeNorm(nu);
    double c = anima::ComputeScalarProduct(dir2, dir1);

    vnl_matrix <double> R;
    vnl_matrix <double> MatNu;
    vnl_matrix <double> SqrMatNu;
    MatNu.set_size(3,3);
    MatNu.fill(0.0);

    MatNu.put(0,1,-nu[2]);
    MatNu.put(0,2, nu[1]);
    MatNu.put(1,0, nu[2]);
    MatNu.put(1,2,-nu[0]);
    MatNu.put(2,0,-nu[1]);
    MatNu.put(2,1, nu[0]);


    //    MatNu[1][0] = -nu[2];
    //    MatNu[2][0] = nu[1];
    //    MatNu[0][1] = nu[2];
    //    MatNu[2][1] = -nu[0];
    //    MatNu[0][2] = -nu[1];
    //    MatNu[1][2] = nu[0];

    SqrMatNu = MatNu*MatNu;

    R.set_size(3,3);
    R.set_identity();
    R = R + MatNu + SqrMatNu*(1-c)/(s*s);
    //    R = R + MatNu + SqrMatNu*(1/1+c);

    return R;

}

} // end namespace anima

