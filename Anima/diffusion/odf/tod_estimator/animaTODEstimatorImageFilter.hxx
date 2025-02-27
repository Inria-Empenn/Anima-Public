#pragma once

#include "animaTODEstimatorImageFilter.h"
#include <animaLogExpMapsUnitSphere.h>
#include <animaMatrixOperations.h>
#include <animaODFFunctions.h>
#include <animaShapesReader.h>
#include <animaSphereOperations.h>

#include <vtkCellArrayIterator.h>

#include <random>

namespace anima
{
    template <typename ScalarType>
    void
    TODEstimatorImageFilter<ScalarType>::BeforeThreadedGenerateData()
    {
        Superclass::BeforeThreadedGenerateData();

        ReferenceImagePointerType refImg = anima::readImage<ReferenceImageType>(m_ReferenceFileName);
        m_LOrder = 8;
        unsigned int vectorLength = (m_LOrder + 1) * (m_LOrder + 2) / 2;
        m_VectorLength = vectorLength;
        OutputImageType *output = this->GetOutput();
        output->SetVectorLength(vectorLength);
        output->SetRegions(refImg->GetLargestPossibleRegion());
        output->SetSpacing(refImg->GetSpacing());
        output->SetDirection(refImg->GetDirection());
        output->SetOrigin(refImg->GetOrigin());
        output->Allocate();

        OutputImagePixelType tmp;
        tmp.SetSize((m_LOrder + 1) * (m_LOrder + 2) / 2);
        tmp.Fill(0);
        output->FillBuffer(tmp);

        m_CstDir[0] = 0;
        m_CstDir[1] = 0;
        m_CstDir[2] = 1;

        m_NbSample = 200;
        anima::GetSphereEvenSampling(m_SphereSampl, m_NbSample);

        m_SpherHarm.set_size(m_NbSample, vectorLength);

        m_ODFSHBasis = new anima::ODFSphericalHarmonicBasis(m_LOrder);
        this->PrecomputeSH();
        this->ComputeCoefs();

        const int Xmax = output->GetLargestPossibleRegion().GetSize()[0];
        const int Ymax = output->GetLargestPossibleRegion().GetSize()[1];
        const int Zmax = output->GetLargestPossibleRegion().GetSize()[2];
        m_ImgDir.resize(Xmax * Ymax * Zmax);

        anima::ShapesReader trackReader;
        trackReader.SetFileName(m_InputFileName);
        trackReader.Update();

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
            lineCount++;

            FiberType fiber = ReadFiber(numberOfPoints, indices, points);
            int nbPoints = fiber.size();

            for (int i = 0; i < nbPoints; i++)
            {
                DirType dir = GetFiberDirection(i, fiber);
                if(output->GetDirection()[0][0] > 0)
                {
                    dir[0] *= -1;
                }
                if(output->GetDirection()[1][1] < 0)
                {
                    dir[1] *= -1;
                }
                if(output->GetDirection()[2][2] < 0)
                {
                    dir[2] *= -1;
                }


                PointType point = GetCenterVoxel(i, fiber);
                itk::Index<3> index;

                output->TransformPhysicalPointToIndex(point, index);

                m_ImgDir[index[0] + Xmax * (index[1] + Ymax * index[2])].push_back(dir);
            }
        }
    }

    template <typename ScalarType>
    void
    TODEstimatorImageFilter<ScalarType>::ComputeCoefs()
    {
        m_GaussCoefs.SetSize(91);
        m_GaussCoefs.Fill(0.0);

        double l1 = 1.71e-4;
        double l2 = 2e-5;
        double d = l1 - l2;
        double spi = sqrt(M_PI);
        double pi = M_PI;
        double at = M_PI + 2 * atan((d - l2) / (2 * sqrt(d * l2)));

        double c0 = 1 / (2 * sqrt(pi));
        double c1 = -sqrt(5) * (3 * sqrt(l2) * (d + l2) * (pi + atan((d - l2) / (2 * sqrt(l2 * d)))) - 4 * sqrt(d) * (2 * d + 3 * l2)) / (16 * sqrt(pi) * pow(d, 1.5));
        double c2 = 3 * (4 * sqrt(d) * (16 * pow(d, 2) + 115 * d * l2 + 105 * pow(l2, 2)) - 15 * sqrt(l2) * (3 * pow(d, 2) + 10 * d * l2 + 7 * pow(l2, 2)) * (pi + 2 * atan((d - l2) / (2 * sqrt(d * l2))))) / (128 * spi * pow(d, 2.5));
        double c3 = -105 * sqrt(13) * (at * sqrt(l2) * (5 * pow(d, 3) + 35 * l2 * pow(d, 2) + 63 * pow(l2, 2) * d + 33 * pow(l2, 3)) - 4 * sqrt(d) * (33 * pow(l2, 3) + 52 * d * pow(l2, 2) + 103 * pow(d, 2) * l2 / 5 + 128 * pow(d, 3) / 105)) / (1024 * spi * pow(d, 3.5));
        double c4 = -315 * sqrt(17) * (at * sqrt(l2) * (35 * pow(d, 4) + 420 * pow(d, 3) * l2 + 1386 * pow(d, 2) * pow(l2, 2) + 1716 * pow(l2, 3) * d + 715 * pow(l2, 4)) - 4 * sqrt(d) * (715 * pow(l2, 4) + 4433 * pow(l2, 3) * d / 3 + 957 * pow(l2, 2) * pow(d, 2) + 6967 * l2 * pow(d, 3) / 35 + 2048 * pow(d, 4) / 315)) / (16384 * spi * pow(d, 4.5));
        double c5 = -24255 * sqrt(21) * (at * sqrt(l2) * (9 * pow(d, 5) + 165 * pow(d, 4) * l2 + 858 * pow(d, 3) * pow(l2, 2) + 12870 * pow(d, 2) * pow(l2, 3) / 7 + 12155 * d * pow(l2, 4) / 7 + 4199 * pow(l2, 5) / 7) - 4 * sqrt(d) * (4199 * pow(l2, 5) / 7 + 32266 * pow(l2, 4) * d / 21 + 20691 * pow(l2, 3) * pow(d, 2) / 15 + 24830 * pow(l2, 2) * pow(d, 3) / 49 + 28799 * l2 * pow(d, 4) / 441 + 32768 * pow(d, 5) / 24255)) / (262144 * pow(d, 5.5) * spi);

        m_GaussCoefs[0] = c0;
        m_GaussCoefs[3] = c1;
        m_GaussCoefs[10] = c2;
        m_GaussCoefs[21] = c3;
        m_GaussCoefs[36] = c4;
        m_GaussCoefs[55] = c5;
    }

    template <typename ScalarType>
    void
    TODEstimatorImageFilter<ScalarType>::DynamicThreadedGenerateData(const TODEstimatorImageFilter::OutputImageRegionType &outputRegionForThread)
    {
        OutputImageType *output = this->GetOutput();
        itk::ImageRegionIterator<OutputImageType> outItr(output, outputRegionForThread);
        typename OutputImageType::IndexType index, index_out;

        OutputImagePixelType nullVector;
        nullVector.SetSize((m_LOrder + 1) * (m_LOrder + 2) / 2);
        nullVector.Fill(0);

        const int Xmax = output->GetLargestPossibleRegion().GetSize()[0];
        const int Ymax = output->GetLargestPossibleRegion().GetSize()[1];
        const int Zmax = output->GetLargestPossibleRegion().GetSize()[2];

        while (!outItr.IsAtEnd())
        {
            index = outItr.GetIndex();

            std::vector<OutputImagePixelType> vecCoefs;
            int tmpSize = m_ImgDir[index[0] + Xmax * (index[1] + Ymax * index[2])].size();
            if (tmpSize != 0)
            {
                DirVectorType mainDirs;
                DirVectorType localDir(tmpSize);
                localDir = m_ImgDir[index[0] + Xmax * (index[1] + Ymax * index[2])];
                this->GetMainDirections(localDir, mainDirs);

                OutputImagePixelType tmpCoefs(m_VectorLength);
                tmpCoefs.Fill(0);
                for (int i = 0; i < mainDirs.size(); i++)
                {
                    if (!isnan(mainDirs[i][0]))
                    {
                        this->GetSHCoef(mainDirs[i], tmpCoefs);

                        vecCoefs.push_back(tmpCoefs);
                    }
                }

                OutputImagePixelType resCoefs(m_VectorLength);
                resCoefs.Fill(0);
                if (vecCoefs.size() != 0)
                {
                    resCoefs = vecCoefs[0];
                    this->AverageODFs(vecCoefs, resCoefs);
                    if (m_UseNormalization)
                    {
                        for (int i = 0; i < m_VectorLength; i++)
                            resCoefs[i] /= vecCoefs.size();
                    }
                    outItr.Set(resCoefs);
                }

                else
                    outItr.Set(nullVector);

                std::vector<double> testSh(m_VectorLength);
                for (int i = 0; i < m_VectorLength; i++)
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

    template <typename ScalarType>
    void
    TODEstimatorImageFilter<ScalarType>::GetMainDirections(DirVectorType inDirs, DirVectorType &mainDirs)
    {

        if (inDirs.size() == 0)
        {
            mainDirs.resize(0);
            return;
        }
        const int numIt = 200;
        int numDirs = inDirs.size(), crit = 0, numClusters = 0;

        if (numDirs < 4)
        {
            numClusters = numDirs;
            mainDirs.resize(numClusters);
            for (int i = 0; i < numClusters; i++)
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

        for (int i = 0; i < numClusters; i++)
        {
            int k = distrib(gen);
            clustersAverage[i] = inDirs[k];
        }

        int It = 0;
        bool stop = false;
        double dist = 0;
        while (It < numIt && !stop)
        {
            std::vector<int> clustersOld(numDirs);
            for (int i = 0; i < numDirs; i++)
                clustersOld[i] = clusters[i];

            for (int i = 0; i < numDirs; i++)
            {
                double minDist = 1e10;
                for (int j = 0; j < numClusters; j++)
                {
                    dist = std::sqrt((clustersAverage[j][0] - inDirs[i][0]) * (clustersAverage[j][0] - inDirs[i][0]) + (clustersAverage[j][1] - inDirs[i][1]) * (clustersAverage[j][1] - inDirs[i][1]) + (clustersAverage[j][2] - inDirs[i][2]) * (clustersAverage[j][2] - inDirs[i][2]));
                    if (dist < minDist)
                    {
                        minDist = dist;
                        clusters[i] = j;
                    }
                }
            }

            crit = 0;
            for (int i = 0; i < numClusters; i++)
            {
                clustersAverage[i] = GetNewClusterAverage(i, inDirs, clusters);
            }

            for (int i = 0; i < numDirs; i++)
                crit += clusters[i] - clustersOld[i];

            if (crit == 0)
                stop = true;

            ++It;
        }

        mainDirs.resize(numClusters);
        for (int i = 0; i < numClusters; i++)
            mainDirs[i] = clustersAverage[i];
    }

    template <typename ScalarType>
    typename TODEstimatorImageFilter<ScalarType>::DirType
    TODEstimatorImageFilter<ScalarType>::GetNewClusterAverage(int numCluster, DirVectorType &indirs, std::vector<int> &cluster)
    {
        DirType sum, sumNorm;
        DirVectorType dirs;
        dirs = indirs;
        for (int j = 0; j < 3; j++)
            sum[j] = 0;

        int size = 0;
        for (int i = 0; i < dirs.size(); i++)
        {
            if (cluster[i] == numCluster)
            {
                for (int j = 0; j < 3; j++)
                    sum[j] = sum[j] + dirs[i][j];
                ++size;
            }
        }
        for (int j = 0; j < 3; j++)
            sum[j] = sum[j] / size;

        anima::Normalize(sum, sumNorm);
        return sumNorm;
    }

    template <typename ScalarType>
    double
    TODEstimatorImageFilter<ScalarType>::GetEuclideanDistance(DirType dir1, DirType dir2)
    {
        return std::sqrt(std::pow((dir2[0] - dir1[0]), 2) + std::pow((dir2[1] - dir1[1]), 2) + std::pow((dir2[2] - dir1[2]), 2));
    }

    template <typename ScalarType>
    void
    TODEstimatorImageFilter<ScalarType>::GetSHCoef(DirType tmpDir, OutputImagePixelType &coefs)
    {
        int pos = 0, T = this->GetOutput()->GetVectorLength();
        OutputImagePixelType tmp, resSH, rotatedModel(m_VectorLength);
        tmp.SetSize(T);
        DirType dir;
        this->GetOutput()->TransformLocalVectorToPhysicalVector(tmpDir, dir);
        for (int k = 0; k <= m_LOrder; k += 2)
        {
            for (int m = -k; m <= k; m++)
            {
                tmp[pos] = m_GaussCoefs[pos];
                pos++;
            }
        }

        itk::Matrix<double> itkRotationMatrix = anima::GetRotationMatrixFromVectors(m_CstDir, dir);
        vnl_matrix<double> rotationMatrix;
        rotationMatrix.set_size(3, 3);
        for (int x = 0; x < 3; ++x)
            for (int y = 0; y < 3; ++y)
                rotationMatrix.put(x, y, itkRotationMatrix[y][x]);

        std::vector<double> eulerAngles;
        anima::GetEulerAnglesFromRotationMatrix(rotationMatrix, eulerAngles);

        double a = rotationMatrix(2, 2);
        rotatedModel = tmp;

        vnl_matrix<double> ODFRotationMatrix;
        for (unsigned int l = 0; l <= m_LOrder; l += 2)
        {
            anima::EstimateLocalODFRotationMatrix(ODFRotationMatrix, l, eulerAngles[0], eulerAngles[1], eulerAngles[2]);

            unsigned int mBaseInd = (l * l + l + 2) / 2 - l - 1;
            for (unsigned int m = 0; m <= 2 * l; ++m)
            {
                rotatedModel[mBaseInd + m] = 0;

                for (unsigned int mp = 0; mp <= 2 * l; ++mp)
                    rotatedModel[mBaseInd + m] += ODFRotationMatrix(m, mp) * tmp[mBaseInd + mp];
            }
        }
        coefs = rotatedModel;
    }

    template <typename ScalarType>
    void
    TODEstimatorImageFilter<ScalarType>::AverageODFs(std::vector<OutputImagePixelType> &vecCoefs, OutputImagePixelType &resOdf)
    {
        int numImages = vecCoefs.size();

        std::vector<std::vector<double>> odfHistos(numImages);
        std::vector<OutputImagePixelType> vecSqrtCoefs(numImages);
        std::vector<double> avgSqrtHisto(m_NbSample);
        OutputImagePixelType avgSqrtCoefs(m_VectorLength);

        for (int i = 0; i < numImages; i++)
        {
            odfHistos[i].resize(m_NbSample);

            this->DiscretizeODF(vecCoefs[i], odfHistos[i]);
            vecSqrtCoefs[i] = GetSquareRootODF(odfHistos[i]);
        }

        this->GetAverageCoefs(vecSqrtCoefs, avgSqrtCoefs);
        this->DiscretizeODF(avgSqrtCoefs, avgSqrtHisto);
        resOdf = this->GetSquareODF(avgSqrtHisto);
    }

    template <typename ScalarType>
    void
    TODEstimatorImageFilter<ScalarType>::PrecomputeSH()
    {
        int k = 0;

        std::vector<double> tmpDir(2), revDir(3);
        for (int i = 0; i < m_NbSample; i++)
        {
            anima::TransformCartesianToSphericalCoordinates(m_SphereSampl[i], tmpDir);
            int c = 0;
            for (double l = 0; l <= m_LOrder; l += 2)
            {
                for (double m = -l; m <= l; m++)
                {
                    m_SpherHarm(k, c++) = m_ODFSHBasis->getNthSHValueAtPosition(l, m, tmpDir[0], tmpDir[1]);
                }
            }
            k++;
        }
    }

    template <typename ScalarType>
    void
    TODEstimatorImageFilter<ScalarType>::DiscretizeODF(OutputImagePixelType ODFCoefs, std::vector<double> &ODFDiscret)
    {
        int k = 0;
        double resVal2 = 0;

        std::vector<double> tmpDir(2), revDir(3);
        for (int i = 0; i < m_NbSample; i++)
        {
            anima::TransformCartesianToSphericalCoordinates(m_SphereSampl[i], tmpDir);
            resVal2 = m_ODFSHBasis->getValueAtPosition(ODFCoefs, tmpDir[0], tmpDir[1]);

            ODFDiscret[i] = resVal2;
            k++;
        }
    }

    template <typename ScalarType>
    typename TODEstimatorImageFilter<ScalarType>::OutputImagePixelType
    TODEstimatorImageFilter<ScalarType>::GetSquareRootODF(std::vector<double> ODFDiscret)
    {
        vnl_matrix<double> Odf(m_NbSample, 1);
        vnl_matrix<double> Coef(m_VectorLength, 1);

        for (int i = 0; i < m_NbSample; ++i)
            Odf(i, 0) = std::sqrt(ODFDiscret[i]);

        Coef = vnl_matrix_inverse<double>(m_SpherHarm.transpose() * m_SpherHarm).as_matrix() * m_SpherHarm.transpose() * Odf;

        OutputImagePixelType ModelValue(m_VectorLength);
        for (int i = 0; i < m_VectorLength; ++i)
        {
            ModelValue[i] = Coef(i, 0);
        }
        return ModelValue;
    }

    template <typename ScalarType>
    typename TODEstimatorImageFilter<ScalarType>::OutputImagePixelType
    TODEstimatorImageFilter<ScalarType>::GetSquareODF(std::vector<double> ODFDiscret)
    {
        vnl_matrix<double> Odf(m_NbSample, 1);
        vnl_matrix<double> Coef(m_VectorLength, 1);

        for (int i = 0; i < m_NbSample; ++i)
            Odf(i, 0) = ODFDiscret[i] * ODFDiscret[i];

        Coef = vnl_matrix_inverse<double>(m_SpherHarm.transpose() * m_SpherHarm).as_matrix() * m_SpherHarm.transpose() * Odf;

        OutputImagePixelType ModelValue(m_VectorLength);
        for (int i = 0; i < m_VectorLength; ++i)
        {
            ModelValue[i] = Coef(i, 0);
        }

        ModelValue[0] = 1 / (2 * sqrt(M_PI));
        return ModelValue;
    }

    template <typename ScalarType>
    void
    TODEstimatorImageFilter<ScalarType>::GetAverageCoefs(std::vector<OutputImagePixelType> &vecCoefs, OutputImagePixelType &avgCoef)
    {
        int numImage = vecCoefs.size();

        std::vector<std::vector<double>> arrayCoef(numImage), arrayLogMap(numImage);
        std::vector<double> expMap(m_VectorLength), mean(m_VectorLength), nextMean(m_VectorLength), tangent(m_VectorLength);

        const int T = 100;

        for (int i = 0; i < numImage; i++)
        {
            arrayCoef[i].resize(m_VectorLength);
            arrayLogMap[i].resize(m_VectorLength);

            for (int n = 0; n < m_VectorLength; n++)
            {
                arrayCoef[i][n] = vecCoefs[i][n];
            }
        }

        nextMean = arrayCoef[0];

        for (int t = 0; t < T; t++)
        {
            mean = nextMean;

            for (int i = 0; i < numImage; i++)
                anima::sphere_log_map(arrayCoef[i], mean, arrayLogMap[i]);

            std::fill(tangent.begin(), tangent.end(), 0);
            for (int n = 0; n < m_VectorLength; n++)
                for (int i = 0; i < numImage; i++)
                    tangent[n] += 1 / (double)numImage * arrayLogMap[i][n];

            anima::sphere_exp_map(tangent, mean, nextMean);
        }

        //    nextMean[0] = 2.82094792e-01;
        for (int i = 0; i < m_VectorLength; ++i)
            avgCoef[i] = nextMean[i];
    }

    template <typename ScalarType>
    typename TODEstimatorImageFilter<ScalarType>::DirType
    TODEstimatorImageFilter<ScalarType>::GetFiberDirection(int index, FiberType &fiber)
    {
        DirType resDir;
        if (index == fiber.size() - 1)
        {
            for (int i = 0; i < 3; i++)
                resDir[i] = fiber[index][i] - fiber[index - 1][i];
        }
        else
        {
            for (int i = 0; i < 3; i++)
                resDir[i] = fiber[index + 1][i] - fiber[index][i];
        }

        anima::Normalize(resDir, resDir);
        return resDir;
    }

    template <typename ScalarType>
    typename TODEstimatorImageFilter<ScalarType>::PointType
    TODEstimatorImageFilter<ScalarType>::GetCenterVoxel(int index, FiberType &fiber)
    {
        PointType resPoint;
        for (int i = 0; i < 3; i++)
        {
            resPoint[i] = fiber[index][i];
        }

        return resPoint;
    }

    template <typename ScalarType>
    typename TODEstimatorImageFilter<ScalarType>::FiberType
    TODEstimatorImageFilter<ScalarType>::ReadFiber(vtkIdType numberOfPoints, const vtkIdType *indices, vtkPoints *points)
    {
        FiberType fiber;
        for (vtkIdType i = 0; i < numberOfPoints; i++)
        {
            double tmpPoint[3];
            points->GetPoint(indices[i], tmpPoint);
            PointType point;
            for (int i = 0; i < 3; i++)
                point[i] = tmpPoint[i];

            fiber.push_back(tmpPoint);
        }

        return fiber;
    }

} // end namespace anima
