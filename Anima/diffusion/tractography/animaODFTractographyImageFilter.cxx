#include "animaODFTractographyImageFilter.h"
#include <cmath>
#include <random>

#include <animaODFMaximaCostFunction.h>
#include <animaNLOPTOptimizers.h>

#include <animaVectorOperations.h>
#include <animaMatrixOperations.h>
#include <animaDistributionSampling.h>
#include <animaVMFDistribution.h>
#include <animaWatsonDistribution.h>

namespace anima
{

  ODFTractographyImageFilter::ODFTractographyImageFilter()
      : BaseTractographyImageFilter()
  {
      m_ODFSHOrder = 4;
      m_GFAThreshold = 0.1;

      m_ODFSHBasis = NULL;

      this->SetModelDimension(15);

  }

  ODFTractographyImageFilter::~ODFTractographyImageFilter()
  {
      if (m_ODFSHBasis)
          delete m_ODFSHBasis;
  }

  void ODFTractographyImageFilter::PrepareTractography()
  {
      // Call base preparation
      BaseTractographyImageFilter::PrepareTractography();

      m_ODFSHOrder = std::round(-1.5 + 0.5 * std::sqrt(8 * this->GetInputImage()->GetNumberOfComponentsPerPixel() + 1));
      this->SetModelDimension((m_ODFSHOrder + 1)*(m_ODFSHOrder + 2)/2);

      // Initialize estimation matrices for Aganj et al based estimation
      if (m_ODFSHBasis)
          delete m_ODFSHBasis;

      m_ODFSHBasis = new anima::ODFSphericalHarmonicBasis(m_ODFSHOrder);
  }

  void ODFTractographyImageFilter::SetInputImage(ModelImageType *input)
  {
      this->Superclass::SetInputImage(input);

      m_interpolator = InterpolatorType::New();
      m_interpolator->SetInputImage(input);
  }


  bool ODFTractographyImageFilter::CheckModelCompatibility(VectorType &modelValue, itk::ThreadIdType threadId)
  {
      double fractionalAnisotropy = this->GetGeneralizedFractionalAnisotropy(modelValue);
      if (fractionalAnisotropy < m_GFAThreshold)
          return false;

      return true;
  }

  bool ODFTractographyImageFilter::CheckIndexInImageBounds(ContinuousIndexType &index)
  {
    return m_interpolator->IsInsideBuffer(index);
  }

  unsigned int ODFTractographyImageFilter::FindODFMaxima(const VectorType &modelValue, DirectionVectorType &maxima, double minVal, bool is2d)
  {
      ListType modelValueList(modelValue.GetSize());
      for (unsigned int i = 0;i < modelValue.GetSize();++i)
          modelValueList[i] = modelValue[i];

      std::vector < std::vector <double> > initDirs(15);

      // Find the max of the ODF v, set max to the value...
      struct XYZ array[] = {
      {-1,0,0},
      {0,-1,0},
      {0,0,-1},
      {0.1789,0.1113,-0.9776},
      {0.0635,-0.3767,-0.9242},
      {-0.7108,-0.0516,-0.7015},
      {-0.6191,0.4385,-0.6515},
      {-0.2424,-0.7843,-0.571},
      {0.2589,0.618,-0.7423},
      {0.8169,-0.1697,-0.5513},
      {0.8438,-0.5261,-0.106},
      {0.2626,-0.9548,-0.1389},
      {-1e-04,-0.9689,0.2476},
      {-0.7453,-0.6663,0.0242},
      {-0.9726,-0.2317,0.0209}
  };

      for (unsigned int i = 0;i < initDirs.size();++i)
      {
          initDirs[i].resize(3);
          initDirs[i][0] = array[i].x;
          initDirs[i][1] = array[i].y;
          initDirs[i][2] = array[i].z;
      }

      typedef anima::ODFMaximaCostFunction CostFunctionType;
      typedef anima::NLOPTOptimizers OptimizerType;

      OptimizerType::Pointer opt = OptimizerType::New();
      opt->SetAlgorithm(NLOPT_LN_BOBYQA);
      opt->SetXTolRel(1.0e-4);
      opt->SetFTolRel(1.0e-6);
      opt->SetMaxEval(200);
      opt->SetVectorStorageSize(2000);

      typedef std::map < double, Vector3DType > MapType;
      MapType dmap;

      std::vector <double> angle;
      OptimizerType::ParametersType tmpValue(2);
      Vector3DType tmpValueVector(1.0);
      Vector3DType cartesianVector(1.0);

      itk::Array<double> lowerBounds(2);
      itk::Array<double> upperBounds(2);

      lowerBounds.fill(0.0);
      upperBounds[0] = M_PI;
      upperBounds[1] = 2.0 * M_PI;

      opt->SetLowerBoundParameters(lowerBounds);
      opt->SetUpperBoundParameters(upperBounds);

      for (unsigned int i = 0; i < initDirs.size();++i)
      {
          anima::TransformCartesianToSphericalCoordinates(initDirs[i],angle);
          tmpValue[0] = angle[0];
          tmpValue[1] = angle[1];

          CostFunctionType::Pointer cost = CostFunctionType::New();
          cost->SetODFSHOrder(m_ODFSHOrder);
          cost->SetBasisParameters(modelValueList);

          opt->SetCostFunction(cost);
          opt->SetMaximize(true);

          opt->SetInitialPosition(tmpValue);
          opt->StartOptimization();

          tmpValue = opt->GetCurrentPosition();
          tmpValueVector[0] = tmpValue[0];
          tmpValueVector[1] = tmpValue[1];

          // Some check needed here to see if we really found a maximum
          anima::TransformSphericalToCartesianCoordinates(tmpValueVector,cartesianVector);
          dmap[opt->GetValue()] = cartesianVector;
      }

      // Find true maximas
      std::vector <bool> usefulMaxima(15,true);
      unsigned int pos = 0;
      for (MapType::reverse_iterator it = dmap.rbegin();it != dmap.rend();++it)
      {
          unsigned int posIn = pos+1;
          MapType::reverse_iterator itin = it;
          ++itin;
          for (;itin != dmap.rend();++itin)
          {
              if (anima::ComputeOrientationAngle((*it).second,(*itin).second) < 15)
                  usefulMaxima[posIn] = false;

              ++posIn;
          }

          ++pos;
      }

      pos = 0;
      maxima.clear();
      for (MapType::reverse_iterator it = dmap.rbegin();it != dmap.rend();++it)
      {
          if ((usefulMaxima[pos])&&((*it).first > minVal))
              maxima.push_back((*it).second);

          ++pos;
      }

      if (is2d)
      {
          std::vector <bool> outOfPlaneDirs(maxima.size(),false);
          for (unsigned int i = 0;i < maxima.size();++i)
          {
              maxima[i][2] = 0;

              double norm = 0;
              for (unsigned int j = 0;j < ModelImageType::ImageDimension - 1;++j)
                  norm += maxima[i][j] * maxima[i][j];
              norm = sqrt(norm);

              outOfPlaneDirs[i] = (std::abs(norm) < 0.5);

              if (!outOfPlaneDirs[i])
              {
                  for (unsigned int j = 0;j < ModelImageType::ImageDimension - 1;++j)
                      maxima[i][j] /= norm;
              }
          }

          DirectionVectorType outMaxima;

          for (unsigned int i = 0;i < maxima.size();++i)
          {
              if (!outOfPlaneDirs[i])
                  outMaxima.push_back(maxima[i]);
          }

          maxima = outMaxima;
      }

      return maxima.size();
  }

  void ODFTractographyImageFilter::GetModelValue(ContinuousIndexType &index, VectorType &modelValue)
  {
//    modelValue.SetSize(this->GetModelDimension());
//    modelValue.Fill(0.0);

//    modelValue = m_interpolator->EvaluateAtContinuousIndex(index);
//    for(int i =0; i < modelValue.GetNumberOfElements(); i++)
//        auto test = modelValue.GetElement(i);
//    IndexType ind;
//    ind[0] = index[0];
//    ind[1] = index[1];
//    ind[2] = index[2];
//    modelValue = this->GetInputImage()->GetPixel(ind);
    if (m_interpolator->IsInsideBuffer(index))
        modelValue = m_interpolator->EvaluateAtContinuousIndex(index);
  }


  double ODFTractographyImageFilter::GetGeneralizedFractionalAnisotropy(VectorType &modelValue)
  {
      double sumSquares = 0;
      for (unsigned int i = 0;i < this->GetModelDimension();++i)
          sumSquares += modelValue[i]*modelValue[i];

      return std::sqrt(1.0 - modelValue[0]*modelValue[0]/sumSquares);

  }

  ODFTractographyImageFilter::PointType ODFTractographyImageFilter::GetModelPrincipalDirection(VectorType &modelValue, bool is2d, itk::ThreadIdType threadId)
  {
    PointType resVec(0.0);

    bool is2d_ = (this->GetInputImage()->GetLargestPossibleRegion().GetSize()[2] == 1);

    DirectionVectorType maximaODF;
    unsigned int numDirs = this->FindODFMaxima(modelValue,maximaODF,m_MinimalDiffusionProbability,is2d_);

    if (numDirs==0)
        return resVec;

    for (unsigned int j = 0;j < ModelImageType::ImageDimension;++j)
        resVec[j] = maximaODF[0][j];

    if (is2d_)
    {
        resVec[2] = 0;
//        resVec.Normalize();
    }
    return resVec;
  }

  ODFTractographyImageFilter::PointType ODFTractographyImageFilter::GetNextDirection(PointType &previousDirection, VectorType &modelValue, bool is2d, itk::ThreadIdType threadId)
  {
    PointType resVec(0.0), tmpVec;

    bool is2d_ = (this->GetInputImage()->GetLargestPossibleRegion().GetSize()[2] == 1);

    DirectionVectorType maximaODF;
    unsigned int numDirs = this->FindODFMaxima(modelValue,maximaODF,m_MinimalDiffusionProbability,is2d_);
    if (numDirs == 0)
        return previousDirection;

    double maxVal = 0;
    for (unsigned int i = 0;i < numDirs;++i)
    {
        for (unsigned int j = 0;j < ModelImageType::ImageDimension;++j)
            tmpVec[j] = maximaODF[i][j];

        double tmpVal = anima::ComputeScalarProduct(previousDirection, tmpVec);

        if (tmpVal < 0)
        {
            for (unsigned int i = 0; i < ModelImageType::ImageDimension; ++i)
                tmpVec[i] *= -1;
            tmpVal *= -1;
        }

        if (tmpVal > maxVal)
        {
            resVec = tmpVec;
            maxVal = tmpVal;
        }
    }

    if (is2d_)
    {
        resVec[2] = 0;
//        resVec.Normalize();
    }
    return resVec;
  }

}
