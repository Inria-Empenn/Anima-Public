#pragma once
#include <animaBlockMatchInitializer.h>

#include <animaMCMImage.h>
#include <animaMultiCompartmentModel.h>

namespace anima {

template <class PixelType, unsigned int NDimensions = 3>
class MCMBlockMatchingInitializer
    : public anima::BlockMatchingInitializer<PixelType, NDimensions> {
public:
  using Superclass = anima::BlockMatchingInitializer<PixelType, NDimensions>;
  using Self = MCMBlockMatchingInitializer<PixelType, NDimensions>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  using MCMImageType = anima::MCMImage<PixelType, NDimensions>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(MCMBlockMatchingInitializer, anima::BlockMatchingInitializer);

  using MCModelType = anima::MultiCompartmentModel;
  using MCModelPointer = typename MCModelType::Pointer;

  using BlockGeneratorThreadStruct =
      typename Superclass::BlockGeneratorThreadStruct;
  using ImageRegionType = typename Superclass::ImageRegionType;

  struct MCMBlockGeneratorThreadStruct : public BlockGeneratorThreadStruct {
    std::vector<std::vector<MCModelPointer>> ref_models;
  };

  void AddReferenceImage(itk::ImageBase<NDimensions> *refImage) ITK_OVERRIDE;

protected:
  MCMBlockMatchingInitializer() : Superclass() {}

  virtual ~MCMBlockMatchingInitializer() {}

  void InitializeThreading(unsigned int maskIndex,
                           BlockGeneratorThreadStruct *&workStr) ITK_OVERRIDE;
  bool CheckOrientedModelVariance(unsigned int imageIndex,
                                  ImageRegionType &region,
                                  double &blockVariance,
                                  BlockGeneratorThreadStruct *workStr,
                                  unsigned int threadId) ITK_OVERRIDE;

private:
  MCMBlockMatchingInitializer(const Self &); // purposely not implemented
  void operator=(const Self &);              // purposely not implemented

  std::vector<MCModelPointer> m_ReferenceModels;
};

} // end of namespace anima

#include "animaMCMBlockMatchInitializer.hxx"
