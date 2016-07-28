#pragma once
#include "animaMCMImage.h"

namespace anima
{

template <typename TPixel, unsigned int VImageDimension>
MCMImage <TPixel,VImageDimension>::MCMImage()
{
    m_DescriptionModel = 0;
}

template <typename TPixel, unsigned int VImageDimension>
void
MCMImage <TPixel,VImageDimension>::SetDescriptionModel(MCMPointer &mcm)
{
    if (!mcm)
        return;

    m_DescriptionModel = mcm->Clone();
}

template <typename TPixel, unsigned int VImageDimension>
typename MCMImage <TPixel,VImageDimension>::MCMPointer &
MCMImage <TPixel,VImageDimension>::GetDescriptionModel()
{
    return m_DescriptionModel;
}

template< typename TPixel, unsigned int VImageDimension >
void
MCMImage <TPixel,VImageDimension>::Graft(const itk::DataObject *data)
{
    if(data == ITK_NULLPTR)
        return;

    Superclass::Graft(data);

    // Attempt to cast data to an Image
    const Self *constImgData = dynamic_cast <const Self *> (data);
    Self *imgData = const_cast <Self *> (constImgData);

    if (!imgData)
    {
        itkExceptionMacro( << "itk::VectorImage::Graft() cannot cast "
                           << typeid( data ).name() << " to "
                           << typeid( const Self * ).name() );
    }

    // Copy from MCMImage< TPixel, dim >
    this->SetDescriptionModel(imgData->GetDescriptionModel());
}


} // end namespace anima
