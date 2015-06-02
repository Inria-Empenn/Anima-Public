#pragma once

#include <itkImageIOBase.h>
#include <itkExceptionObject.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkVector.h>
/**
  * Retrieve the type of the components of one image from an itk::ImageIOBase,
  * then call a given function templated on that type.
  *
  * imageIO     -> An itk::ImageIOBase setted for your image file.
  * function    -> The name of the function to call: function<ComponentType>(args, args...)
  * args...     -> A list of arguments to pass to the function.
  */
#define ANIMA_RETRIEVE_COMPONENT_TYPE(imageIO, function, ...) \
switch (imageIO->GetComponentType())\
{\
    case itk::ImageIOBase::UCHAR:\
        std::cerr << "Component type detected is 'unsigned char', anima only use 'char', overflow may happen" << std::endl;\
    case itk::ImageIOBase::CHAR:\
        function<char>(__VA_ARGS__);\
        break;\
    case itk::ImageIOBase::UINT:\
        std::cerr << "Component type detected is 'unsigned int', anima only use 'unsigned short'.\n"\
                     "Undetermine behavor may happen." << std::endl;\
    case itk::ImageIOBase::USHORT:\
        function<unsigned short>(__VA_ARGS__);\
        break;\
    case itk::ImageIOBase::INT:\
        std::cerr << "Component type detected is 'int', anima only use 'short'."\
                     "Undetermine behavor may happen" << std::endl;\
    case itk::ImageIOBase::SHORT:\
        function<short>(__VA_ARGS__);\
        break;\
    case itk::ImageIOBase::FLOAT:\
        function<float>(__VA_ARGS__);\
        break;\
    case itk::ImageIOBase::DOUBLE:\
        function<double>(__VA_ARGS__);\
        break;\
    default:\
        itk::ExceptionObject excp(__FILE__, __LINE__, "Component type not supported.", ITK_LOCATION);\
        throw excp;\
}

/**
  * Retrieve the number dimension of one image from an itk::ImageIOBase,
  * then call a given function templated on that a Component type and the found number of dimension.
  *
  * imageIO         -> An itk::ImageIOBase setted for your image file.
  * ComponentType   -> The type used to templatized the function to call.
  * function        -> The name of the function to call.
  *                    function<ComponentType Dimensions>(args, args...)
  * args...         -> A list of arguments to pass to the function.
  */
#define ANIMA_RETRIEVE_NUMBER_OF_DIMENSIONS(imageIO, ComponentType, function, ...) \
switch (imageIO->GetNumberOfDimensions())\
{\
    case 2:\
        function<ComponentType, 2>(__VA_ARGS__);\
        break;\
    case 3:\
        function<ComponentType, 3>(__VA_ARGS__);\
        break;\
    case 4:\
        function<ComponentType, 4>(__VA_ARGS__);\
        break;\
    default:\
        itk::ExceptionObject excp(__FILE__, __LINE__, "Number of type not supported.", ITK_LOCATION);\
        throw excp;\
}


/**
  * Check if one image is a Scalar image or not from an itk::ImageIOBase,
  * then call a given function templated on either a itk::Image or itk::VectorImage
  *
  * imageIO         -> An itk::ImageIOBase setted for your image file.
  * ComponentType   -> The type used to templatized the function to call.
  * Dimension       -> The number of dimension on used to templatized the function yo call.
  * function        -> The name of the function to call.
  *                    function<itk::Image<ComponentType, Dimension> >(args, args...)
  *                     or
  *                    function<itk::VectorImage<ComponentType, Dimension> >(args, args...)
  * args...     -> A list of arguments to pass to the function
  */
#define ANIMA_CHECK_IF_COMPONENTS_ARE_VECTORS(imageIO, ComponentType, Dimension, function, ...) \
switch (imageIO->GetNumberOfComponents())\
{\
case 0:\
    {\
        itk::ExceptionObject excp(__FILE__, __LINE__, "Number of Component not supported.", ITK_LOCATION);\
        throw excp;\
        break;\
    }\
case 1:\
    {\
        typedef itk::Image<ComponentType, Dimension> ImageType;\
        function<ImageType>(__VA_ARGS__);\
        break;\
    }\
default:\
    {\
        typedef itk::VectorImage<ComponentType, Dimension> VectorImageType;\
        function<VectorImageType>(__VA_ARGS__);\
    }\
}

/**
  * Retrieve the number of components of one image from an itk::ImageIOBase,
  * then call a given function templated on a vector image of the right size.
  *
  * imageIO         -> An itk::ImageIOBase setted for your image file.
  * ComponentType   -> The type used to templatized the function to call.
  * Dimension       -> The number of dimension on used to templatized the function yo call.
  * function        -> The name of the function to call.
  *                    function<itk::Image<ik::Vector<ComponentType, ComponentSize>, Dimension> >(args, args...)
  * args...         -> A list of arguments to pass to the function
  */
#define ANIMA_RETRIEVE_NUMBER_OF_COMPONENTS(imageIO, ComponentType, dimension, function, ...) \
switch (imageIO->GetNumberOfComponents())\
{\
case 1:\
    function<itk::Image<ComponentType, dimension>(__VA_ARGS__);\
    break;\
case 3:\
    function<itk::Image<itk::Vector<ComponentType, 3>, 1>(__VA_ARGS__);\
case 6:\
    function<itk::Image<itk::Vector<ComponentType, 6>, 1>(__VA_ARGS__);\
default:\
    itk::ExceptionObject excp(__FILE__, __LINE__, "Number of Component not supported.", ITK_LOCATION);\
    throw excp;\
}
