#include <animaLogExpMapsUnitSphere.h>
#include <animaDistributionSampling.h>
#include <animaRealUniformDistribution.h>

#include <itkVariableLengthVector.h>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);

    TCLAP::ValueArg<unsigned int> numPointsArg("n", "nb-points", "Number of points (default: 20)", false, 20, "number of points", cmd);
    TCLAP::SwitchArg equalWeightsArg("E", "equal-weights", "Use equal weights", cmd, false);

    try
    {
        cmd.parse(argc, argv);
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    using RandomGeneratorType = std::mt19937;
    RandomGeneratorType generator(time(0));

    // This is what will actually supply drawn values.
    std::vector<std::vector<double>> random_sphere_points(numPointsArg.getValue());
    std::vector<double> random_sphere_point(3, 0);
    std::vector<double> average_point(3);
    std::vector<double> weights(numPointsArg.getValue());

    for (unsigned int i = 0; i < numPointsArg.getValue(); ++i)
    {
        anima::SampleFromUniformDistributionOn2Sphere(generator, random_sphere_point);
        random_sphere_points[i] = random_sphere_point;
    }

    if (equalWeightsArg.isSet())
    {
        for (unsigned int i = 0; i < numPointsArg.getValue(); ++i)
            weights[i] = 1.0 / numPointsArg.getValue();
    }
    else
    {
        anima::RealUniformDistribution unifDistr;
        unifDistr.SetLowerBoundParameter(0.01);
        unifDistr.SetUpperBoundParameter(100.0);
        unifDistr.Random(weights, generator);

        double sumWeights = 0;
        for (unsigned int i = 0; i < numPointsArg.getValue(); ++i)
            sumWeights += weights[i];

        for (unsigned int i = 0; i < numPointsArg.getValue(); ++i)
            weights[i] /= sumWeights;
    }

    anima::ComputeSphericalCentroid(random_sphere_points, average_point, random_sphere_points[0], weights);

    std::cout << "Initial points: " << std::endl;

    for (unsigned int i = 0; i < numPointsArg.getValue(); ++i)
        std::cout << random_sphere_points[i][0] << " " << random_sphere_points[i][1] << " " << random_sphere_points[i][2] << std::endl;

    std::cout << "Weights: ";
    for (unsigned int i = 0; i < weights.size(); ++i)
        std::cout << weights[i] << " ";
    std::cout << std::endl;

    std::cout << "Average point: " << average_point[0] << " " << average_point[1] << " " << average_point[2] << std::endl;

    return EXIT_SUCCESS;
}
