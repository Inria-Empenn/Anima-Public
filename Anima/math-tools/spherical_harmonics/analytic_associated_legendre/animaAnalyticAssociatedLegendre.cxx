#include <tclap/CmdLine.h>
#include <iostream>
#include <vector>
#include <cmath>

unsigned int gcd(unsigned int a, unsigned int b)
{
    if (a < b)
    {
        unsigned int tmp = a;
        a = b;
        b = tmp;
    }

    while ((a != b) && (b != 0))
    {
        unsigned int newB = a % b;
        a = b;
        b = newB;
    }

    return a;
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<unsigned int> orderArg("o","order","Polynomial maximal order",false,4,"polynomial order",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    unsigned int orderAsked = orderArg.getValue();

    std::cout << "Order 0: $P_0^0(x) = 1$" << std::endl;
    if (orderAsked > 0)
        std::cout << "Order 1: $P_1^0(x) = x$" << std::endl;

    if (orderAsked < 2)
        return EXIT_SUCCESS;

    unsigned int maxIterNumber = orderAsked / 2;

    unsigned int currentAlpha = 0, currentGamma = 0;
    std::vector <int> currentBeta(1,1), currentLambda(1,1);

    for (unsigned int i = 0;i < maxIterNumber;++i)
    {
        unsigned int newIndex = i + 1;
        unsigned int newSize = newIndex + 1;
        std::vector <int> newBetaNumerators(newSize,1), newLambdaNumerators(newSize,1);

        unsigned int pow2Value = std::pow(2,currentGamma - currentAlpha);

        // Deal with P_{2*newIndex}^0
        unsigned int newAlpha = currentGamma + 1;
        newBetaNumerators[0] = - (2 * i + 1) * pow2Value * currentBeta[0];
        std::vector <unsigned int> newBetaDenominators(newSize, i + 1);
        std::vector <unsigned int> newLambdaDenominators(newSize, (i + 1) * (2 * i + 3));

        for (unsigned int j = 1;j < newIndex;++j)
            newBetaNumerators[j] = (4 * i + 3) * currentLambda[j - 1] - (2 * i + 1) * pow2Value * currentBeta[j];

        newBetaNumerators[newIndex] = (4 * i + 3) * currentLambda[i];

        unsigned int maxPow2Index = 0;
        for (unsigned int j = 0;j < newSize;++j)
        {
            unsigned int test = std::abs(newBetaNumerators[j]);
            int gcdValue = gcd(test,newBetaDenominators[j]);

            if (gcdValue > 1)
            {
                newBetaNumerators[j] /= gcdValue;
                newBetaDenominators[j] /= gcdValue;
            }

            unsigned int powIndex = 0;
            unsigned int testedDenominator = newBetaDenominators[j];
            while (testedDenominator % 2 == 0)
            {
                ++powIndex;
                testedDenominator /= 2;
            }

            if (powIndex > maxPow2Index)
                maxPow2Index = powIndex;
        }

        newAlpha += maxPow2Index;
        for (unsigned int j = 0;j < newSize;++j)
        {
            newBetaNumerators[j] *= std::pow(2,maxPow2Index);
            if (newBetaNumerators[j] % newBetaDenominators[j] != 0)
            {
                std::cerr << "Problem beta " << newBetaNumerators[j] << " " << newBetaDenominators[j] << std::endl;
                return EXIT_FAILURE;
            }

            newBetaNumerators[j] /= newBetaDenominators[j];
        }

        // Now deal with P_{2*newIndex + 1}^0
        unsigned int newGamma = currentGamma + 1;

        newLambdaNumerators[0] = - (4 * i + 5) * (2 * i + 1) * pow2Value * currentBeta[0] - 4 * (i + 1) * (i + 1) * currentLambda[0];

        for (unsigned int j = 1;j < newIndex;++j)
            newLambdaNumerators[j] = (4 * i + 5) * ((4 * i + 3) * currentLambda[j - 1] - (2 * i + 1) * pow2Value * currentBeta[j]) - 4 * (i + 1) * (i + 1) * currentLambda[j];

        newLambdaNumerators[newIndex] = (4 * i + 5) * (4 * i + 3) * currentLambda[i];

        maxPow2Index = 0;
        for (unsigned int j = 0;j < newSize;++j)
        {
            unsigned int test = std::abs(newLambdaNumerators[j]);
            int gcdValue = gcd(test,newLambdaDenominators[j]);

            if (gcdValue > 1)
            {
                newLambdaNumerators[j] /= gcdValue;
                newLambdaDenominators[j] /= gcdValue;
            }

            unsigned int powIndex = 0;
            unsigned int testedDenominator = newLambdaDenominators[j];
            while (testedDenominator % 2 == 0)
            {
                ++powIndex;
                testedDenominator /= 2;
            }

            if (powIndex > maxPow2Index)
                maxPow2Index = powIndex;
        }

        newGamma += maxPow2Index;
        for (unsigned int j = 0;j < newSize;++j)
        {
            newLambdaNumerators[j] *= std::pow(2,maxPow2Index);
            if (newLambdaNumerators[j] % newLambdaDenominators[j] != 0)
            {
                std::cerr << "Problem lambda " << j << " " << newLambdaNumerators[j] << " " << newLambdaDenominators[j] << std::endl;
                return EXIT_FAILURE;
            }

            newLambdaNumerators[j] /= newLambdaDenominators[j];
        }

        currentAlpha = newAlpha;
        currentGamma = newGamma;
        currentBeta = newBetaNumerators;
        currentLambda = newLambdaNumerators;

        // Print polynomials
        std::cout << "Order " << 2 * newIndex << ": $P_" << 2 * newIndex << "^0(x) = \\frac{1}{" << std::pow(2,currentAlpha) << "} \\left(";
        for (int j = newIndex;j >= 0;--j)
        {
            unsigned int power = 2 * j;
            if (currentBeta[j] > 0)
            {
                if (j != newIndex)
                    std::cout << " + ";
                else
                    std::cout << " ";

                std::cout << currentBeta[j];
            }
            else
                std::cout << " - " << std::abs(currentBeta[j]);

            if (j > 0)
                std::cout << " x^" << power;
        }
        std::cout << " \\right)$" << std::endl;

        if (orderAsked > 2 * newIndex)
        {
            std::cout << "Order " << 2 * newIndex + 1 << ": $P_" << 2 * newIndex + 1 << "^0(x) = \\frac{1}{" << std::pow(2,currentGamma) << "} \\left(";
            for (int j = newIndex;j >= 0;--j)
            {
                unsigned int power = 2 * j + 1;
                if (currentLambda[j] > 0)
                {
                    if (j != newIndex)
                        std::cout << " + ";
                    else
                        std::cout << " ";

                    std::cout << currentLambda[j];
                }
                else
                    std::cout << " - " << std::abs(currentLambda[j]);

                if (j > 0)
                    std::cout << " x^" << power;
                else
                    std::cout << " x";
            }
            std::cout << " \\right)$" << std::endl;
        }
    }

    return EXIT_SUCCESS;
}
