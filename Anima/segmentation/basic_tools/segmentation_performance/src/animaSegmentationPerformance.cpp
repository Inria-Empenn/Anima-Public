/**
 * @file main.cpp
 * @brief Contains the main function
 * @author Florent  LERAY
 * @author Baptiste LAURENT
 * @date 13/04/2016
 * @version 2.1
 */

#include "SegPerfApp.h"

/**
* @brief This function detect and applies the "about" command line switch before TCLAP.
* @param    [in] argc is the number of command line arguments.
* @param    [in] argv is the table  of command line arguments.
* @return true if "About switch" is find on command line
*/
bool detectAboutSwitch(int argc, char * *argv)
{
   bool bFound = false;
   for (int i = 0; (i < argc && !bFound); ++i)
   {
      if (argv[i] && argv[i][0])
      {
         if (!strcmp("--About", argv[i]))
         {
         }
         else if (argv[i][0] == '-' && argv[i][1] != '-')
         {
            int iLen = strlen(argv[i]);
            for (int j = 1; (j < iLen && !bFound); ++j)
            {
               bFound = argv[i][j] == 'A';
            }
         }
      }
   }
   if (bFound)
   {
      CSegPerfApp::about();
   }

   return bFound;
}

/**
 * @fn int main (int argc, char * *argv)
 * @brief Application entry point.
 * @param    [in] argc is the number of command line arguments.
 * @param    [in] argv is the table  of command line arguments.
 * @return positive value in success or -1 if field to write all measures .
 */
int main(int argc, char * *argv)
{
   int iRes = -1;
   if (!(detectAboutSwitch(argc, argv) && argc == 2))
   {
      CSegPerfApp oSegPerfApp;
      if(oSegPerfApp.init(argc, argv))
      {
         if (oSegPerfApp.checkParamsChoerancy())
         {
            oSegPerfApp.checkOutputChoerancy();
            oSegPerfApp.preparOutput();
            iRes = (int)oSegPerfApp.play();
         }
      }
   }

   return iRes;
}