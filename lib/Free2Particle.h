#ifndef FREE2PARTICLE_H_INCLUDED
#define FREE2PARTICLE_H_INCLUDED

#include <vector>
#include <string>

namespace FreeSpec
{
    class Data {
        public:
            // Constructors and Destructor
            Data() {}
            Data(std::string infile_name) { 
                infile = infile_name;
                Lmin = Lmax = dL = Emax = msq[0] = -1.0;
                readData(infile_name); 
                setOutfile();
            }
            Data(std::string infile_name, std::string outfile_name) { 
                infile = infile_name;
                outfile = outfile_name;
                Lmin = Lmax = dL = Emax = msq[0] = -1.0;
                readData(infile);
                setOutfile();
            }
            ~Data() {}

            // Inline functions

            // Functions defined in read_Ecm_dat.cpp
            void readData(std::string filename);
            void setOutfile();
            void printParams();
            void free_2_out();

        private:
            int nP[3];
            double Lmax, Lmin, dL, Emax, msq[3];
            std::vector<double> Lvals;
            std::string infile, outfile, nP_str;
    };
}

#endif