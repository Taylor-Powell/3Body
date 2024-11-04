#ifndef COUPLED22_H
#define COUPLED22_H

#include <string>
#include <vector>

namespace coupled22 {
    class Data {
        public:
            Data() {}
            Data(std::string infile_name) {
                infile = infile_name;
                Lmin = Lmax = dL = Emax = m[0] = m[1] = -1.0;
                at = g[0] = g[1] = -1.0;
                nDivs = 10;
                for (int i = 0; i < 3; i++) nP[i] = -1.0;
                readData(infile);
                setOutfile();
            }
            ~Data() {}

            void readData(std::string filename);
            void setOutfile();
            void printParams();
            void freeSpecOut();
            void intSpecOut();

        private:
            int nDivs, nShell, nStates, nP[3];
            double at, m[2], mFlatte, g[2];
            double alpha, dL, Lmin, Lmax, Emax;
            std::string infile, outfile, nP_str;
            std::vector<double> Lvals;
    };
}

#endif