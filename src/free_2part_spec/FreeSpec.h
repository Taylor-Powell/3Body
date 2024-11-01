#include <fstream>
#include <string>

namespace FreeSpec
{
    class Data {
        public:
            // Constructors and Destructor
            Data() {}
            Data(std::string infile_name) { 
                infile = infile_name;
                readData(infile_name); 
            }
            Data(std::string infile_name, std::string outfile_name) { 
                infile = infile_name;
                outfile = outfile_name;
                readData(infile);
            }
            ~Data() {}

            // Inline functions
            void setOutputFile(std::string filename) { outfile = filename; }

            // Functions defined in read_Ecm_dat.cpp
            void readData(std::string filename);
            void printParams();
            void free_2_out();
            double free_spec_2part(double L, int n[]);            
            std::vector<double> gen_free_spec_2part(double L, char flag = 'A');



        private:
            int nP[3];
            double Lmax, Lmin, dL, Emax, msq[3];
            std::string infile, outfile, nP_str;
    };
}