#include "ampl/ampl.h"
#include <iostream>
using namespace ampl;

int main() {
  try {
    std::string ampl_path = AMPLAPI_DIR;
    std::string ampl_files_path = MPC_AMPL_PATH;

    ampl::Environment env(ampl_path);
    AMPL ampl(env);

    // Load model file
    try {
      ampl.read(ampl_files_path + "/examples/example1_files/main.mod");
      ampl_files_path += "/examples/example1_files/";
    } catch (const std::exception &e) {
      ampl_files_path += "/example1_files/";
      ampl.read(ampl_files_path + "/main.mod");
    }

    ampl.read(ampl_files_path + "/cost.mod");
    ampl.read(ampl_files_path + "/eq_con.mod");
    ampl.read(ampl_files_path + "/ineq_con.mod");
    ampl.readData(ampl_files_path + "/Data.dat");

    // Solve with default solver (set AMPL's solver separately)
    ampl.setOption("solver", "ipopt");
    ampl.eval("drop restriction2;"); // to make the problem feasible
    ampl.eval("drop restriction3;"); // to make the problem feasible

    switch (1) {
    case 0:
      ampl.eval("restore restriction2;");
      ampl.eval("solve;");
      break;
    case 1:
      ampl.eval("restore restriction3;");
      ampl.solve();
      break;
    }

    // Get variables
    Variable x = ampl.getVariable("x");

    std::cout << "x = " << x.value() << std::endl;
    // Get objective
  } catch (const std::exception &e) {
    std::cout << "Error: " << e.what() << std::endl;
  }
}
